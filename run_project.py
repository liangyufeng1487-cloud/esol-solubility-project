"""CPU-only ESOL benchmark. Run: python run_project.py"""
from pathlib import Path
import hashlib, itertools, json, platform, time, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn import __version__ as sklearn_version
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.base import clone

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)
SEED = 42
FEATURES = ['MolWt','LogP','RotatableBonds','AromaticProportion','HBD','HBA','TPSA','RingCount','HeavyAtoms','FractionCSP3']
TARGET = 'measured log solubility in mols per litre'

def features(m):
    heavy = m.GetNumHeavyAtoms()
    return [Descriptors.MolWt(m), Crippen.MolLogP(m), Lipinski.NumRotatableBonds(m),
            sum(a.GetIsAromatic() for a in m.GetAtoms()) / max(heavy, 1),
            Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m), rdMolDescriptors.CalcTPSA(m),
            rdMolDescriptors.CalcNumRings(m), heavy, rdMolDescriptors.CalcFractionCSP3(m)]

def metrics(y, p):
    return dict(RMSE=float(np.sqrt(mean_squared_error(y,p))), MAE=float(mean_absolute_error(y,p)), R2=float(r2_score(y,p)))

def split_indices(frame, mode):
    ids = np.arange(len(frame))
    if mode == 'random':
        tr, rest = train_test_split(ids, test_size=.30, random_state=SEED)
        va, te = train_test_split(rest, test_size=.5, random_state=SEED)
    else:
        # Greedy, deterministic 70/15/15 scaffold grouping; all acyclic molecules
        # share the empty scaffold. Fractions are approximate, not exact.
        groups = [list(v) for v in frame.groupby('scaffold',sort=True).indices.values()]
        groups.sort(key=lambda g:(-len(g),min(g)))
        tr, va, te = [], [], []
        for group in groups:
            if len(tr)+len(group) <= .70*len(frame): tr.extend(group)
            elif len(va)+len(group) <= .15*len(frame): va.extend(group)
            else: te.extend(group)
        tr, va, te = map(np.array,(tr,va,te))
        sets=[set(frame.iloc[idx].scaffold) for idx in (tr,va,te)]
        assert not (sets[0]&sets[1] or sets[0]&sets[2] or sets[1]&sets[2])
    assert len(set(tr)&set(va)) == len(set(tr)&set(te)) == len(set(va)&set(te)) == 0
    assert len(tr)+len(va)+len(te) == len(frame)
    return tr, va, te

def main():
    start = time.perf_counter()
    source = ROOT/'data/delaney-processed.csv'
    if not source.exists():
        import urllib.request
        source.parent.mkdir(exist_ok=True)
        urllib.request.urlretrieve('https://raw.githubusercontent.com/deepchem/deepchem/master/datasets/delaney-processed.csv',source)
    raw = pd.read_csv(source)
    records=[]
    invalid=0
    for _, row in raw.iterrows():
        mol = Chem.MolFromSmiles(str(row['smiles']).strip())
        if mol is None or not np.isfinite(row[TARGET]):
            invalid += 1
            continue
        records.append({'canonical_smiles':Chem.MolToSmiles(mol), 'logS':float(row[TARGET])})
    data = pd.DataFrame(records).groupby('canonical_smiles',sort=True,as_index=False).agg(logS=('logS','mean'),measurements=('logS','size'))
    mols = [Chem.MolFromSmiles(s) for s in data.canonical_smiles]
    data['scaffold'] = [MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False) for m in mols]
    X = np.array([features(m) for m in mols]); y = data.logS.to_numpy()
    assert np.isfinite(X).all()
    # CSV-provided predicted ESOL values and descriptors are NEVER used as features.
    pd.concat([data,pd.DataFrame(X,columns=FEATURES)],axis=1).to_csv(OUT/'clean_features.csv',index=False)
    metadata = {'raw_rows':len(raw),'invalid_rows':invalid,'unique_molecules':len(data),
                'duplicate_rows_merged':len(records)-len(data),'seed':SEED,
                'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                'versions':{'python':platform.python_version(),'rdkit':rdBase.rdkitVersion,'sklearn':sklearn_version,'numpy':np.__version__,'pandas':pd.__version__},'splits':{}}
    all_scores=[]; trials=[]; predictions=[]
    for mode in ['scaffold','random']:
        tr,va,te=split_indices(data,mode)
        assign=data[['canonical_smiles','scaffold','logS']].copy()
        assign['split']=''; assign.loc[tr,'split']='train'; assign.loc[va,'split']='validation'; assign.loc[te,'split']='test'
        assign.to_csv(OUT/f'{mode}_split.csv',index=False)
        best_score=float('inf'); best_params=None
        # All 12 candidate models see TRAIN only; validation selects parameters.
        for depth,leaf,mf in itertools.product([None,8],[1,2,4],[1.0,.7]):
            params=dict(n_estimators=200,max_depth=depth,min_samples_leaf=leaf,max_features=mf,random_state=SEED,n_jobs=2)
            model=RandomForestRegressor(**params).fit(X[tr],y[tr])
            score=metrics(y[va],model.predict(X[va]))['RMSE']
            trials.append({'split':mode,'max_depth':str(depth),'min_samples_leaf':leaf,'max_features':mf,'validation_RMSE':score})
            if score < best_score: best_score,best_params=score,params
        metadata['splits'][mode]={'train':len(tr),'validation':len(va),'test':len(te),
            'best_params':best_params,'best_validation_RMSE':best_score,
            'scaffold_overlap_train_test':len(set(data.iloc[tr].scaffold)&set(data.iloc[te].scaffold))}
        configs=[('Mean',DummyRegressor(),10),('Ridge4',make_pipeline(StandardScaler(),Ridge(alpha=1.0)),4),
                 ('Ridge10',make_pipeline(StandardScaler(),Ridge(alpha=1.0)),10),
                 ('RF_default',RandomForestRegressor(n_estimators=200,random_state=SEED,n_jobs=2),10),
                 ('RF_tuned',RandomForestRegressor(**best_params),10)]
        fig,ax=plt.subplots(figsize=(5.8,5.2))
        for name,model,k in configs:
            # No test data are used for fitting, scaling, or parameter selection.
            dev=np.r_[tr,va]; model.fit(X[dev,:k],y[dev]); pred=model.predict(X[te,:k])
            score=metrics(y[te],pred)
            all_scores.append({'split':mode,'model':name,**score})
            for ix,p in zip(te,pred): predictions.append({'split':mode,'model':name,'canonical_smiles':data.iloc[ix].canonical_smiles,'actual':y[ix],'predicted':p})
            print(mode,name,{k:round(v,4) for k,v in score.items()},flush=True)
            if name=='RF_tuned':
                ax.scatter(y[te],pred,s=18,alpha=.65,color='#176B78')
                pd.DataFrame({'feature':FEATURES,'impurity_importance':model.feature_importances_}).sort_values('impurity_importance',ascending=False).to_csv(OUT/f'{mode}_feature_importance.csv',index=False)
        ax.plot([-12,2],[-12,2],'--',color='#555555',linewidth=1)
        ax.set(xlabel='Measured logS (mol/L)',ylabel='Predicted logS (mol/L)',title=f'Tuned random forest: {mode} test')
        fig.tight_layout();fig.savefig(OUT/f'{mode}_parity.png',dpi=190);plt.close(fig)
    scores=pd.DataFrame(all_scores);scores.to_csv(OUT/'metrics.csv',index=False)
    pd.DataFrame(trials).to_csv(OUT/'tuning_trials.csv',index=False)
    preds=pd.DataFrame(predictions);preds.to_csv(OUT/'test_predictions.csv',index=False)
    # Paired bootstrap of test molecules; conditional on each fixed split/model.
    cis=[]
    for mode in ['scaffold','random']:
        a=preds[(preds.split==mode)&(preds.model=='RF_default')].reset_index(drop=True)
        b=preds[(preds.split==mode)&(preds.model=='RF_tuned')].reset_index(drop=True)
        assert a.canonical_smiles.equals(b.canonical_smiles)
        rng=np.random.default_rng(SEED); values=[]
        for _ in range(1000):
            idx=rng.integers(0,len(a),len(a))
            values.append(np.sqrt(np.mean((b.actual.to_numpy()[idx]-b.predicted.to_numpy()[idx])**2))-np.sqrt(np.mean((a.actual.to_numpy()[idx]-a.predicted.to_numpy()[idx])**2)))
        cis.append({'split':mode,'delta_RMSE_tuned_minus_default_low':float(np.quantile(values,.025)),'delta_RMSE_tuned_minus_default_high':float(np.quantile(values,.975))})
    pd.DataFrame(cis).to_csv(OUT/'paired_bootstrap.csv',index=False)
    fig,ax=plt.subplots(figsize=(8.6,4.6))
    names=scores.model.unique();xx=np.arange(len(names))
    for offset,mode,color in [(-.18,'scaffold','#176B78'),(.18,'random','#D98B35')]:
        vals=scores[scores.split==mode].set_index('model').loc[names,'RMSE']
        bars=ax.bar(xx+offset,vals,.34,label=mode,color=color)
        ax.bar_label(bars,fmt='%.2f',fontsize=8,padding=3)
    ax.set_xticks(xx,names);ax.set(ylabel='Test RMSE (log mol/L)',title='Model comparison on two fixed partitions')
    ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(OUT/'model_comparison.png',dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,3.5));ax.hist(y,bins=30,color='#176B78',edgecolor='white')
    ax.set(xlabel='Measured logS (mol/L)',ylabel='Molecules',title='ESOL target distribution after deduplication')
    fig.tight_layout();fig.savefig(OUT/'target_distribution.png',dpi=190);plt.close(fig)
    metadata['runtime_seconds']=round(time.perf_counter()-start,2)
    (OUT/'run_metadata.json').write_text(json.dumps(metadata,indent=2))
    print('Completed in',metadata['runtime_seconds'],'seconds. Results:',OUT,flush=True)

if __name__=='__main__': main()
