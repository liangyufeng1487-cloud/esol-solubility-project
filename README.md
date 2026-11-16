# Molecular aqueous solubility prediction

CPU ESOL benchmark for AIDD8122. Start with `ESOL_Quick_Run.ipynb` in Google Colab, or install `requirements.txt` and run `python run_project.py`.

Full Chinese instructions: [README_CN.md](README_CN.md). Dataset and original publications: [REFERENCES.md](REFERENCES.md).

The included run uses 1117 unique structures after merging 11 duplicate rows. The primary scaffold test RMSE is 0.9248 for the tuned random forest. All validation trials, partitions and predictions are under `results/`.

Report and presentation drafts are in `submission/`. Before submission, fill in the student ID and the actual accessible GitHub repository URL in the Word report, then export a fresh PDF. Current PDF is a draft with those fields pending.

This project is an evaluation of established methods on public data. It does not claim a new model or prospective drug-discovery validation.
