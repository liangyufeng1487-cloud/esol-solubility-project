# ESOL 水溶解度预测课程项目

## 最快运行方法 Google Colab
1. 打开 https://colab.research.google.com/ ，选择“文件 → 上传笔记本”。
2. 上传 ESOL_Quick_Run.ipynb。运行时保持 CPU，不需要 GPU。
3. 选择“运行时 → 全部运行”。首次运行会安装依赖，安装耗时取决于网络，可能需要数分钟。
4. 看到 Completed in ... seconds 即训练完成。本次环境中计算、调参和出图共 9.13 秒，安装时间未计入，其他机器会不同。
5. 查看下方结果表与图片；最后一个单元格会下载你自己运行的结果 ZIP。

## 本地命令行方法
需要 Python 3.11 或 3.12；在本项目文件夹打开终端：
```bash
python -m pip install -r requirements.txt
python run_project.py
```
出现 rdkit 安装失败时，优先使用上面的 Colab。项目已包含原始 CSV，通常不需要下载数据。

## 项目在研究什么
根据分子结构 SMILES 预测 log10 水溶解度，属于药物发现中的分子性质预测。
这是回归任务，报告 RMSE、MAE、R²；RMSE 越低越好，不能把 R² 写成分类准确率。
原始 1128 条数据，按 canonical SMILES 合并 11 条重复后剩 1117 个结构。

## 实验设计和优化参数
- Mean：均值基线；Ridge4：4 描述符；Ridge10：10 描述符；RF_default：默认结构参数的森林；RF_tuned：验证集选参数的森林。
- 基础 4 特征：分子量、LogP、可旋转键、芳香重原子比例。扩展加入 HBD、HBA、TPSA、环数、重原子数、sp3 碳比例。
- 森林固定 200 棵树。搜索 max_depth=None/8，min_samples_leaf=1/2/4，max_features=1.0/0.7，共 12 组。
- 两种划分均选择 None、1、0.7。0.7 表示每个节点考虑约 70% 的特征。
- 先在训练集拟合、验证集选择；最后在训练加验证集重新拟合，独立测试集评分。
- 主实验按分子骨架分组；次实验随机划分。骨架实验训练/验证/测试为 781/167/169。

## 已实际运行的关键结果
骨架测试：Ridge4 RMSE=1.4148；RF_default=0.9340；RF_tuned=0.9248，MAE=0.7150，R²=0.8162。
随机测试：RF_tuned RMSE=0.7516。
随机森林优于线性基线；调参在骨架测试上的提升只有约 1%，不应夸大成显著优化。
骨架调参差值的 95% 配对 bootstrap 区间包含零。详见 results/paired_bootstrap.csv。
不同硬件或库版本可能产生细小差异；报告数值来自包内 run_metadata.json 记录的那次运行。

## 结果文件
results/metrics.csv：全部模型测试指标。
results/tuning_trials.csv：每一组参数的验证 RMSE。
results/test_predictions.csv：逐个分子的测试预测。
results/scaffold_split.csv 和 random_split.csv：固定划分。
results/run_metadata.json：软件版本、参数、数据校验值和运行时间。
results/*.png：报告和演示使用的真实计算图。
submission/：英文报告 Word/PDF 草稿、英文演示 PDF/PPTX 和中文讲解。

## 提交前必须完成
1. 自己运行一次笔记本并确认理解结果；若改参数导致结果变化，同时更新报告与演示。
2. 打开 submission/ESOL_Report_Draft.docx，将 [STUDENT ID] 改为学号，核对英文姓名。
3. 在 GitHub 创建可访问的仓库（例如 esol-solubility-project）。进入仓库 Add file → Upload files，上传 run_project.py、requirements.txt、README_CN.md、REFERENCES.md、笔记本，以及 data 和 results 文件夹。不要只上传尚未解压的 ZIP。
4. 复制实际仓库主页网址，替换报告目录页的 [GITHUB REPOSITORY URL]。用 Word“另存为/导出 PDF”。用未登录窗口确认老师能访问仓库。
5. 核对演示封面的姓名；按学校要求提交报告 PDF 和演示 PDF。GitHub 仓库尚未由本项目包创建。
6. 按截图，截止时间为 2026 年 11 月 16 日 14:30（以课程页面为准），汇报 10 分钟加 5 分钟问答。

## 使用边界
这是小型公开数据集上的方法比较，不是新药发现，也没有前瞻实验验证。
只使用一次随机划分和一次自定义骨架划分，不宣称达到官方基准最优成绩。
bootstrap 区间仅描述固定测试集和模型条件下的分子重采样，不包含换种子或重新划分的变动。
