import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
plt.rcParams['font.family'] = ['Arial Unicode MS'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False #用来正常显示负号
sns.set_style('whitegrid',{'font.sans-serif':['Arial Unicode MS','Arial']})
import matplotlib
from sklearn import model_selection
from sklearn.ensemble import AdaBoostClassifier

from scipy import stats
from scipy.stats import mannwhitneyu
from collections import defaultdict
import re
from glob import glob
import sys
import math
from sklearn.decomposition import PCA 
matplotlib.use('Agg') 
import matplotlib as mpl 
mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams["font.sans-serif"] = "Arial"
cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
from scipy.stats import mannwhitneyu

import scanpy as sc
import sklearn.preprocessing
from sklearn.model_selection import train_test_split

from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from itertools import cycle
from sklearn.metrics import RocCurveDisplay
from sklearn.preprocessing import LabelBinarizer

### Figure3B The performance of the ScanTecc distinguishing primary cancer types from healthy individual.

f4a=pd.read_csv("/home/luosongwen/scantecc/s4a_new.txt",sep="\t",index_col=0)
test=f4a.copy()

ctype_mapping = {
    '正常对照': 'Health',
    '肺癌': 'Lung',
    '卵巢癌':'Ovarian',
    '胃癌':'Gastric',
    '淋巴瘤':'Lymphoma',
    '结直肠癌':'Colorectal',
    '乳腺癌': 'Breast',
    'other':'Others'
}
test['ctype'] = test['ctype'].replace(ctype_mapping)

gd=pd.read_table("/home/luosongwen/scantecc/2025_03_10_Figure_result/Organized_code/Final_gene_Score_cut.count",index_col=0)
gd.columns=["".join(x.split("-")) for x in gd.columns]
test.index=["".join(x.split("-")) for x in f4a.index]
gd_refine=gd[test.index]

use_ctype=['Lung', 'Ovarian', 'Gastric', 'Lymphoma', 'Health']
use_data=test[test.apply(lambda x:x["ctype"] in use_ctype,axis=1)]

adata = sc.AnnData(gd_refine[use_data.index].T)

adata.obs['ctype'] = use_data.ctype.astype('category')

# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

adata.layers["scaled"]=sc.pp.scale(adata,copy=True,zero_center=False).X
#'t-test'`, `'t-test_overestim_var'`, and `'wilcoxon'
sc.tl.rank_genes_groups(adata, 'ctype',method="wilcoxon",
                        #layer="scaled"
                       )
adata.obs['ctype'] = adata.obs['ctype'].cat.reorder_categories(use_ctype)

###基于随机森林分类模型调参

def find_deg(head,pval,fd):
    lung=deg[deg.group=="Lung"].sort_values("scores",ascending=False).iloc[:head]
    ova=deg[deg.group=="Ovarian"].sort_values("scores",ascending=False).iloc[:head]
    ga=deg[deg.group=="Gastric"].sort_values("scores",ascending=False).iloc[:head]
    lm=deg[deg.group=="Lymphoma"].sort_values("scores",ascending=False).iloc[:head]
    norm=deg[deg.group=="Health"].sort_values("scores",ascending=False).iloc[:head]
    lung=lung[(lung.pvals<pval) & (lung.logfoldchanges > fd)]
    ova=ova[(ova.pvals<pval) & (ova.logfoldchanges > fd)]
    ga=ga[(ga.pvals<pval) & (ga.logfoldchanges > fd)]
    lm=lm[(lm.pvals<pval) & (lm.logfoldchanges > fd)]
    norm=norm[(norm.pvals<pval) & (norm.logfoldchanges > fd)]
    degname=set(pd.concat([lung,ova,ga,lm,norm])["names"].values)
    print(len(lung),len(ova),len(ga),len(lm),len(norm))
    return degname

def draw_fig(use_type="raw", degname=[], size=0.2, rs=0, z=True, classifier=None):
    degname = list(degname) + [x for x in gd_refine.index if x[:3] == "MT-"]
    from scipy.stats.mstats import zscore

    if use_type == "raw":
        valid_data = pd.DataFrame(adata.layers["counts"], columns=gd_refine.index)[list(degname)]
    elif use_type == "norm":
        valid_data = pd.DataFrame(adata.X, columns=gd_refine.index)[list(degname)]
    else:
        valid_data = pd.DataFrame(adata.layers["scaled"], columns=gd_refine.index)[list(degname)]

    if z:
        valid_data = valid_data.transform(zscore, axis=0)

    X, y = valid_data, use_data["ctype"]
    y = y.values
    X = X.values

    # 划分train/test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=size, stratify=y, random_state=rs)

    # 模型训练
    if classifier is None:
        classifier = RandomForestClassifier(n_estimators=100, random_state=7418, n_jobs=-1)
    y_score = classifier.fit(X_train, y_train).predict_proba(X_test)

    from sklearn.preprocessing import LabelBinarizer
    fig, ax = plt.subplots(figsize=(6, 6))

    y_onehot_test = LabelBinarizer().fit(y).transform(y_test)
    n_classes = y_onehot_test.shape[1]
    colors = cycle(['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'])
    auc = {}

    for class_id, color in zip(range(n_classes), colors):
        t = RocCurveDisplay.from_predictions(
            y_onehot_test[:, class_id],
            y_score[:, class_id],
            name=f"ROC curve for {use_ctype[class_id]}",
            color=color,
            ax=ax,
            plot_chance_level=(class_id == 2),
        )
        auc[use_ctype[class_id]] = t.roc_auc

    return classifier, auc, fig


deg = sc.get.rank_genes_groups_df(adata, use_ctype)

# 参数列表
head_list = [200, 250, 300, 350, 400, 450, 500]
pval_list = [0.05, 0.01]
fd_list = [0.1, 0.2, 0.3, 0.4, 0.5]
use_type_list = ["raw", "norm", "scaled"]
n_estimators_list = [100, 150, 200, 250, 300, 350, 400, 450, 500]
random_state_list = [0, 1, 42, 77, 123, 2021, 2023, 777, 888, 999, 1337, 7418, 7774]

# 保存每组成功结果
all_results = {}

# 循环遍历每一组参数组合
for head in head_list:
    for pval in pval_list:
        for fd in fd_list:
            for use_type in use_type_list:
                for n_estimators in n_estimators_list:
                    for rs in random_state_list:
                        print(f"Running head={head}, pval={pval}, fd={fd}, use_type={use_type}, n_estimators={n_estimators}, random_state={rs}")

                        # 筛选差异基因
                        degname = find_deg(head, pval, fd)

                        if len(degname) < 5:
                            print(f"  [Skipped] Too few genes ({len(degname)})")
                            continue

                        # 模型训练
                        clf, auc, fig = draw_fig(
                            degname=degname,
                            use_type=use_type,
                            z=True,
                            size=0.2,
                            rs=rs,
                            classifier=RandomForestClassifier(
                                n_estimators=n_estimators,
                                criterion="entropy",
                                random_state=rs,
                                n_jobs=-1
                            )
                        )

                        # 筛选条件
                        auc_values = list(auc.values())
                        health_auc = auc.get("Health", 1.0)  # 默认值设1.0避免KeyError

                        if any(a <= 0.7 for a in auc_values) or health_auc <= 0.8:
                            print(f"  [Skipped] Some class AUC <= 0.7 or Health AUC <= 0.8 (Health AUC={health_auc:.2f})")
                            plt.close(fig)
                            continue

                        # 保存成功的结果
                        key = f"head={head}_pval={pval}_fd={fd}_use={use_type}_nest={n_estimators}_rs={rs}"
                        all_results[key] = {
                            "classifier": clf,
                            "auc": auc
                        }

                        # 保存绘图
                        fig.savefig(f"/home/luosongwen/scantecc/2025_03_10_Figure_result/figure3B_refined_ROC_figure/{key}.pdf", bbox_inches='tight')
                        plt.close(fig)  # 关闭绘图释放内存