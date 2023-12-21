import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import os,re
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC, LinearSVC,NuSVC
import sklearn.metrics as sm
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
import graphviz
from sklearn.ensemble import VotingClassifier
from sklearn.naive_bayes import GaussianNB,BernoulliNB,MultinomialNB
## 数据模拟
from sdv.tabular import GaussianCopula
from sdv.evaluation import evaluate
import ctgan
from ctgan import CTGANSynthesizer
os.environ['PATH'] = os.pathsep + r'C:\Program Files\Graphviz\bin'

pwd = "C:/Users/lucc/Desktop/微创文件/9.ERA/LH59"
count_f = pwd + "/" + "gene_TPM.csv"
data =  pd.read_table(count_f,sep=",")
data.index = data['gene_id']
data2 = data.drop(columns=['gene_id'])

## 目标gene list
gene_f = "C:/Users/lucc/Desktop/微创文件/9.ERA/LH59/target_gene.list.txt"
gene_d = pd.read_table(gene_f,header=None)
gene_list = gene_d[0].to_list()

##   数据整理
data_T = pd.DataFrame(data2.values.T,columns=data2.index,index=data2.columns)
xxx = list(data_T.index.str.split("LH",expand=True))
d_xxx = list(pd.DataFrame(xxx)[1])
data_T.index = d_xxx
data_T2 = data_T[gene_list]
#data_T2 = data_T2[~(data_T2==0).all(axis=1)]

####  数据构造
LH2 = data_T2[data_T2.index == "plus2"]
LH7 = data_T2[data_T2.index == "plus7"]
discrete_columns = list(data_T2.columns)
#
ctgan = CTGANSynthesizer(batch_size=50,epochs=5,verbose=False)
ctgan.fit(LH2,discrete_columns)
## 生成200条数据集
samples_LH2 = ctgan.sample(200)
evaluate(samples_LH2, LH2)
#
ctgan = CTGANSynthesizer(batch_size=30,epochs=6,verbose=False)
ctgan.fit(LH7,discrete_columns)
samples_LH7 = ctgan.sample(200)
evaluate(samples_LH7, LH7)

X_samples = pd.concat([samples_LH2,samples_LH7])
Y_samples = ['plus2']*200 + ['plus7']*200
## 划分数据集
x_train,x_test,y_train,y_test=train_test_split(data_T2,data_T2.index,test_size=0.3)
## PCA降维
#pca = PCA(n_components=20)#降低的维度#是否写入
#x_train_pca = pca.fit_transform(x_train,y_train)
#x_test_pca = pca.fit_transform(x_test,y_test)

#### SVM
modelSVM = svm.SVC(kernel='rbf')
modelSVM.fit(x_train, y_train)
#
print('分类准确度：{:.4f}'.format(modelSVM.score(X, Y)))  # 对训练集的分类准确度
# 预测
pred_test_y = modelSVM.predict(X_samples)
# 计算模型精度
bg = sm.classification_report(Y_samples, pred_test_y)
print('分类报告：', bg, sep='\n')

### 决策树
clf = tree.DecisionTreeClassifier(criterion="entropy")
clf = clf.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(clf.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = clf.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
print('分类报告：', bg, sep='\n')
### plot
dot_data = tree.export_graphviz(clf, out_file=None, feature_names=data_T2.columns, class_names=['plus2','plus7'], rounded=True, filled=True) # rounded和字体有关，filled设置颜色填充
graph = graphviz.Source(dot_data)
graph.render('决策树可视化')

## 随机森林
rfc = RandomForestClassifier(n_estimators=25)
rfc = rfc.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(rfc.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = rfc.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
print('分类报告：', bg, sep='\n')

## 梯度提升
GBDT=GradientBoostingClassifier()
GBDT = GBDT.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(GBDT.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = GBDT.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
print('分类报告：', bg, sep='\n')

## 贝叶斯
#高斯分布型
gnb = GaussianNB()
gnb = gnb.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(gnb.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = gnb.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
print('分类报告：', bg, sep='\n')
#多项式型
mnb = MultinomialNB()
mnb = gnb.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(mnb.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = mnb.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
print('分类报告：', bg, sep='\n')
#伯努利型
bnb = BernoulliNB()
bnb = gnb.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(bnb.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = bnb.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
print('分类报告：', bg, sep='\n')

##### 投票
clf1 = GradientBoostingClassifier()
clf2 = RandomForestClassifier(n_estimators=25)
clf3 = tree.DecisionTreeClassifier(criterion="entropy")
eclf = VotingClassifier(estimators=[('GBDT',clf1),('rf',clf2),('tree',clf3)], voting='hard')
#使用投票法将三个模型结合在以前，estimotor采用 [(name1,clf1),(name2,clf2),...]这样的输入，和Pipeline的输入相同 voting='hard'表示硬投票
eclf.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(eclf.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = eclf.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
print('分类报告：', bg, sep='\n')

### 合并
X = pd.concat([x_train,x_test])
Y = list(y_train) + list(y_test)

######### SVM组交叉验证
LinearSVC_l = []
SVM_poly_l = []
SVM_rbf_l = []
gnb_l = []
mnb_l = []
for i in range(30):
	LinearSVC = SVC(kernel='linear')
	LinearSVC_s = cross_val_score(rfc, X, Y, cv=10,n_jobs=10).mean()
	LinearSVC_l.append(LinearSVC_s)
	SVM_poly = svm.SVC(kernel='poly')
	SVM_poly_s = cross_val_score(SVM_poly, X, Y, cv=10,n_jobs=10).mean()
	SVM_poly_l.append(SVM_poly_s)
	SVM_rbf = svm.SVC(kernel='rbf')
	SVM_rbf_s = cross_val_score(SVM_rbf, X, Y, cv=10,n_jobs=10).mean()
	SVM_rbf_l.append(SVM_rbf_s)
	gnb = GaussianNB()
	gnb_s = cross_val_score(gnb, X, Y, cv=10,n_jobs=10).mean()
	gnb_l.append(gnb_s)
	mnb = MultinomialNB()
	mnb_s = cross_val_score(mnb, X, Y, cv=10,n_jobs=10).mean()
	mnb_l.append(mnb_s)

plt.plot(range(1, 31), LinearSVC_l, label="LinearSVC")
plt.plot(range(1, 31), SVM_poly_l, label="poly")
plt.plot(range(1, 31), SVM_rbf_l, label="rbf")
plt.plot(range(1, 31), gnb_l, label="GaussianNB")
plt.plot(range(1, 31), mnb_l, label="MultinomialNB")
plt.legend()
plt.show()


