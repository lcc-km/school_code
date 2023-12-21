import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import os,re
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import svm
import sklearn.metrics as sm
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
import graphviz
os.environ['PATH'] = os.pathsep + r'C:\Program Files\Graphviz\bin'

pwd = "C:/Users/lucc/Desktop/微创文件/9.ERA"
count_f = pwd + "/" + "GSE114157_p15_Expression_Matrix.csv"
data =  pd.read_table(count_f,sep=",")
data.index = data['sample']
data2 = data.drop(columns=['sample'])
data2 = data2[~(data2==0).all(axis=1)]
####   数据整理
data_T = pd.DataFrame(data2.values.T,columns=data2.index,index=data2.columns)
xxx = list(data_T.index.str.split("_",expand=True))
d_xxx = list(pd.DataFrame(xxx)[0])
data_T.index = d_xxx
columnNames = data_T.iloc[0]
data_T2 = data_T[1:]
data_T2.columns = columnNames

## 划分数据集
x_train,x_test,y_train,y_test=train_test_split(data_T2,data_T2.index,test_size=0.30)

#### SVM
modelSVM = svm.SVC(kernel='rbf')
modelSVM.fit(x_train, y_train)
#
print('分类准确度：{:.4f}'.format(modelSVM.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = modelSVM.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
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
dot_data = tree.export_graphviz(clf, out_file=None, feature_names=data_T2.columns, class_names=['Deiter', 'IHC','OHC'], rounded=True, filled=True) # rounded和字体有关，filled设置颜色填充
graph = graphviz.Source(dot_data)
graph.render('决策树可视化')


## 随机森林
rfc = RandomForestClassifier(random_state=0)
rfc = rfc.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(rfc.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = rfc.predict(x_test)
# 计算模型精度
bg = sm.classification_report(y_test, pred_test_y)
print('分类报告：', bg, sep='\n')
###
X = pd.concat([x_train,x_test])
Y = list(y_train) + list(y_test)
# 预测
pred_test_y = rfc.predict(X)
# 计算模型精度
bg = sm.classification_report(Y, pred_test_y)
print('分类报告：', bg, sep='\n')
######### 随机森林和决策树在10组交叉验证
rfc_l = []
clf_l = []
GBDT_l = []
modelSVM_l = []
for i in range(10):
	rfc = RandomForestClassifier(n_estimators=25)
	rfc_s = cross_val_score(rfc, X, Y, cv=10,n_jobs=10).mean()
	rfc_l.append(rfc_s)
	clf = DecisionTreeClassifier()
	clf_s = cross_val_score(clf, X, Y, cv=10,n_jobs=10).mean()
	clf_l.append(clf_s)
	GBDT = GradientBoostingClassifier()
	GBDT_s = cross_val_score(GBDT, X, Y, cv=10,n_jobs=10).mean()
	GBDT_l.append(GBDT_s)
	modelSVM = svm.SVC(kernel='rbf')
	modelSVM_s = cross_val_score(modelSVM, X, Y, cv=10,n_jobs=10).mean()
	modelSVM_l.append(modelSVM_s)

plt.plot(range(1, 11), rfc_l, label="Random Forest")
plt.plot(range(1, 11), clf_l, label="Decision Tree")
plt.plot(range(1, 11), GBDT_l, label="gradien boosting decision tree")
plt.plot(range(1, 11), modelSVM_l, label="SVM_rbf")
plt.legend()
plt.show()

## 梯度提升
GBDT=GradientBoostingClassifier()
GBDT = GBDT.fit(x_train, y_train)
print('分类准确度：{:.4f}'.format(GBDT.score(x_train, y_train)))  # 对训练集的分类准确度
# 预测
pred_test_y = GBDT.predict(X)
# 计算模型精度
bg = sm.classification_report(Y, pred_test_y)
print('分类报告：', bg, sep='\n')


GBDT_l = []
for i in range(10):
	GBDT = GradientBoostingClassifier()
	GBDT_s = cross_val_score(GBDT, X, Y, cv=10, n_jobs=10).mean()
	GBDT_l.append(GBDT_s)

plt.plot(range(1, 11), GBDT_l, label="gradien boosting decision tree")
plt.legend()
plt.show()