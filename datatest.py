import numpy as np
import pandas as pd
import keras
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense
from sklearn.metrics import classification_report
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, mutual_info_classif,chi2,f_classif
from sklearn.metrics import confusion_matrix
from sklearn.cluster import KMeans
import tensorflow as tf
#加入忽略
import os
from tensorflow.python.keras import Sequential

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'


def feature_select(x_data, tar, fea_num):
    selector = SelectKBest(score_func=chi2, k=fea_num)
    selector.fit(x_data, tar)
    selected_index = selector.get_support(True)
    print("selected index:", selected_index)
    part_data = selector.transform(x_data)
    return part_data, selected_index


def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range


def output_csv(list,csv_name):
    #mat = np.matrix(list)
    #tmat = mat.T
    test = pd.DataFrame(data=list)
    test.to_csv(csv_name, encoding='gbk', index=0, header=0)




def draw_acc(acc,loss,epochs):
    plt.plot(epochs, acc, color='r', label='acc')  # r表示红色
    plt.plot(epochs, loss, color='b', label='loss')  # 也可以用RGB值表示颜色
    #####非必须内容#########
    plt.xlabel('epochs')  # x轴表示
    plt.ylabel('value')  # y轴表示
    plt.title("loss$acc")  # 图标标题表示
    plt.legend()  # 每条折线的label显示
    plt.show()


# def local_test_sample(testdir, all_sample):
#     test_index_list = []
#     for fold in range(1, 3):
#         test_index =[]
#         test_file_name = testdir + str(fold) + '.csv'
#         test_data = pd.read_csv(test_file_name, index_col=0)
#         test_sam_name = test_data.index
#         for sam in test_sam_name:
#             tmp=np.where(all_sample == sam)
#             test_index.extend(tmp)
#         test_index_list.append(test_index)
#     return test_index_list


dataset = pd.read_csv('STAD_latter/STAD_all.csv', index_col=0)
gene_name = dataset.columns
fea_number=100
#all_sample = dataset.index
#testdir="BRCAsubtype/testset1/testdata/test_data"
#test_sample_list = local_test_sample(testdir, all_sample)
#print(len(test_sample_list))
print(dataset.head(5))
X_OR = dataset.iloc[:, :-1].values

# x_nor = normalization(X_OR)
y = dataset.iloc[:, -1].values

x, fea_index = feature_select(X_OR, y, fea_number)
feature_gene = gene_name[fea_index, ]
print(feature_gene)
output_csv(feature_gene, csv_name="C:/Users/孙佩硕的magic/Desktop/STAD_latter/Chi-square_test/feature"+str(fea_number)+".csv")
X = normalization(x)

features = len(fea_index)
print("features_num:", features)
print("X.shape:", X.shape)
print("y.shape:", y.shape)
encoder = LabelEncoder()
y1 = encoder.fit_transform(y)
Y = pd.get_dummies(y).values
print(Y.shape)
splits = StratifiedShuffleSplit(n_splits=10, test_size=0.1, random_state=5)
pre_class = []
true_class = []
# all_index = range(0, len(all_sample))

#for test_index in test_sample_list:
    # test_index = np.array(test_index, dtype=int)
    # test_index = test_index.tolist()
    # print(type(test_index))
    # print(test_index)
    # X2 = X.copy()
    # print(type(X2))
    # X_test = X[test_index,]
    # X_train =X2.drop(test_index, axis=0)
    # print(X_test)
    #test_index = test_sample_list
    #train_index = set(all_index)-set(test_index)
#     print("test_index:", test_index)
#     print("train_index:", train_index)
for train_index, test_index in splits.split(X, Y):
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = Y[train_index], Y[test_index]
    model = Sequential()
    model.add(Dense(50, input_shape=(features,), activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(5, activation='softmax'))
    model.compile(optimizer =tf.keras.optimizers.SGD(lr = 0.1), loss='categorical_crossentropy', metrics=['accuracy'])
    model.summary()
    early_stopping = keras.callbacks.EarlyStopping(monitor='accuracy', min_delta=0.0001,
                                                   patience=50, verbose=1, mode='max',
                                                   baseline=None, restore_best_weights=True)

    model.fit(X_train, y_train, epochs=300, callbacks=early_stopping)
    # histiry_sta=model.fit(X_train, y_train, epochs=500)
    #
    # acc = histiry_sta.history['accuracy']
    # epo = histiry_sta.epoch
    # loss = histiry_sta.history['loss']
    # draw_acc(acc, loss, epo)

    y_pred = model.predict(X_test)
    y_pred_class = np.argmax(y_pred, axis=1)
    y_test_class = np.argmax(y_test, axis=1)
    pre_class.extend(list(y_pred_class))
    true_class.extend(list(y_test_class))

print("pre", pre_class)
print("true", true_class)
con_mat = confusion_matrix(true_class, pre_class)
print(con_mat)
co = pd.DataFrame(con_mat)
co.to_csv("C:/Users/孙佩硕的magic/Desktop/STAD_latter/Chi-square_test/con_mat"+str(fea_number)+".csv")

# 其实就是记录每个数组中值最大的数的index
report = classification_report(true_class, pre_class,output_dict=True)
print(report)
df = pd.DataFrame(report).transpose()
df.to_csv("C:/Users/孙佩硕的magic/Desktop/STAD_latter/Chi-square_test/result"+str(fea_number)+".csv", index=True)
#

#
#
#
