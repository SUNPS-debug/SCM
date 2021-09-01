import numpy as np
import pandas as pd
import tensorflow as tf
import keras
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
from keras.models import Sequential
from keras.layers import Dense, Dropout
from sklearn.metrics import classification_report
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix
from pandas.core.frame import DataFrame
from sklearn import svm
from sklearn.metrics import accuracy_score

import os
from functools import reduce

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'


def ZscoreNormalization(x):
    """Z-score normaliaztion"""
    x = (x - np.mean(x)) / np.std(x)
    return x


def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range





# 多次测试最终产生投票结果
def voting(predict_class):
    vote_input = np.matrix(predict_class)
    print(vote_input.shape)
    ncol = vote_input.shape[1]
    vote_result = []
    for col in range(0, ncol):
        vote_all = vote_input[:, col]
        vote_list = vote_all.tolist()
        v = max(vote_list, key=vote_list.count)
        vote_result.append(v)
    print(type(vote_result))
    return vote_result


# def soft_voting(predict_class):
#     sum_pred = reduce(lambda x, y: x + y, predict_class)
#     print('sum_pred', sum_pred)
#     print('sum_pred type:', type(sum_pred))
#     soft_vote = np.argmax(sum_pred, axis=1)
#     soft_vote_array = soft_vote.A
#     soft_vote_list = list(soft_vote_array)
#     return soft_vote_list

def output_csv(list,csv_name):
    #mat = np.matrix(list)
    #tmat = mat.T
    test = pd.DataFrame(data=list)
    test.to_csv(csv_name, encoding='gbk', index=0, header=0)


# def draw_acc(acc,loss,epochs):
#     plt.plot(epochs, acc, color='r', label='acc')  # r表示红色
#     plt.plot(epochs, loss, color='b', label='loss')  # 也可以用RGB值表示颜色
#     #####非必须内容#########
#     plt.xlabel('epochs')  # x轴表示
#     plt.ylabel('value')  # y轴表示
#     plt.title("loss$acc")  # 图标标题表示
#     plt.legend()  # 每条折线的label显示
#     plt.show()


def dnn_classification(fold,vt):
    test_data_list = []
    train_file = 'BRCA_latter/100feature/trainset/train_data'+str(fold)+'.csv'
    train_data = pd.read_csv(train_file, index_col=0)
    #train_data = pd.read_csv(train_file)
    print('train_data:', train_data.head(5))
    for i in range(1, (vt+1)):
        test_file_name = 'BRCA_latter/100feature/testset'+str(i)+'/test_data' + str(fold) + '.csv'
        print(test_file_name)
        test_data = pd.read_csv(test_file_name,  index_col=0)
        #test_data = pd.read_csv(test_file_name)
        test_data_list.append(test_data)
    # # encoder Y
    encoder = LabelEncoder()
    train = train_data.iloc[:, :-1].values
    # sns.heatmap(train, vmax=0.1, vmin=0)
    # plt.show()
    train_x = normalization(train)
    y_train = train_data.iloc[:, -1].values
    print("train_shape:", train.shape)
    y1_train = encoder.fit_transform(y_train)
    train_y = pd.get_dummies(y1_train).values
    test_x = []
    test_y = []
    for te in range(0, vt):
        test_tmp = test_data_list[te]
        x_test = DataFrame(test_tmp)
        test = x_test.iloc[:, :-1].values
        # sns.heatmap(test, vmax=0.1, vmin=0)
        # plt.show()
        test_x.append(normalization(test))
        y_test = x_test.iloc[:, -1].values
        y1_test = encoder.fit_transform(y_test)
        #test_y.append(y1_test)
        test_y.append(pd.get_dummies(y1_test).values)
    # classifier = svm.SVC(C=3, tol=0.001, kernel='rbf', gamma='auto',
    #                          decision_function_shape='ovr', random_state=2)  # ovr:一对多策略
    # classifier.fit(train_x, y1_train)
    # tra_label = classifier.predict(train_x)
    # print("训练集：", accuracy_score(y_train, tra_label))
    features = train_x.shape[1]
    model = Sequential()
    model.add(Dense(50, input_shape=(features,), activation='relu'))
    model.add(Dense(10, activation='relu'))
    model.add(Dense(4, activation='softmax'))
    model.compile(optimizer ='sgd', loss='categorical_crossentropy', metrics=['accuracy'])
    model.summary()
    early_stopping = keras.callbacks.EarlyStopping(monitor='accuracy', min_delta=0.001,
                                                   patience=10, verbose=1, mode='max',
                                                   baseline=None, restore_best_weights=False)

    model.fit(train_x, train_y, epochs=300, callbacks=early_stopping)
    # acc = histiry_sta.history['accuracy']
    # epo = histiry_sta.epoch
    # loss = histiry_sta.history['loss']
    # #draw_acc(acc, loss, epo)
    y_test_class = []
    y_pred_class = []
    for test_t in range(0, vt):
        # y_pred=classifier.predict(test_x[test_t])
        # y_pred_class.append(y_pred)
        # print(y_pred_class)
        # y_test_class.append(test_y[test_t])
        y_pred = model.predict(test_x[test_t])
        #y_pred_mat = np.matrix(y_pred)
        #y_pred_class.append(y_pred_mat)
        y_pred_class.append(np.argmax(y_pred, axis=1))
        y_test_class.append(np.argmax(test_y[test_t], axis=1))
    return y_test_class, y_pred_class

#y_test_class, y_pred_class = dnn_classification(1)
# vote_true = voting(y_test_class)
# vote_pred = voting(y_pred_class)
# print('vote_true:', vote_true)
# print('vote_pre:', vote_pred)
pre_class = []
true_class = []

for f in range(1, 11):
    y_test_class, y_pred_class = dnn_classification(f, vt=5)
    vote_true = voting(y_test_class)
    vote_pred = voting(y_pred_class)
    print('vote_true:', vote_true)
    print('vote_pre:', vote_pred)
    pre_class.extend(list(vote_pred))
    true_class.extend(list(vote_true))
print('predict:', pre_class)
print('true:', true_class)
#y_test_mat = np.matrix(y_pred_class)
#test_mat = np.append(test_mat, y_test_mat, axis=1)
#     vote_pred = voting(y_pred_class)
#     vote_true = voting(y_test_class)
#     pre_class.extend(list(vote_pred))
#     true_class.extend(list(vote_true))
# output_csv(true_class, 'true_type.csv')
# output_csv(pre_class, 'vote_type.csv')
# output_csv(test_mat.T, 'predict_5times.csv')
# true_class, pre_class= vote_true, vote_pred

fea_number=100
con_mat = confusion_matrix(true_class, pre_class)
print(con_mat)
co = pd.DataFrame(con_mat)
co.to_csv("C:/Users/孙佩硕的magic/Desktop/BRCA_latter/new/con_mat"+str(fea_number)+".csv")

report = classification_report(true_class, pre_class, output_dict=True)
print(report)
df = pd.DataFrame(report).transpose()
df.to_csv("C:/Users/孙佩硕的magic/Desktop/BRCA_latter/new/result"+str(fea_number)+".csv", index=True)

#



