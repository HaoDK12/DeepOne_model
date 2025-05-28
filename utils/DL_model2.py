import os
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, Dropout, AveragePooling1D, Flatten, Dense, concatenate


# 定义模型结构的函数
def build_model(eLENGTH, eDEPTH, TYPE="CG"):
    inputs = []
    if 'C' in TYPE:
        input_c = Input(shape=(eLENGTH, eDEPTH), name="input_onehot")
        inputs.append(input_c)
    if 'G' in TYPE:
        input_g = Input(shape=(1,), name="input_dGB")
        inputs.append(input_g)

    for_dense = []
    if 'C' in TYPE:
        # 第一个卷积层
        conv1_out = Conv1D(100, 3, activation='relu', name="conv_3")(input_c)
        conv1_dropout_out = Dropout(0.3, name="drop_3")(conv1_out)
        conv1_pool_out = AveragePooling1D(2, padding='same', name="pool_3")(conv1_dropout_out)
        conv1_flatten_out = Flatten(name="flatten_3")(conv1_pool_out)
        for_dense.append(conv1_flatten_out)

        # 第二个卷积层
        conv2_out = Conv1D(70, 5, activation='relu', name="conv_5")(input_c)
        conv2_dropout_out = Dropout(0.3, name="drop_5")(conv2_out)
        conv2_pool_out = AveragePooling1D(2, padding='same', name="pool_5")(conv2_dropout_out)
        conv2_flatten_out = Flatten(name="flatten_5")(conv2_pool_out)
        for_dense.append(conv2_flatten_out)

        # 第三个卷积层
        conv3_out = Conv1D(40, 7, activation='relu', name="conv_7")(input_c)
        conv3_dropout_out = Dropout(0.3, name="drop_7")(conv3_out)
        conv3_pool_out = AveragePooling1D(2, padding='same', name="pool_7")(conv3_dropout_out)
        conv3_flatten_out = Flatten(name="flatten_7")(conv3_pool_out)
        for_dense.append(conv3_flatten_out)

    # 拼接卷积层的输出
    if len(for_dense) == 1:
        concat_out = for_dense[0]
    else:
        concat_out = concatenate(for_dense)

    for_dense1 = []
    # 第一个全连接层
    dense0_out = Dense(80, activation='relu', name="dense_0")(concat_out)
    dense0_dropout_out = Dropout(0.3, name="drop_d0")(dense0_out)
    for_dense1.append(dense0_dropout_out)

    # 添加delta Gb输入
    if 'G' in TYPE:
        for_dense1.append(input_g)

    if len(for_dense1) == 1:
        concat1_out = for_dense1[0]
    else:
        concat1_out = concatenate(for_dense1)

    # 第二个全连接层
    dense1_out = Dense(80, activation='relu', name="dense_1")(concat1_out)
    dense1_dropout_out = Dropout(0.3, name="drop_d1")(dense1_out)

    # 第三个全连接层
    dense2_out = Dense(60, activation='relu', name="dense_2")(dense1_dropout_out)
    dense2_dropout_out = Dropout(0.3, name="drop_d2")(dense2_out)

    # 输出层
    output = Dense(1, name="output")(dense2_dropout_out)

    # 构建模型
    model = Model(inputs=inputs, outputs=[output])
    return model
