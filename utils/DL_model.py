from tensorflow.keras.layers import *
from tensorflow.keras.models import *
import numpy as np
import tensorflow as tf

def get_model(max_len_en,max_dep):


    sequences = Input(shape=(max_len_en,max_dep,),name="input_onehot")
#    crispron= Input(shape=(1,))
    energy = Input(shape=(1,),name="input_energy")
    for_dense = list()
    conv1_out = Conv1D(100, 3, activation='relu', padding="valid", name="conv_3")(sequences)
    conv1_dropout_out = Dropout(0.3, name="drop_3")(conv1_out)
    conv1_pool_out = AveragePooling1D(2, padding='valid', name="pool_3")(conv1_dropout_out)
    conv1_flatten_out = Flatten(name="flatten_3")(conv1_pool_out)
    for_dense.append(conv1_flatten_out)

    conv2_out = Conv1D(70, 5, activation='relu', padding="valid", name="conv_5")(sequences)
    conv2_dropout_out = Dropout(0.3, name="drop_5")(conv2_out)
    conv2_pool_out = AveragePooling1D(2, padding='valid', name="pool_5")(conv2_dropout_out)
    conv2_flatten_out = Flatten(name="flatten_5")(conv2_pool_out)
    for_dense.append(conv2_flatten_out)

    conv3_out = Conv1D(40, 7, activation='relu', padding="valid", name="conv_7")(sequences)
    conv3_dropout_out = Dropout(0.3, name="drop_7")(conv3_out)
    conv3_pool_out = AveragePooling1D(2, padding='valid', name="pool_7")(conv3_dropout_out)
    conv3_flatten_out = Flatten(name="flatten_7")(conv3_pool_out)
    for_dense.append(conv3_flatten_out)

    focus_sequence = Lambda(lambda x: x[:, 18:28, :], name="focus_seq")(sequences)
    focus_conv_out = Conv1D(10, 2, activation='relu', padding="valid", name="focus_conv")(focus_sequence)
    focus_dropout_out = Dropout(0.3, name="focus_drop")(focus_conv_out)
    focus_pool_out = AveragePooling1D(2, padding='valid', name="focus_pool")(focus_dropout_out)
    focus_flatten_out = Flatten(name="focus_flatten")(focus_pool_out)
    for_dense.append(focus_flatten_out)

    if len(for_dense) == 1:
        concat_out = for_dense[0]
    else:
         concat_out = concatenate(for_dense)

    for_dense1 = list()
#first dense (fully connected) layer
    dense0_out = Dense(80, activation='relu',name="dense_0")(concat_out)
    dense0_dropout_out = Dropout(0.3, name="drop_d0")(dense0_out)
    for_dense1.append(dense0_dropout_out)

    #Gb input used raw
    for_dense1.append(energy)
#    for_dense1.append(crispron)
    if len(for_dense1) == 1:
        concat1_out = for_dense1[0]
    else:
        concat1_out = concatenate(for_dense1)


#first dense (fully connected) layer
    dense1_out = Dense(80, activation='relu',name="dense_1")(concat1_out)
    dense1_dropout_out = Dropout(0.3, name="drop_d1")(dense1_out)

#second dense (fully connected) layer
    dense2_out = Dense(60, activation='relu',name="dense_2")(dense1_dropout_out)
    dense2_dropout_out = Dropout(0.3, name="drop_d2")(dense2_out)

#output layer
    output = Dense(1,name="output")(dense2_dropout_out)

#model construction
    model= Model([sequences,energy], output)

    return model
