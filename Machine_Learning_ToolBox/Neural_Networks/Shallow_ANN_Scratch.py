#Author: Ciara Conway
#Data: 12/12/21

import numpy as np
import pandas as pd
import math
from random import random
from random import seed
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
import matplotlib.pyplot as plt
import seaborn as sns
import time

def initializeNN(num_input, num_hidd, num_out):
    network = list()
    hid_lay = [{'weights':[random() for i in range(num_input + 1)]} for i in range(num_hidd)]
    network.append(hid_lay)
    out_lay = [{'weights':[random() for i in range(num_hidd + 1)]} for i in range(num_out)]
    network.append(out_lay)
    return network

###FP###
def activation(w, inputs): #The input could be a row from our training dataset, as in the case of the hidden layer. It may also be the outputs from each neuron in the hidden layer, in the case of the output layer.
    active = w[-1]
    for i in range(len(w)-1):
        active += w[i]*inputs[i]
    return active

def transfer(activation):
    return 1.0/(1.0 + math.exp(-activation))

def f_prop(network, row): #row from training set or from a hidden layer
    input = row
    for layer in network:
        new_in = []
        for neuron in layer:
            active = activation(neuron['weights'], input)
            neuron['output'] = transfer(active)
            new_in.append(neuron['output'])
        input = new_in
    return input #output from the nn fp

###BP###
def der_transfer(output):
    return output*(1.0 - output)

def b_prop(network, expected):
    for i in reversed(range(len(network))):
        layer = network[i]
        errors = []
        if i !=len(network)-1:
            for j in range(len(layer)):
                error = 0.0
                for neuron in network[i+1]:
                    error += (neuron['weights'][j]*neuron['delta'])
                    errors.append(error)
        else:
            for j in range(len(layer)):
                neuron = layer[j]
                errors.append(neuron['output'] - expected[j])
        for j in range(len(layer)):
            neuron = layer[j]
            neuron['delta'] = errors[j]*der_transfer(neuron['output'])

#update network weights from BP
def up_weights(network, row, l_rate):
    for i in range(len(network)):
        input = row[:-1]
        if i != 0:
            input = [neuron['output'] for neuron in network[i - 1]]
        for neuron in network[i]:
            for j in range(len(input)):
                neuron['weights'][j] -= l_rate*neuron['delta']*input[j]
            neuron['weights'][-1] -= l_rate*neuron['delta']

#epoch means training the neural network with all the training data for one cycle. In an epoch, we use all of the data exactly once.
# A forward pass and a backward pass together are counted as one pass:
# An epoch is made up of one or more batches, where we use a part of the dataset to train the neural network

#training network
def train_net(network, data, l_rate, num_epoch, num_out):
    data = np.array(data)
    data = data.tolist()
    epoch_plot = []
    sum_err_plot = []
    for epoch in range(num_epoch):
        sum_err = 0
        for row in data:
            outputs = f_prop(network, row)
            expect = [0 for i in range(num_out)]
            expect[int(row[-1])] = 1
            sum_err += sum([(expect[i] - outputs[i])**2 for i in range(len(expect))])
            b_prop(network, expect)
            up_weights(network, row, l_rate)
            epoch_plot.append(epoch)
            sum_err_plot.append(sum_err)
        print('>epoch=%d, lrate=%.3f, error=%.3f' % (epoch, l_rate, sum_err))
    plt.scatter(epoch_plot, sum_err_plot)
    plt.show()
    return network, epoch_plot, sum_err_plot


#testing predictions with NN
def predict(network, row):
    outputs = f_prop(network, row)
    return outputs.index(max(outputs))

#BP with stochastic gradient descent
def b_prop_sto(data, test, l_rate, num_epoch, num_hid):
    num_in = len(data[0] - 1)
    num_out = 1
    network = initializeNN(num_in, num_hid, num_out)
    train_net(network, data, l_rate, num_epoch, num_out)
    predictions = []
    for row in test:
        prediction = predict(network, row)
        predictions.append(prediction)
    return predictions

#test predictions of network
seed(1)
#normalize data
def minmax_data(data):
    data_norm = (data - data.min(axis=0))/(data.max(axis=0)-data.min(axis=0))
    return data_norm

def accuracy(actual, predicted):
    correct = 0
    for i in range(len(actual)):
        if actual[i] == predicted[i]:
            correct += 1
    return correct/float(len(actual))*100
def training_test_split(data_norm, target, seed):
    X_tr, X_te, y_tr, y_te = train_test_split(data_norm, target, test_size=0.2, random_state=seed)
    return X_tr, X_te, y_tr, y_te
def nn_algorithm(data, target, hidd, num_inputs, num_outputs):
    network_new = None
    network_list = []
    cycle = 1
    training_list = []
    accuracy_list = []
    test_actual = []
    test_predict = []
    X_tr, X_te, y_tr, y_te = training_test_split(data, target, 1)
    X_tr = minmax_data(X_tr)
    X_te = minmax_data(X_te)
    data_norm_ac_te = np.column_stack((X_te, y_te))
    data_norm_ac_te = data_norm_ac_te.tolist()
    while cycle < 6:
        train_list = []
        actual = []
        predicted = []
        VX_tr, VX_te, Vy_tr, Vy_te = training_test_split(X_tr, y_tr, cycle)
        data_norm_tr = np.column_stack((VX_tr, Vy_tr))
        data_norm_te = np.column_stack((VX_te, Vy_te))
        data_norm_te = data_norm_te.tolist()
        network = initializeNN(num_inputs,hidd,num_outputs)
        network, epoch_plot, sum_err_plot = train_net(network,data_norm_tr,0.1,20,num_outputs)
        for layer in network:
            train_list.append(layer)
            print(layer)
        for row in data_norm_te:
            prediction = predict(network, row)
            actual.append(int(row[-1]))
            predicted.append(prediction)
            print('Expected=%d, Got=%d' % (int(row[-1]), prediction))
        accur = accuracy(actual, predicted)
        training_list.append((sum_err_plot[-1]))
        accuracy_list.append(accur)
        network_new = network
        network_list.append(network)
        cycle += 1
    max_acc = max(accuracy_list)
    best_network = accuracy_list.index(max_acc)
    for row in data_norm_ac_te:
        prediction = predict(network_list[best_network], row)
        test_actual.append(int(row[-1]))
        test_predict.append(prediction)
        print('Expected=%d, Got=%d' % (int(row[-1]), prediction))
    accur_test = accuracy(test_actual, test_predict)
    return network_new, training_list, accuracy_list, test_actual, test_predict, accur_test

def feature_selection(selection, features):
    if selection == "univariate":
        data, target = datasets.load_breast_cancer(return_X_y=True, as_frame=True)
        best_feat = SelectKBest(score_func=chi2, k=10)
        fit = best_feat.fit(data, target)
        scores = pd.DataFrame(fit.scores_)
        cols = pd.DataFrame(data.columns)
        score_feat = pd.concat([cols, scores], axis=1)
        score_feat.columns = ["Col", "Score"]
        largest = score_feat.nlargest(features,"Score")
        largest = largest.index.values.tolist()
        all_largest = score_feat.nlargest(30, "Score")
        #all_largest.to_csv("/Users/ciaraconway/Documents/mlu.csv")
        #print(score_feat.nlargest(10,"Score"))
        ###clean data###
        data, target = datasets.load_breast_cancer(return_X_y=True)
        data = data[:, largest]
    elif selection == "correlation":
        data, target = datasets.load_breast_cancer(return_X_y=True, as_frame=True)
        corrmat = data.corr()
        top_corr_features = corrmat.index
        plt.figure(figsize=(20,20))
        g = sns.heatmap(data[top_corr_features].corr(), annot=True, cmap="RdYlGn")
        data, target = datasets.load_breast_cancer(return_X_y=True)
    else:
        data, target = datasets.load_breast_cancer(return_X_y=True)
    return data, target

def write_out_file(textFile, file):
    out_file = open(file, "w")
    out_file.write(textFile)
    out_file.close()
    print("Outfile: ", file, " has been saved")

def main():
    textFile = ""
    ###All you have to do is change file path for text file and click run####
    #User parameters to set
    selection = "correlation" #can set to "correlation" to generate heatmap and run with 30 features
    #User parameters to set
    alphas = [2, 6, 8, 10]
    features = [5, 10, 15, 20, 30]
    for i, alpha in enumerate(alphas):
        for j, feature in enumerate(features):
            start = time.time()
            data, target = feature_selection(selection, feature)
            num_sam = len(data[:,0])
            num_inputs = len(data[0])
            num_outputs = 2
            hidd = round(num_sam/(alpha*(num_inputs+num_outputs)))
            network_new, training_list, accuracy_list, test_actual, test_predict, accur_test = nn_algorithm(data, target, hidd, num_inputs, num_outputs)
            #df = pd.DataFrame(list(zip(training_list, accuracy_list)),columns =['validation training error', 'accuracy'])
            #df.to_csv("/Users/ciaraconway/Documents/accuracyerror" + str(feature) + ".csv")
            print(feature)
            print(accur_test)
            end = time.time()
            total = end - start
            print(total)
            textFile += "alpha: " + str(alpha) + "\n" + "num features: " + str(feature) + "\n" + "Accuracy: " + str(accur_test) + "\n" + "time: " + str(total) + "\n\n"
    write_out_file(textFile, "/Users/ciaraconway/Documents/NN_final_data_10_0.5.txt")

if __name__ == '__main__':
    main()
