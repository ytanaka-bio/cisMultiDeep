#module load cuda
#export LD_LIBRARY_PATH=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/cudacore/11.0.2/targets/x86_64-linux/lib/stubs/:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/home/ytanaka/program/cellranger-7.0.0/external/anaconda/lib/:$LD_LIBRARY_PATH

#0 initialization
import sys
import argparse
import warnings
warnings.filterwarnings("ignore")

#0.1 parse argument
parser = argparse.ArgumentParser(description='Run deep learning classifier',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#0.2 parameter
parser.add_argument('-i',nargs=1,help='data filename',type=argparse.FileType('r'),dest="data_file",required=True)
parser.add_argument('-p',nargs=1,help='phenotype filename',type=argparse.FileType('r'),dest="phenotype_file",required=True)
parser.add_argument('-o',nargs=1,help='output prefix (<prefix>_test.csv (with -l option), <prefix>_predict.csv (with -l option), <prefix>_result.csv, and <prefix>_model/)',type=str,dest="out_file",required=True)
parser.add_argument('-c', help='Scale input data?',nargs = "?", default=True,dest="scale")
parser.add_argument('-v', help='ratio of validation dataset', default=0.2, type=float, dest="ratio_val")
parser.add_argument('-s', help='ratio of test dataset', default=0.2, type=float, dest="ratio_test")
parser.add_argument('-e', help='epoch',type=int, default=10,dest='epoch')
parser.add_argument('-x', help="Max trial",type=int, default=50,dest="max_trial")
parser.add_argument('-l', help='Evaluate model?',action='store_true',default=False,dest="eval")
parser.add_argument('-n', help='Size of background dataset',type=int, default=1000,dest="size")
parser.add_argument('-m', help='metrics',default='accuracy',type=str,dest="metrics")
parser.add_argument('-u', help='activation function in output layer', nargs=1,type=str,default='softmax',dest="act_out")
parser.add_argument('-t', help="number of threads",default=1,type=int,dest="thread")
parser.add_argument('-d', help='directory for model tuning', default='model_autoencoder_hyperparam_tuning_trans', type=str, dest="dir")


#0.3. get arguments of selected mode
args = parser.parse_args()

#0.4. import required packages
import numpy as np
import pandas as pd

import tensorflow as tf
tf.config.threading.set_inter_op_parallelism_threads(args.thread)
tf.config.threading.set_intra_op_parallelism_threads(args.thread)

import shap
from scipy.stats import zscore
from tensorflow import feature_column, keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Model
from sklearn.model_selection import train_test_split
import keras_tuner as kt
from keras.models import load_model

#1 define functions
#1.1. deep learning model
def deephypermodel(hp):
  depth = hp.Int('layer', min_value = 2, max_value = 10, step = 1)
  hidden_dim = hp.Int('units', min_value = 50, max_value = 500, step = 50)
  act_fct = hp.Choice('act', ['relu','sigmoid','tanh'])
  drop = hp.Float('dropout', 0, 0.5, step=0.1, default=0.5)
  opt = hp.Choice('opt', ['adam','sgd'])
  loss_fct = hp.Choice('loss', ['binary_crossentropy','categorical_crossentropy','mse'])
  learning_rate = hp.Choice('learning_rate', values=[0.1, 0.01, 0.001, 0.0001, 0.00001])
  model = build_model(depth, hidden_dim, act_fct, opt, loss_fct, drop, learning_rate, feature_layer, out_dim)
  return model

#1.2. Build model
def build_model(depth, hidden_dim, act_fct, opt, loss_fct, drop, learning_rate, feature_layer, out_dim):
  model = tf.keras.Sequential(feature_layer) 
  for i in range(1, depth):
    model.add(tf.keras.layers.Dense(hidden_dim, activation = act_fct))
  model.add(tf.keras.layers.Dropout(drop))
  model.add(tf.keras.layers.Dense(out_dim, activation = 'sigmoid')) # CHANGE ACCORDING TO OUTPUT WE WANT TO PREDICT
  if opt == "adam":
    optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate)
  elif opt == "sgd":
    optimizer=tf.keras.optimizers.SGD(learning_rate=learning_rate)
  model.compile(
    optimizer = optimizer,
    loss = loss_fct,
    metrics = [args.metrics]) #hp.Choice('metric', ['accuracy', 'precision', 'recall'])
  return model

#1.3. convert dataset format
def df_to_dataset(dataframe, dindex, feature, batch_size=2048):
  dataframe = dataframe.copy()
  labels = dataframe.loc[:,feature]
  ds = tf.data.Dataset.from_tensor_slices((dict(dataframe.loc[:,dindex]), labels))
  ds = ds.batch(batch_size)
  return ds

#2 main
#2.1. read files
#2.1.1. data file
print("Read",args.data_file[0])
dataframe = pd.read_csv(args.data_file[0],index_col=0)
if args.scale == True:
  dataframe = dataframe.apply(zscore)

dindex = dataframe.columns
reduced_dim = len(dindex)
output = args.out_file[0]

#2.1.2. phenotype or model file
#2.1.2.1 phenotype file
print("Read",args.phenotype_file[0])
dataout = pd.read_csv(args.phenotype_file[0],index_col=0)
pheno = dataout.columns
out_dim = len(pheno)

#2.2. Prepare dataset
if args.eval == True:
  X_train, X_test, Y_train, Y_test = train_test_split(dataframe, dataout, test_size=args.ratio_test)
  X_train, X_val, Y_train, Y_val = train_test_split(X_train, Y_train, test_size=args.ratio_val)

X_train2, X_val2, Y_train2, Y_val2 = train_test_split(dataframe, dataout, test_size=args.ratio_val)
      
#2.3. model ceation
#2.3.1. model parameter optimization
print("Model optimization")
feature_layer = tf.keras.layers.InputLayer(input_shape=(len(dindex),))
tuner = kt.RandomSearch(deephypermodel, objective='val_accuracy', executions_per_trial = 3, max_trials=args.max_trial, directory = args.dir, project_name = args.out_file[0])
tuner.search(X_train2,Y_train2, epochs=args.epoch, validation_data=(X_val2,Y_val2)) 
best_hps=tuner.get_best_hyperparameters(num_trials=1)[0]
  
#2.3.2. model build
print("Build the model")
model = tuner.hypermodel.build(best_hps)

if args.eval == True:
  print("Assess the model")
  model.fit(X_train, Y_train, validation_data=(X_val, Y_val), epochs=args.epoch)
  y_score = pd.DataFrame(model.predict(test_ds))
  y_test = test.loc[:,pheno]
  y_score.index = y_test.index
  y_score.columns = y_test.columns
  output1 = output + "_test.csv"
  y_test.to_csv(output1)
  output2 = output + "_predict.csv"
  y_score.to_csv(output2)
  
print("Output the prediction")
model.fit(X_train2, Y_train2, validation_data=(X_val2, Y_val2), epochs=args.epoch)
all_score = pd.DataFrame(model.predict(dataframe))
all_score.index = dataframe.index
all_score.columns = pheno
output5 = output + "_result.csv"
all_score.to_csv(output5)

print("Save the model")
output3 = output + "_model"
model.save(filepath=output3,overwrite=True)

print("Calculate SHAP values")
output4 = output + "_expected.csv"
background = X_train2.iloc[np.random.choice(X_train2.shape[0], args.size, replace=False)]
explainer = shap.DeepExplainer(model,np.array(background))
shap_values = explainer.shap_values(np.array(dataframe))
#expected_value = explainer.expected_value
#expected_value = pd.DataFrame(expected_value)
#expected_value.to_csv(output4)

print("Calculate mean abolute SHAP values")
output6 = output + "_meanshap.csv"
shap_values = np.absolute(shap_values)
mean_shap = pd.DataFrame(0, index=dindex,columns=pheno)
for i in range(0,len(pheno)):
  df = pd.DataFrame(shap_values[i],columns=dindex)
  mean_shap.iloc[:,i] = mean_shap.iloc[:,i] + df.sum()

mean_shap = mean_shap / args.size
mean_shap.to_csv(output6)

