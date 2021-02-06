################################
# Example: python3 11 GSE20194.RANK.tsv  path/to/result(folder)
################################

import sys
import os
os.environ["KERAS_BACKEND"] = "tensorflow"

import pandas as pd

from keras.models import load_model
import keras.backend as K
import tensorflow as tf
config = tf.ConfigProto()
config.gpu_options.allow_growth=True
sess = tf.Session(config=config)
K.set_session(sess)

ksigmoid = sys.argv[1]

####
input_file = sys.argv[2]
encoder_hdf5 = os.path.join(sys.argv[3], ksigmoid+".encoder.hdf5")

print(encoder_hdf5)

print("[INFO] Loading model...")
encoder = load_model(encoder_hdf5,compile=False)

print("[INFO] Loading file...")
input_df = pd.read_table(input_file, index_col=0)
#input_df = input_df.T
print("--->Input Shape:",input_df.shape)

encoded_df = encoder.predict(input_df)
encoded_df = pd.DataFrame(encoded_df, index=input_df.index)

encoded_file = os.path.join(input_file[:-4]+".latent.tsv")
encoded_df.to_csv(encoded_file, sep='\t')

print('[INFO] Done !')
