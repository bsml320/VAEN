import os
import sys
import random as rn
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
os.environ['KERAS_BACKEND']='tensorflow'

import tensorflow as tf
from keras import backend as K

## obtain reproducible results - START
os.environ['PYTHONHASHSEED'] = '0'
# The below is necessary for starting Numpy generated random numbers
# in a well-defined initial state.
#np.random.seed(40)

np.random.seed(int(sys.argv[1]))


# The below is necessary for starting core Python generated random numbers
# in a well-defined state.
##rn.seed(12345)

rn.seed(int(sys.argv[1]))

# The below tf.set_random_seed() will make random number generation
# in the TensorFlow backend have a well-defined initial state.
# For further details, see:
# https://www.tensorflow.org/api_docs/python/tf/set_random_seed
session_conf = tf.ConfigProto(intra_op_parallelism_threads=1,
                              inter_op_parallelism_threads=1)

#tf.set_random_seed(1234)

tf.set_random_seed(int(sys.argv[1]))

sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
K.set_session(sess)
## obtain reproducible results - END

from keras.layers import Input, Dense, Lambda, Layer, Activation, Dropout
from keras.layers.normalization import BatchNormalization
from keras.models import Model

from keras import metrics, optimizers
from keras.callbacks import Callback
import keras

import pydot
import graphviz
from keras.utils import plot_model
from keras_tqdm import TQDMNotebookCallback
from IPython.display import SVG
from keras.utils.vis_utils import model_to_dot

print(keras.__version__)
tf.__version__

work_dir = sys.argv[2]
train_file_path = sys.argv[3] #ccle.zeroone_5000_0.2.tsv
val_file_1_path = sys.argv[4] #PANCAN.zeroone_5000_0.2.tsv
train_latent_file = sys.argv[5] #CCLE_latent_5000_0.2.tsv
train_weight_file = sys.argv[6] #CCLE_gene_weights_5000_0.2.tsv
val_1_pred_file = sys.argv[7] #PANCAN_prediction_5000_0.2.tsv
encoder_file = sys.argv[8] # CCLE_encoder_onehidden_vae.hdf5
decoder_file = sys.argv[9] # CCLE_decoder_onehidden_vae.hdf5

print ("work_dir: %s" % work_dir)

plt.style.use('seaborn-notebook')
sns.set(style="white", color_codes=True)
sns.set_context("paper", rc={"font.size":14,"axes.titlesize":15,"axes.labelsize":20,'xtick.labelsize':14, 'ytick.labelsize':14})


def sampling(args):
    #import tensorflow as tf
    z_mean, z_log_var = args
    epsilon = K.random_normal(shape=tf.shape(z_mean), mean=0.,stddev=epsilon_std)
    z = z_mean + K.exp(z_log_var / 2) * epsilon
    return z


class CustomVariationalLayer(Layer):
    def __init__(self, **kwargs):
        self.is_placeholder = True
        super(CustomVariationalLayer, self).__init__(**kwargs)
    def vae_loss(self, x_input, x_decoded):
        reconstruction_loss = original_dim * metrics.binary_crossentropy(x_input, x_decoded)
        kl_loss = - 0.5 * K.sum(1 + z_log_var_encoded - K.square(z_mean_encoded) - 
                                K.exp(z_log_var_encoded), axis=-1)
        return K.mean(reconstruction_loss + (K.get_value(beta) * kl_loss))
    def call(self, inputs):
        x = inputs[0]
        x_decoded = inputs[1]
        loss = self.vae_loss(x, x_decoded)
        self.add_loss(loss, inputs=inputs)
        return x


class WarmUpCallback(Callback):
    def __init__(self, beta, kappa):
        self.beta = beta
        self.kappa = kappa
    def on_epoch_end(self, epoch, logs={}):
        if K.get_value(self.beta) <= 1:
            K.set_value(self.beta, K.get_value(self.beta) + self.kappa)



####rnaseq_file = os.path.join('/data1_2/jiap/projects/18-CCLE-VAE/SCALE', 'ccle.zeroone_5000_0.2.tsv')
rnaseq_file = os.path.join(work_dir, train_file_path)
rnaseq_df = pd.read_table(rnaseq_file, index_col=0)
print(rnaseq_df.shape)
rnaseq_df.head(2)

####val_file = os.path.join('/data1_2/jiap/projects/18-CCLE-VAE/SCALE', 'PANCAN.zeroone_5000_0.2.tsv')
val_file_1 = os.path.join(work_dir, val_file_1_path)
val_df_1 = pd.read_table(val_file_1, index_col=0)
print(val_df_1.shape)
val_df_1.head(2)


test_set_percent = 0.1
rnaseq_test_df = rnaseq_df.sample(frac=test_set_percent)
rnaseq_train_df = rnaseq_df.drop(rnaseq_test_df.index)


original_dim = rnaseq_df.shape[1]
latent_dim = 100

batch_size = 100
epochs = 100
learning_rate = 0.0005

epsilon_std = 1.0
beta = K.variable(0)
kappa = 1


rnaseq_input = Input(shape=(original_dim, ))

z_mean_dense_linear = Dense(latent_dim, kernel_initializer='glorot_uniform')(rnaseq_input)
z_mean_dense_batchnorm = BatchNormalization()(z_mean_dense_linear)
z_mean_encoded = Activation('sigmoid')(z_mean_dense_batchnorm)

z_log_var_dense_linear = Dense(latent_dim, kernel_initializer='glorot_uniform')(rnaseq_input)
z_log_var_dense_batchnorm = BatchNormalization()(z_log_var_dense_linear)
z_log_var_encoded = Activation('sigmoid')(z_log_var_dense_batchnorm)

z = Lambda(sampling, output_shape=(latent_dim, ))([z_mean_encoded, z_log_var_encoded])

drop_layer = Dropout(rate = 0.2, noise_shape = None)(z)			####### modified added

decoder_to_reconstruct = Dense(original_dim, kernel_initializer='glorot_uniform', activation='sigmoid')

rnaseq_reconstruct = decoder_to_reconstruct(drop_layer)			####### modified

adam = optimizers.Adam(lr=learning_rate)
vae_layer = CustomVariationalLayer()([rnaseq_input, rnaseq_reconstruct])
vae = Model(rnaseq_input, vae_layer)
vae.compile(optimizer=adam, loss=None, loss_weights=[beta])

vae.summary()

#output_model_file = os.path.join('/data1_2/jiap/projects/18-CCLE-VAE/SCALE', 'onehidden_vae_architecture.png')
#plot_model(vae, to_file=output_model_file)
#SVG(model_to_dot(vae).create(prog='dot', format='svg'))

hist = vae.fit(np.array(rnaseq_train_df),
               shuffle=True,
               epochs=epochs,
               verbose=0,
               batch_size=batch_size,
               validation_data=(np.array(rnaseq_test_df), None),
               callbacks=[WarmUpCallback(beta, kappa),
               TQDMNotebookCallback(leave_inner=True, leave_outer=True)])



history_df = pd.DataFrame(hist.history)
loss_log_file = os.path.join(work_dir+"/result/", sys.argv[1]+"."+"sampling.K.loss.log.txt")
history_df.to_csv(loss_log_file, sep='\t')


encoder = Model(rnaseq_input, z_mean_encoded)
encoded_rnaseq_df = encoder.predict_on_batch(rnaseq_df)
encoded_rnaseq_df = pd.DataFrame(encoded_rnaseq_df, index=rnaseq_df.index)

encoded_rnaseq_df.columns.name = 'sample_id'
encoded_rnaseq_df.columns = encoded_rnaseq_df.columns + 1
####encoded_file = os.path.join('/data1_2/jiap/projects/18-CCLE-VAE/SCALE', 'CCLE_latent_5000_0.2.tsv')
encoded_file = os.path.join(work_dir+"/result/", sys.argv[1]+"."+train_latent_file)
encoded_rnaseq_df.to_csv(encoded_file, sep='\t')

decoder_input = Input(shape=(latent_dim, ))  # can generate from any sampled z vector
_x_decoded_mean = decoder_to_reconstruct(decoder_input)
decoder = Model(decoder_input, _x_decoded_mean)


#encoder_model_file = os.path.join('/data1_2/jiap/projects/18-CCLE-VAE/SCALE', 'CCLE_encoder_onehidden_vae.hdf5')
#decoder_model_file = os.path.join('/data1_2/jiap/projects/18-CCLE-VAE/SCALE', 'CCLE_decoder_onehidden_vae.hdf5')

encoder_model_file = os.path.join(work_dir+"/result/", sys.argv[1]+"."+encoder_file)
decoder_model_file = os.path.join(work_dir+"/result/", sys.argv[1]+"."+decoder_file)

encoder.save(encoder_model_file)
decoder.save(decoder_model_file)


sum_node_activity = encoded_rnaseq_df.sum(axis=0).sort_values(ascending=False)
print(sum_node_activity.head(10))
sum_node_activity.tail(10)



weights = []
for layer in decoder.layers:
	weights.append(layer.get_weights())


weight_layer_df = pd.DataFrame(weights[1][0], columns=rnaseq_df.columns, index=range(1, 101))
weight_layer_df.index.name = 'encodings'
weight_layer_df.head(2)

####weight_file = os.path.join('/data1_2/jiap/projects/18-CCLE-VAE/SCALE', 'CCLE_gene_weights_5000_0.2.tsv')
weight_file = os.path.join(work_dir+"/result/", sys.argv[1]+"."+train_weight_file)
weight_layer_df.to_csv(weight_file, sep='\t')

###################################
encoded_val_df = encoder.predict_on_batch(val_df_1)
encoded_val_df = pd.DataFrame(encoded_val_df, index=val_df_1.index)

encoded_val_df.columns.name = 'sample_id'
encoded_val_df.columns = encoded_val_df.columns + 1
####encoded_file = os.path.join('/data1_2/jiap/projects/18-CCLE-VAE/SCALE', 'PANCAN_prediction_5000_0.2.tsv')
encoded_file = os.path.join(work_dir+"/result/", sys.argv[1]+"."+val_1_pred_file)
encoded_val_df.to_csv(encoded_file, sep='\t')
encoded_val_df.head(2)


