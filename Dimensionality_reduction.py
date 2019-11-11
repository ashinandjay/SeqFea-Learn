# This code is used for feature Reduction methods, 
# The output file will be generated in same folder of this code
# -------------------------------------------------------------

# load packages:
import sys,os
import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.decomposition import KernelPCA
from sklearn.manifold import LocallyLinearEmbedding
from sklearn.decomposition import NMF
from sklearn.decomposition import TruncatedSVD
from sklearn.manifold import MDS
from sklearn.decomposition import FactorAnalysis
from sklearn.decomposition import FastICA
from sklearn.cluster import FeatureAgglomeration
from sklearn.random_projection import GaussianRandomProjection
from sklearn.random_projection import SparseRandomProjection

# find the path
Script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

# choose the method
option = sys.argv[1]

# read number of clustering
selected_number = sys.argv[2]

# read feature vectors
feature = sys.argv[3]
feature_o = pd.read_csv(feature,delim_whitespace=True,index_col=0,float_precision='round_trip')


# Feature selection tsne
def tsneselection(k,X):
    X=np.array(X)
    Tsneresult =TSNE(n_components=k).fit_transform(X)
    np.savetxt("TSNE_out.csv", Tsneresult, delimiter=",")    
    return None 

# Feature selection PCA
def pcaselection(k,X):
    X=np.array(X)
    pcaresult = PCA(n_components=k).fit_transform(X)
    np.savetxt("PCA_out.csv", pcaresult, delimiter=",")
    return None

# Feature selection kernelPCA
def kernelPCA(k,X):
    X=np.array(X)
    kpca=KernelPCA(n_components=k)
    kpcaresult=kpca.fit_transform(X)
    np.savetxt("KernelPCA_out.csv", kpcaresult, delimiter=",")
    return None

# Feature selection Locally-linear embedding
def lle(k,X):
    X=np.array(X)
    embedding = LocallyLinearEmbedding(n_components=k)
    lleresult = embedding.fit_transform(X)
    np.savetxt("LLE_out.csv", lleresult, delimiter=",")
    return None

# Feature selection SVD
def svd(k,X):
    X=np.array(X)
    SVDMethod = TruncatedSVD(n_components=k)
    svdresult=SVDMethod.fit_transform(X)
    np.savetxt("SVD_out.csv", svdresult, delimiter=",")
    return None

# Feature selection NMF
def nmf(k,X):
    X=np.array(X)
    NMFmodel=NMF(n_components=k)
    nmfresult= NMFmodel.fit_transform(X)
    np.savetxt("NMF_out.csv", nmfresult, delimiter=",")
    return None

# Feature selection MDS
def mds(k,X):
    MDSmodel = MDS(n_components=k)
    mdsresult = MDSmodel.fit_transform(X)
    np.savetxt("MDS_out.csv", mdsresult, delimiter=",")
    return None

# Feature selection ICA
def ICA(data,n_components):
    ica = FastICA(n_components=n_components)
    X_transformed = ica.fit_transform(data)
    np.savetxt("ICA_out.csv", X_transformed, delimiter=",")
    return None

# Feature selection FA
def FA(data,n_components):
    fa= FactorAnalysis(n_components=n_components)
    FAresult=fa.fit_transform(data)
    np.savetxt("FA_out.csv", FAresult, delimiter=",")
    return None

# Feature selection Agglomerate features
def Aggselection(k,X):
    X=np.array(X)
    Aggresult =FeatureAgglomeration(n_clusters=k).fit_transform(X)
    np.savetxt("Agglomeration_out.csv", Aggresult, delimiter=",")    
    return None 

# Feature selection Gaussian random projection
def Gaussianselection(k,X):
    X=np.array(X)
    gaussiansresult=GaussianRandomProjection(n_components=k).fit_transform(X)
    np.savetxt("Gaussian_out.csv", gaussiansresult, delimiter=",")
    return None

# Feature selection Sparse random projection
def Sparseselection(k,X):
    X=np.array(X)
    sparsesresult=SparseRandomProjection(n_components=k).fit_transform(X)
    np.savetxt("Sparse_out.csv", sparsesresult, delimiter=",")
    return None

class Autoencoder(object):

    def __init__(self, n_input, n_hidden, transfer_function=tf.nn.softplus, optimizer = tf.train.AdamOptimizer()):
        self.n_input = n_input
        self.n_hidden = n_hidden
        self.transfer = transfer_function

        network_weights = self._initialize_weights()
        self.weights = network_weights

        # model
        self.x = tf.placeholder(tf.float32, [None, self.n_input])
        self.hidden = self.transfer(tf.add(tf.matmul(self.x, self.weights['w1']), self.weights['b1']))
        self.reconstruction = tf.add(tf.matmul(self.hidden, self.weights['w2']), self.weights['b2'])

        # cost
        self.cost = 0.5 * tf.reduce_sum(tf.pow(tf.subtract(self.reconstruction, self.x), 2.0))
        self.optimizer = optimizer.minimize(self.cost)

        init = tf.global_variables_initializer()
        self.sess = tf.Session()
        self.sess.run(init)


    def _initialize_weights(self):
        all_weights = dict()
        all_weights['w1'] = tf.get_variable("w1", shape=[self.n_input, self.n_hidden],
            initializer=tf.contrib.layers.xavier_initializer())
        all_weights['b1'] = tf.Variable(tf.zeros([self.n_hidden], dtype=tf.float32))
        all_weights['w2'] = tf.Variable(tf.zeros([self.n_hidden, self.n_input], dtype=tf.float32))
        all_weights['b2'] = tf.Variable(tf.zeros([self.n_input], dtype=tf.float32))
        return all_weights

    def partial_fit(self, X):
        cost, opt = self.sess.run((self.cost, self.optimizer), feed_dict={self.x: X})
        return cost

    def calc_total_cost(self, X):
        return self.sess.run(self.cost, feed_dict = {self.x: X})

    def transform(self, X):
        return self.sess.run(self.hidden, feed_dict={self.x: X})

    def generate(self, hidden = None):
        if hidden is None:
            hidden = self.sess.run(tf.random_normal([1, self.n_hidden]))
        return self.sess.run(self.reconstruction, feed_dict={self.hidden: hidden})

    def reconstruct(self, X):
        return self.sess.run(self.reconstruction, feed_dict={self.x: X})

    def getWeights(self):
        return self.sess.run(self.weights['w1'])

    def getBiases(self):
        return self.sess.run(self.weights['b1'])
    
def get_random_block_from_data(data, batch_size):
    start_index = np.random.randint(0, len(data) - batch_size)
    return data[start_index:(start_index + batch_size)]

class AdditiveGaussianNoiseAutoencoder(object):
    def __init__(self, n_input, n_hidden, transfer_function = tf.nn.softplus, optimizer = tf.train.AdamOptimizer(),
                 scale = 0.1):
        self.n_input = n_input
        self.n_hidden = n_hidden
        self.transfer = transfer_function
        self.scale = tf.placeholder(tf.float32)
        self.training_scale = scale
        network_weights = self._initialize_weights()
        self.weights = network_weights

        # model
        self.x = tf.placeholder(tf.float32, [None, self.n_input])
        self.hidden = self.transfer(tf.add(tf.matmul(self.x + scale * tf.random_normal((n_input,)),
                self.weights['w1']),
                self.weights['b1']))
        self.reconstruction = tf.add(tf.matmul(self.hidden, self.weights['w2']), self.weights['b2'])

        # cost
        self.cost = 0.5 * tf.reduce_sum(tf.pow(tf.subtract(self.reconstruction, self.x), 2.0))
        self.optimizer = optimizer.minimize(self.cost)

        init = tf.global_variables_initializer()
        self.sess = tf.Session()
        self.sess.run(init)

    def _initialize_weights(self):
        all_weights = dict()
        all_weights['w1'] = tf.get_variable("w1", shape=[self.n_input, self.n_hidden],
            initializer=tf.contrib.layers.xavier_initializer())
        all_weights['b1'] = tf.Variable(tf.zeros([self.n_hidden], dtype = tf.float32))
        all_weights['w2'] = tf.Variable(tf.zeros([self.n_hidden, self.n_input], dtype = tf.float32))
        all_weights['b2'] = tf.Variable(tf.zeros([self.n_input], dtype = tf.float32))
        return all_weights

    def partial_fit(self, X):
        cost, opt = self.sess.run((self.cost, self.optimizer), feed_dict = {self.x: X,
                                                                            self.scale: self.training_scale
                                                                            })
        return cost

    def calc_total_cost(self, X):
        return self.sess.run(self.cost, feed_dict = {self.x: X,
                                                     self.scale: self.training_scale
                                                     })

    def transform(self, X):
        return self.sess.run(self.hidden, feed_dict = {self.x: X,
                                                       self.scale: self.training_scale
                                                       })

    def generate(self, hidden=None):
        if hidden is None:
            hidden = self.sess.run(tf.random_normal([1, self.n_hidden]))
        return self.sess.run(self.reconstruction, feed_dict = {self.hidden: hidden})

    def reconstruct(self, X):
        return self.sess.run(self.reconstruction, feed_dict = {self.x: X,
                                                               self.scale: self.training_scale
                                                               })

    def getWeights(self):
        return self.sess.run(self.weights['w1'])

    def getBiases(self):
        return self.sess.run(self.weights['b1'])
    
class VariationalAutoencoder(object):

    def __init__(self, n_input, n_hidden, optimizer = tf.train.AdamOptimizer()):
        self.n_input = n_input
        self.n_hidden = n_hidden

        network_weights = self._initialize_weights()
        self.weights = network_weights

        # model
        self.x = tf.placeholder(tf.float32, [None, self.n_input])
        self.z_mean = tf.add(tf.matmul(self.x, self.weights['w1']), self.weights['b1'])
        self.z_log_sigma_sq = tf.add(tf.matmul(self.x, self.weights['log_sigma_w1']), self.weights['log_sigma_b1'])

        # sample from gaussian distribution
        eps = tf.random_normal(tf.stack([tf.shape(self.x)[0], self.n_hidden]), 0, 1, dtype = tf.float32)
        self.z = tf.add(self.z_mean, tf.multiply(tf.sqrt(tf.exp(self.z_log_sigma_sq)), eps))

        self.reconstruction = tf.add(tf.matmul(self.z, self.weights['w2']), self.weights['b2'])

        # cost
        reconstr_loss = 0.5 * tf.reduce_sum(tf.pow(tf.subtract(self.reconstruction, self.x), 2.0))
        latent_loss = -0.5 * tf.reduce_sum(1 + self.z_log_sigma_sq
                                           - tf.square(self.z_mean)
                                           - tf.exp(self.z_log_sigma_sq), 1)
        self.cost = tf.reduce_mean(reconstr_loss + latent_loss)
        self.optimizer = optimizer.minimize(self.cost)

        init = tf.global_variables_initializer()
        self.sess = tf.Session()
        self.sess.run(init)

    def _initialize_weights(self):
        all_weights = dict()
        all_weights['w1'] = tf.get_variable("w1", shape=[self.n_input, self.n_hidden],
            initializer=tf.contrib.layers.xavier_initializer())
        all_weights['log_sigma_w1'] = tf.get_variable("log_sigma_w1", shape=[self.n_input, self.n_hidden],
            initializer=tf.contrib.layers.xavier_initializer())
        all_weights['b1'] = tf.Variable(tf.zeros([self.n_hidden], dtype=tf.float32))
        all_weights['log_sigma_b1'] = tf.Variable(tf.zeros([self.n_hidden], dtype=tf.float32))
        all_weights['w2'] = tf.Variable(tf.zeros([self.n_hidden, self.n_input], dtype=tf.float32))
        all_weights['b2'] = tf.Variable(tf.zeros([self.n_input], dtype=tf.float32))
        return all_weights

    def partial_fit(self, X):
        cost, opt = self.sess.run((self.cost, self.optimizer), feed_dict={self.x: X})
        return cost

    def calc_total_cost(self, X):
        return self.sess.run(self.cost, feed_dict = {self.x: X})

    def transform(self, X):
        return self.sess.run(self.z_mean, feed_dict={self.x: X})

    def generate(self, hidden = None):
        if hidden is None:
            hidden = self.sess.run(tf.random_normal([1, self.n_hidden]))
        return self.sess.run(self.reconstruction, feed_dict={self.z: hidden})

    def reconstruct(self, X):
        return self.sess.run(self.reconstruction, feed_dict={self.x: X})

    def getWeights(self):
        return self.sess.run(self.weights['w1'])

    def getBiases(self):
        return self.sess.run(self.weights['b1'])
  
def autoencoder_dimension(X_train,dimension=100,epochs=20,batch_size = 128):
    tf.reset_default_graph() 
    n_samples=X_train.shape[0]
    n_input=X_train.shape[1]
    autoencoder = Autoencoder(
               n_input=n_input,
               n_hidden=dimension,
               transfer_function=tf.nn.softplus,
               optimizer=tf.train.AdamOptimizer(learning_rate=0.001))
    training_epochs = epochs
    batch_size = batch_size
    for epoch in range(training_epochs):
        avg_cost=0
        total_batch=int(n_samples/batch_size)
        for i in range(total_batch):
            batch_xs=get_random_block_from_data(X_train,batch_size)
            cost=autoencoder.partial_fit(batch_xs)
            avg_cost += cost / n_samples * batch_size
    X_train_transform=autoencoder.transform(X_train)
    return X_train_transform

def autoencoder_noise_dimension(X_train,dimension=100,epochs=20,batch_size = 128):
    tf.reset_default_graph()  
    n_samples=X_train.shape[0]
    n_input=X_train.shape[1]
    autoencoder = AdditiveGaussianNoiseAutoencoder(
               n_input=n_input,
               n_hidden=dimension,
               transfer_function=tf.nn.softplus,
               optimizer=tf.train.AdamOptimizer(learning_rate=0.001))
    training_epochs = epochs
    batch_size = batch_size
    for epoch in range(training_epochs):
        avg_cost=0
        total_batch=int(n_samples/batch_size)
        for i in range(total_batch):
            batch_xs=get_random_block_from_data(X_train,batch_size)
            cost=autoencoder.partial_fit(batch_xs)
            avg_cost += cost / n_samples * batch_size
    X_train_transform=autoencoder.transform(X_train)
    return X_train_transform  
 
   
def autoencoder_variation_dimension(X_train,dimension=100,epochs=20,batch_size = 128):
    tf.reset_default_graph()  
    n_samples=X_train.shape[0]
    n_input=X_train.shape[1]
    autoencoder = VariationalAutoencoder(
               n_input=n_input,
               n_hidden=dimension,
               optimizer=tf.train.AdamOptimizer(learning_rate=0.001))
    training_epochs = epochs
    batch_size = batch_size
    for epoch in range(training_epochs):
        avg_cost=0
        total_batch=int(n_samples/batch_size)
        for i in range(total_batch):
            batch_xs=get_random_block_from_data(X_train,batch_size)
            cost=autoencoder.partial_fit(batch_xs)
            avg_cost += cost / n_samples * batch_size
    X_train_transform=autoencoder.transform(X_train)
    return X_train_transform

# main code

if(option == "1"):
    # Feature selection T-SNE
    tsneselection(selected_number,feature_o)
elif(option == "2"):
    # Feature selection PCA
    pcaselection(selected_number,feature_o)
elif(option == "3"):
    # Feature selection KernelPCA
    kernelPCA(selected_number,feature_o)
elif(option == "4"):
    # Feature selection LLE
    lle(selected_number,feature_o)
elif(option == "5"):
    # Feature selection SVD
    svd(selected_number,feature_o)
elif(option == "6"):
    # Feature selection NMF
    nmf(selected_number,feature_o)
elif(option == "7"):
    # Feature selection MDS
    mds(selected_number,feature_o)
elif(option == "8"):
    # Feature selection ICA
    X=np.array(feature_o)
    ICA(X,n_components=selected_number)
elif(option == "9"):
    # Feature selection FA
    X=np.array(feature_o)
    FA(X,n_components=selected_number)
elif(option == "10"):
    # Feature selection Agglomerate feature
    Aggselection(selected_number,feature_o)
elif(option == "11"):
    # Feature selection Gaussian Random Projection
    Gaussianselection(selected_number,feature_o)
elif(option == "12"):
    # Feature selection Sparse Random Projection
    Sparseselection(selected_number,feature_o) 
elif(option == "13"):
     vector=autoencoder_dimension(feature_o,dimension=selected_number)
     vector.to_csv("Autoencoder_out.csv")
elif(option == "14"):
    vector=autoencoder_noise_dimension(feature_o,dimension=selected_number)
    vector.to_csv("Autoencoder_noise_out.csv")
elif(option == "15"):
    vector=autoencoder_variation_dimension(feature_o,dimension=selected_number)
    vector.to_csv("Autoencoder_variation_out.csv")
else:
    print("Invalid method number. Please check the method table!")



