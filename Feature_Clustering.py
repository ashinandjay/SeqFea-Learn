# This code is used for clustering methods, 
# The output file will be generated in same folder of this code
# -------------------------------------------------------------

# load packages:
import sys,os
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import MeanShift
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
from sklearn.cluster import SpectralClustering
from sklearn import mixture
# find the path
Script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

# choose the method
option = sys.argv[1]

# read feature vectors
feature = sys.argv[2]
feature_o = pd.read_csv(feature,delim_whitespace=True,index_col=0,float_precision='round_trip')
k=sys.argv[3]

# Feature clustering Kmeans
def kmeansselection(k,X):
    X=np.array(X)
    kmeansresult=KMeans(n_clusters = k).fit_transform(X)
    np.savetxt("KMeans_out.csv", kmeansresult, delimiter=",")    
    return None 

# Feature clustering Affinty Propagation
def clustering_AP(X):
    X=np.array(X)
    APresult=AffinityPropagation().fit_predict(X)
    np.savetxt("AP.csv", APresult, delimiter=",")    
    return None

# Feature clustering MeanShift
def clustering_MS(X):
    X=np.array(X)
    MSresult=MeanShift().fit_predict(X)
    np.savetxt("MS.csv", MSresult, delimiter=",")    
    return None

# Feature clustering SpectralClustering
def clustering_DBSCAN(X):
    X=np.array(X)
    DBSCANresult=DBSCAN().fit_predict(X)
    np.savetxt("DBSCAN.csv", DBSCANresult, delimiter=",")    
    return None

# Feature clustering OPTICS
def clustering_OPTICS(X):
    X=np.array(X)
    OPTICSresult=OPTICS().fit_predict(X)
    np.savetxt("OPTICS.csv", OPTICSresult, delimiter=",")    
    return None

# Spectral Clustering
def Spectral(k,X):
    X=np.array(X)
    spectralsresult=SpectralClustering(n_components=k).fit_predict(X)
    np.savetxt("SpectralClustering_out.csv", spectralsresult, delimiter=",")

# GaussianMixture Clustering
def GMClustering(k,X):
    X=np.array(X)
    GMresulte=mixture.GaussianMixture(n_components=k).fit_predict(X)
    np.savetxt("GaussianMixture_out.csv", GMresulte, delimiter=",")
    
if(option == "1"):
   # Kmeans cluster
   kmeansselection(k,feature_o)
elif(option == "2"):
    # AP cluster
    clustering_AP(feature_o)
elif(option == "3"):
    # MS cluster
    clustering_MS(feature_o)
elif(option == "4"):
    # DBSCAN cluster
    clustering_DBSCAN(feature_o)
elif(option == "5"):
    # OPTICS cluster
    clustering_OPTICS(feature_o)
elif(option == "6"):
    # Spectral cluster
    Spectral(k,feature_o)
elif(option == "7"):
    # GMC cluster
    GMClustering(k,feature_o)
else:
    print("Invalid method number. Please check the method table!")
