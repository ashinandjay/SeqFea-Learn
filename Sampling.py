import sys,os
import re
import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from imblearn.over_sampling import SMOTE, ADASYN
from imblearn.under_sampling import RandomUnderSampler
from imblearn.under_sampling import NeighbourhoodCleaningRule
from imblearn.over_sampling import RandomOverSampler

Script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

filename = sys.argv[1]

feature = pd.read_csv(filename,header=None,index_col=0,float_precision='round_trip')

# read feature labels
labelfile = sys.argv[2]
#labelfile = 'S_label.txt'

file = open(labelfile,'r')
index=[]
label=[]
for line in file.readlines():
    if(re.match(">",line)):
        array=re.split('>',line)
        index.append(array[1])
    else:
        label.append(int(line[0]))
data={'Index':index,'Label':label}
label = pd.DataFrame(data)
label = label.set_index(['Index'])

#######
X_= np.array(feature)[:,1:]

X = scale(X_)

y= np.array(label.values.ravel())

# choose the method
option = sys.argv[1]
# the input sequence
file = sys.argv[2]

if(option == "1"):
    #Random over sampling method
    ros = RandomOverSampler(random_state=0)
    X_resampled, y_resampled = ros.fit_resample(X, y)
    csv_X=pd.DataFrame(data=X_resampled)
    csv_y=pd.DataFrame(data=y_resampled)
    csv_X.to_csv('ros_feature.csv',header=False,index=False)
    csv_y.to_csv('ros_label.csv',header=False,index=False)
    
if(option == "2"):
    #SMOTE method
    X_resampled, y_resampled = SMOTE().fit_resample(X, y) 
    csv_X=pd.DataFrame(data=X_resampled)
    csv_y=pd.DataFrame(data=y_resampled)
    csv_X.to_csv('ros_feature.csv',header=False,index=False)
    csv_y.to_csv('ros_label.csv',header=False,index=False)
    
if(option == "3"):
    #ADASYN method
    X_resampled, y_resampled = ADASYN().fit_resample(X, y) 
    csv_X=pd.DataFrame(data=X_resampled)
    csv_y=pd.DataFrame(data=y_resampled)
    csv_X.to_csv('ros_feature.csv',header=False,index=False)
    csv_y.to_csv('ros_label.csv',header=False,index=False)

if(option == "4"):
    #Random under sampling method
    rus = RandomUnderSampler(random_state=0)
    X_resampled, y_resampled = rus.fit_resample(X, y)
    csv_X=pd.DataFrame(data=X_resampled)
    csv_y=pd.DataFrame(data=y_resampled)
    csv_X.to_csv('ros_feature.csv',header=False,index=False)
    csv_y.to_csv('ros_label.csv',header=False,index=False)
    
if(option == "5"):
    #Neighbourhood cleaning rule method
    ncr = NeighbourhoodCleaningRule()
    X_resampled, y_resampled = ncr.fit_resample(X, y)
    csv_X=pd.DataFrame(data=X_resampled)
    csv_y=pd.DataFrame(data=y_resampled)
    csv_X.to_csv('ros_feature.csv',header=False,index=False)
    csv_y.to_csv('ros_label.csv',header=False,index=False)