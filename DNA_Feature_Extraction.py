import os
import sys,re
import pickle
from collections import Counter
from functools import reduce
import itertools
import pandas as pd
import numpy as np
import platform

ALPHABET='ACGT'

def readDNAFasta(file):
	with open(file) as f:
		records = f.read()
	if re.search('>', records) == None:
		print('Error,the input DNA sequence must be fasta format.')
		sys.exit(1)
	records = records.split('>')[1:]
	myFasta = []
	for fasta in records:
		array = fasta.split('\n')
		name, sequence = array[0].split()[0], re.sub('[^ACGT-]', '-', ''.join(array[1:]).upper())
		myFasta.append([name, sequence])
	return myFasta

def frequency(t1_str, t2_str):

    i, j, tar_count = 0, 0, 0
    len_tol_str = len(t1_str)
    len_tar_str = len(t2_str)
    while i < len_tol_str and j < len_tar_str:
        if t1_str[i] == t2_str[j]:
            i += 1
            j += 1
            if j >= len_tar_str:
                tar_count += 1
                i = i - j + 1
                j = 0
        else:
            i = i - j + 1
            j = 0
    return tar_count

def generate_list(k, alphabet):
    ACGT_list=["".join(e) for e in itertools.product(alphabet, repeat=k)]
    return ACGT_list

def generate_property(raw_property, new_property):
    if new_property is None or len(new_property) == 0:
       return raw_property
    for key in list(raw_property.keys()):
        raw_property[key].extend(new_property[key])
    return raw_property

def obtain_property(k, property_list):
    property_value = {}
    if 0 == len(property_list):
        for nucleotide in generate_list(k, ALPHABET):
            property_value[nucleotide] = []
        return property_value

    nucleotide_value = get_property_dic(k)
    for nucleotide in generate_list(k, ALPHABET):
        if nucleotide not in property_value:
           property_value[nucleotide] = []
        for e in nucleotide_value[nucleotide]:
            if e[0] in property_list:
               property_value[nucleotide].append(e[1])
    return property_value

def get_property_psednc():

    raw_property = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11],
                   'AC': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04],
                   'AG': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62],
                   'AT': [1.07, 0.22, 0.62, -1.02, 2.51, 1.17],
                   'CA': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25],
                   'CC': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24],
                   'CG': [-1.66, -1.22, -0.44, -0.82, -0.29, -1.39],
                   'CT': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62],
                   'GA': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71],
                   'GC': [-0.08, 0.22, 1.33, -0.35, 0.65, 1.59],
                   'GG': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24],
                   'GT': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04],
                   'TA': [-1.23, -2.37, -0.44, -2.24, -1.51, -1.39],
                   'TC': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71],
                   'TG': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25],
                   'TT': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11]}
    extra_phyche_index={}
    property_value = generate_property(raw_property, extra_phyche_index)
    return property_value

def get_physicochemical_list():
    diphyche_list = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content',
                  'A-philicity', 'Propeller twist', 'Duplex stability:(freeenergy)',
                  'Duplex tability(disruptenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
                  'Stabilising energy of Z-DNA', 'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH',
                  'Breslauer_dS', 'Electron_interaction', 'Hartman_trans_free_energy', 'Helix-Coil_transition',
                  'Ivanov_BA_transition', 'Lisser_BZ_transition', 'Polar_interaction', 'SantaLucia_dG',
                  'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Stacking_energy',
                  'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS', 'Watson-Crick_interaction', 'Twist', 'Tilt',
                  'Roll', 'Shift', 'Slide', 'Rise']
    triphyche_list = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
                  'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
                  'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']
    return diphyche_list,triphyche_list
def generate_property_value(k):
    
    diphyche_list,triphyche_list=get_physicochemical_list()
    if k==2:
       phyche_index=diphyche_list
    if k==3:
       phyche_index=triphyche_list
    extend_index={}
    phyche_index_dic=generate_property(obtain_property(k, phyche_index),extend_index)
    return phyche_index_dic

def get_property_dic(k):

    if 2 == k:
        file_path = os.path.split(os.path.realpath(__file__))[0] + r'\data\phy1.data' if platform.system() == 'Windows' else  os.path.split(os.path.realpath(__file__))[0] + '/data/phy1.data'
    elif 3 == k:
        file_path =os.path.split(os.path.realpath(__file__))[0] + r'\data\phy2.data' if platform.system() == 'Windows' else  os.path.split(os.path.realpath(__file__))[0] + '/data/phy1.data'
    else:
        sys.stderr.write("The k can just be 2 or 3.")
        sys.exit(0)
    try:
        with open(file_path, 'rb') as f:
            property_factor_dic = pickle.load(f)
    except:
        with open(file_path, 'r') as f:
            property_factor_dic = pickle.load(f)
    return property_factor_dic

def parallel_cor_function(nucleotide1, nucleotide2, phyche_index):
    temp_sum = 0.0
    phyche_index_values = list(phyche_index.values())
    len_phyche_index = len(phyche_index_values[0])
    for u in range(len_phyche_index):
        temp_sum += pow(float(phyche_index[nucleotide1][u]) - float(phyche_index[nucleotide2][u]), 2)
    parallel_value=temp_sum / len_phyche_index   
    return parallel_value

def get_parallel_factor(k, lamada, sequence, phyche_value):

    theta = []
    l = len(sequence)
    for i in range(1, lamada + 1):
        temp_sum = 0.0
        for j in range(0, l - k - i + 1):
            nucleotide1 = sequence[j: j+k]
            nucleotide2 = sequence[j+i: j+i+k]
            temp_sum += parallel_cor_function(nucleotide1, nucleotide2, phyche_value)
        theta.append(temp_sum / (l - k - i + 1))
    return theta

def get_parallel_factor_psednc(lamada, sequence, phyche_value):

    theta = []
    l = len(sequence)

    for i in range(1, lamada + 1):
        temp_sum = 0.0
        for j in range(0, l - 1 - lamada):
            nucleotide1 = sequence[j] + sequence[j + 1]
            nucleotide2 = sequence[j + i] + sequence[j + i + 1]
            temp_sum += parallel_cor_function(nucleotide1, nucleotide2, phyche_value)
        theta.append(temp_sum / (l - i - 1))
    return theta

def make_pseknc_vector(sequence_list, lamada, w, k, phyche_value):

    kmer = generate_list(k, ALPHABET)
    header = ['#']
    for f in range((16+lamada)):
        header.append('pseknc.'+str(f))
    vector=[]
    vector.append(header)
    for sequence_ in sequence_list:
        name,sequence=sequence_[0],sequence_[1]
        if len(sequence) < k or lamada + k > len(sequence):
            error_info = "Error, the sequence length must be larger than " + str(lamada + k)
            sys.stderr.write(error_info)
        fre_list = [frequency(sequence, str(key)) for key in kmer]
        fre_sum = float(sum(fre_list))
        fre_list = [e / fre_sum for e in fre_list]
        theta_list = get_parallel_factor(k, lamada, sequence, phyche_value)
        theta_sum = sum(theta_list)
        denominator = 1 + w * theta_sum
        
        temp_vec = [round(f / denominator, 3) for f in fre_list]      
        for theta in theta_list:
            temp_vec.append(round(w * theta / denominator, 4))
        sample=[name]
        sample=sample+temp_vec
        vector.append(sample)
    return vector

def make_old_pseknc_vector(sequence_list, lamada, w, k, phyche_value):

    kmer = generate_list(k, ALPHABET)
    header=['#']
    for f in range((4**k+lamada)):
        header.append('pseknc.'+str(f))
    vector=[]
    vector.append(header)
    for sequence_ in sequence_list:
        name,sequence=sequence_[0],sequence_[1]
        if len(sequence) < k or lamada + k > len(sequence):
            error_info = "error, the sequence length must be larger than " + str(lamada + k)
            sys.stderr.write(error_info)
            sys.exit(0)
        fre_list = [frequency(sequence, str(key)) for key in kmer]
        fre_sum = float(sum(fre_list))
        fre_list = [e / fre_sum for e in fre_list]
        theta_list = get_parallel_factor_psednc(lamada, sequence, phyche_value)
        theta_sum = sum(theta_list)
        denominator = 1 + w * theta_sum
        temp_vec = [round(f / denominator, 3) for f in fre_list]
        for theta in theta_list:
            temp_vec.append(round(w * theta / denominator, 4))
        sample=[name]
        sample=sample+temp_vec
        vector.append(sample)
    return vector

def make_ac_vector(sequence_list, lag, phyche_value, k):
    phyche_values = list(phyche_value.values())
    len_phyche_value = len(phyche_values[0])
    vector = []
    ac_vector=[]
    header=['#']
    for f in range(lag*len_phyche_value):
        header.append('AC.'+str(f))
    vector.append(header)
    ac_vector.append(header)
    for sequence_ in sequence_list:
        name,sequence=sequence_[0],sequence_[1]
        len_seq = len(sequence)
        each_vec = []
        ac_vec=[]
        for temp_lag in range(1, lag + 1):
            for j in range(len_phyche_value):

                ave_phyche_value = 0.0
                for i in range(len_seq - temp_lag - k + 1):
                    nucleotide = sequence[i: i + k]
                    ave_phyche_value += float(phyche_value[nucleotide][j])
                ave_phyche_value /= len_seq
                # Calculate the vector.
                temp_sum = 0.0
                for i in range(len_seq - temp_lag - k + 1):
                    nucleotide1 = sequence[i: i + k]
                    nucleotide2 = sequence[i + temp_lag: i + temp_lag + k]
                    temp_sum += (float(phyche_value[nucleotide1][j]) - ave_phyche_value) * (
                        float(phyche_value[nucleotide2][j]))
                each_vec.append(round(temp_sum / (len_seq - temp_lag - k + 1), 3))
        sample=[name]
        sample=sample+each_vec
        vector.append(sample)
        ac_vec=ac_vec+each_vec
        ac_vector.append(ac_vec)
    return vector,ac_vector

def make_cc_vector(sequence_list, lag, phyche_value, k):
    phyche_values = list(phyche_value.values())
    len_phyche_value = len(phyche_values[0])
    vector = []
    cc_vector=[]
    header=['#']
    for f in range(lag*len_phyche_value*(len_phyche_value-1)):
        header.append('CC.'+str(f))
    vector.append(header)
    cc_vector.append(header[1:])
    for sequence_ in sequence_list:
        name,sequence=sequence_[0],sequence_[1]
        len_seq = len(sequence)
        each_vec = []
        cc_vec=[]
        for temp_lag in range(1, lag + 1):
            for i1 in range(len_phyche_value):
                for i2 in range(len_phyche_value):
                    if i1 != i2:
                        # Calculate average phyche_value for a nucleotide.
                        ave_phyche_value1 = 0.0
                        ave_phyche_value2 = 0.0
                        for j in range(len_seq - temp_lag - k + 1):
                            nucleotide = sequence[j: j + k]
                            ave_phyche_value1 += float(phyche_value[nucleotide][i1])
                            ave_phyche_value2 += float(phyche_value[nucleotide][i2])
                        ave_phyche_value1 /= len_seq
                        ave_phyche_value2 /= len_seq
                        # Calculate the vector.
                        temp_sum = 0.0
                        for j in range(len_seq - temp_lag - k + 1):
                            nucleotide1 = sequence[j: j + k]
                            nucleotide2 = sequence[j + temp_lag: j + temp_lag + k]
                            temp_sum += (float(phyche_value[nucleotide1][i1]) - ave_phyche_value1) * \
                                        (float(phyche_value[nucleotide2][i2]) - ave_phyche_value2)
                        each_vec.append(round(temp_sum / (len_seq - temp_lag - k + 1), 3))
        sample=[name]
        sample=sample+each_vec
        vector.append(sample)
        cc_vec=cc_vec+each_vec
        cc_vector.append(cc_vec)
    return vector,cc_vector

def acc_property(k):

    property_value = generate_property_value(k) 
    return property_value

def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer

def RC(kmer):
    myDict = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    return ''.join([myDict[nc] for nc in kmer[::-1]])

def generateRCKmer(kmerList):
    rckmerList = set()
    myDict = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    for kmer in kmerList:
        rckmerList.add(sorted([kmer, ''.join([myDict[nc] for nc in kmer[::-1]])])[0])
    return sorted(rckmerList)

#########################################################################
def Kmer(input_data, k=2, normalize=True):
    
    fastas=readDNAFasta(input_data)
    vector = []
    header = ['#']
    if k < 1:
        print('error, the k must be positive integer.')
        return 0
    for kmer in itertools.product(ALPHABET, repeat=k):
        header.append(''.join(kmer))
    vector.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        kmers = kmerArray(sequence, k)
        count = Counter()
        count.update(kmers)
        if normalize == True:
           for key in count:
               count[key] = count[key] / len(kmers)
        code = [name]
        for j in range(1, len(header)):
            if header[j] in count:
               code.append(count[header[j]])
            else:
                code.append(0)
        vector.append(code)
    return vector

#vector=Kmer('DNA_data.txt',k=4)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('Kmer_out.csv',header=False,index=False)

def RCKmer(input_data, k=2, normalize=True):
    
    fastas=readDNAFasta(input_data)
    vector = []
    header = ['#']
    if k < 1:
        print('error, the k must be positive integer.')
        return 0
    tmpHeader = []
    for kmer in itertools.product(ALPHABET, repeat=k):
        tmpHeader.append(''.join(kmer))
    header = header + generateRCKmer(tmpHeader)
    myDict = {}
    for kmer in header[2:]:
        rckmer = RC(kmer)
        if kmer != rckmer:
           myDict[rckmer] = kmer
    vector.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        kmers = kmerArray(sequence, k)
        for j in range(len(kmers)):
            if kmers[j] in myDict:
               kmers[j] = myDict[kmers[j]]
        count = Counter()
        count.update(kmers)
        if normalize == True:
           for key in count:
               count[key] = count[key] / len(kmers)
        code = [name]
        for j in range(1, len(header)):
            if header[j] in count:
                code.append(count[header[j]])
            else:
                code.append(0)
        vector.append(code)
    return vector

#vector=RCKmer('DNA_data.txt',k=2)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('RCKmer_out.csv',header=False,index=False)

def Psednc(input_data,lamada=10, w=0.05, k = 2):
    #Psednc
    phyche_value = get_property_psednc()    
    fastas=readDNAFasta(input_data)
    vector = make_pseknc_vector(fastas, lamada, w, k, phyche_value)
    return vector

#vector=Psednc('DNA_data.txt',lamada=10)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('Psednc_out.csv',header=False,index=False)

def Pseknc(input_data, k=3, lamada=10, w=0.05):
    #Pseknc
    phyche_value = get_property_psednc()        
    fastas=readDNAFasta(input_data)  
    vector=make_old_pseknc_vector(fastas, lamada, w, k, phyche_value)    
    return vector

#vector=Pseknc('DNA_data.txt',lamada=10)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('Pseknc_out.csv',header=False,index=False)

def DAC(input_data,k=2,lag=5):
    #DAC
    phyche_value = acc_property(k)   
    fastas=readDNAFasta(input_data)   
    vector,_=make_ac_vector(fastas, lag, phyche_value, k)    
    return vector
#vector=DAC('DNA_data.txt',lag=5)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('DAC_out.csv',header=False,index=False)

def DCC(input_data,k=2,lag=2):
    #DCC
    phyche_value = acc_property(k)   
    fastas=readDNAFasta(input_data)   
    vector,_=make_cc_vector(fastas, lag, phyche_value, k)    
    return vector

#vector=DCC('DNA_data.txt',lag=1)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('DCC_out.csv',header=False,index=False)

def DACC(input_data, k=2,lag=1):
    #DACC
    phyche_value = acc_property(k)
    fastas=readDNAFasta(input_data)
    vector1,ac_vector=make_ac_vector(fastas, lag, phyche_value, k)
    vector2,cc_vector=make_cc_vector(fastas, lag, phyche_value, k)    
    zipped = list(zip(vector1,cc_vector))   
    vector = [reduce(lambda x, y: x + y, e) for e in zipped]
    return vector

#vector=DACC('DNA_data.txt',lag=1)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('DACC_out.csv',header=False,index=False)

def TAC(input_data,k=3,lag=1):
    #TAC
    phyche_value = acc_property(k) 
    fastas=readDNAFasta(input_data)    
    vector,_=make_ac_vector(fastas, lag, phyche_value, k)    
    return vector

#vector=TAC('DNA_data.txt',lag=1)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('TAC_out.csv',header=False,index=False)

def TCC(input_data,k=3,lag=1):
    #TCC
    phyche_value = acc_property(k)   
    fastas=readDNAFasta(input_data)    
    vector,_=make_cc_vector(fastas, lag, phyche_value, k)   
    return vector

#vector=TCC('DNA_data.txt',lag=1)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('TCC_out.csv',header=False,index=False)

def TACC(input_data, k=3,lag=1):
    #TACC
    phyche_value = acc_property(k)
    fastas=readDNAFasta(input_data)
    vector1,ac_vector=make_ac_vector(fastas, lag, phyche_value, k)
    vector2,cc_vector=make_cc_vector(fastas, lag, phyche_value, k)    
    zipped = list(zip(vector1,cc_vector))   
    vector = [reduce(lambda x, y: x + y, e) for e in zipped]
    return vector

#vector=TACC('DNA_data.txt',lag=1)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('TACC_out.csv',header=False,index=False)

def NAC(input_data):
    fastas=readDNAFasta(input_data)
    vector = []
    header = ['#']
    for f in range(4):
        header.append('NAC.'+str(f))
    vector.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = [name]
        for na in ALPHABET:
            code.append(count[na])
        vector.append(code)
    return vector

#vector=NAC('DNA_data.txt')
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('NAC_out.csv',header=False,index=False)

def DNC(input_data):
    fastas=readDNAFasta(input_data)
    vector = []
    header = ['#']
    for f in range(16):
        header.append('DNC.'+str(f))
    vector.append(header)
    AADict = {}
    for i in range(len(ALPHABET)):
        AADict[ALPHABET[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        tmpCode = [0] * 16
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        vector.append(code)
    return vector

#vector=DNC('DNA_data.txt')
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('DNC_out.csv',header=False,index=False)
 
def TNC(input_data):
    fastas=readDNAFasta(input_data)
    vector = []
    header = ['#']
    for f in range(64):
        header.append('TNC.'+str(f))
    vector.append(header)
    AADict = {}
    for i in range(len(ALPHABET)):
        AADict[ALPHABET[i]] = i
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        tmpCode = [0] * 64
        for j in range(len(sequence) - 3 + 1):
            tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j+1]]*4 + AADict[sequence[j+2]]] = tmpCode[AADict[sequence[j]] * 16 + AADict[sequence[j+1]]*4 + AADict[sequence[j+2]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        vector.append(code)
    return vector

#vector=TNC('DNA_data.txt')
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('TNC_out.csv',header=False,index=False)  

def zCurve(x):
    t=[]
    A = x.count('A'); C = x.count('C'); G = x.count('G');TU=x.count('U')
    x_ = (A + G) - (C + TU)
    y_ = (A + C) - (G + TU)
    z_ = (A + TU) - (C + G)
            # print(x_, end=','); print(y_, end=','); print(z_, end=',')
    t.append(x_); t.append(y_); t.append(z_)

    return t
def zCurve_vector(input_data):   
    fastas=readDNAFasta(input_data)
    header=['#']
    for f in range (3):
        header.append('zCurve.'+str(f))           
    vector=[] 
    vector.append(header) 
    sample=[]
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        sample = [name]
        each_vec=zCurve(sequence)
        sample=sample+each_vec
        vector.append(sample)
    return vector
#vector=zCurve_vector('DNA_data.txt')
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('zCurve_out.csv',header=False,index=False)  
  
def kmers(seq, k):
    v = []
    for i in range(len(seq) - k + 1):
        v.append(seq[i:i + k])
    return v
def MonoKGap(x, g):  # 1___1
    t=[]
    m = list(itertools.product(ALPHABET, repeat=2))
    L_sequence=(len(x)-g-1)
    for i in range(1, g + 1, 1):
        V = kmers(x, i + 2)
        for gGap in m:
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[-1] == gGap[1]:
                    C += 1
            t.append(C/L_sequence)    
    return t

def MonoKGap_vector(input_data,g):   
    fastas=readDNAFasta(input_data)
    vector=[] 
    header=['#']
    for f in range((g)*16):
        header.append('Mono.'+str(f))
    vector.append(header)   
    sample=[]
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        sample = [name]
        each_vec=MonoKGap(sequence,g)
        sample=sample+each_vec
        vector.append(sample)
    return vector
#vector=MonoKGap_vector('DNA_data.txt',g=2)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('MonoKGap_out.csv',header=False,index=False)   

def MonoDiKGap(x, g):     
    t=[]
    L_sequence=(len(x)-g-2)
    m = list(itertools.product(ALPHABET, repeat=3))
    for i in range(1, g + 1, 1):
        V = kmers(x, i + 3)
        for gGap in m:
            C = 0
            for v in V:
                if v[0] == gGap[0] and v[-2] == gGap[1] and v[-1] == gGap[2]:
                    C += 1
            t.append(C/L_sequence) 
    return t
def MonoDiKGap_vector(input_data,g):   
    fastas=readDNAFasta(input_data)
    vector=[] 
    header=['#']
    for f in range((g)*32):
        header.append('MonoDi.'+str(f))
    vector.append(header)
    sample=[]
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        sample = [name]
        each_vec=MonoDiKGap(sequence,g)
        sample=sample+each_vec
        vector.append(sample)
    return vector
#vector=MonoDiKGap_vector('DNA_data.txt',g=2)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('MonoDiKGap_out.csv',header=False,index=False) 
    
def binary(input_data):
    fastas=readDNAFasta(input_data)
    vector = []
    header = ['#']
    for i in range(1, len(fastas[0][1]) * 4 + 1):
        header.append('BINARY.F'+str(i))
    vector.append(header)
    for i in fastas:
        name, sequence = i[0], i[1]
        code = [name]
        for aa in sequence:
            if aa == '-':
                code = code + [0, 0, 0, 0]
                continue
            for aa1 in ALPHABET:
                tag = 1 if aa == aa1 else 0
                code.append(tag)
        vector.append(code)
    return vector
#vector=binary('DNA_data.txt')
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('binary_out.csv',header=False,index=False) 
    
def CalculateMatrix_Di(data, order):
    matrix = np.zeros((len(data[0]) - 1, 16))
    for i in range(len(data[0]) - 1): # position
        for j in range(len(data)):
            if re.search('-', data[j][i:i+2]):
                pass
            else:
                matrix[i][order[data[j][i:i+2]]] += 1
    return matrix

def PSDNP(input_data,label):
    
    fastas=readDNAFasta(input_data)
    for i in fastas:
        if re.search('[^ACGT-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT[U]" are allowed by this encoding scheme.')
            return 0

    encodings = []
    header = ['#']
    for pos in range(len(fastas[0][1])-1):
        header.append('PSDNP.%d' %(pos+1))
    encodings.append(header)

    positive = []
    negative = []
    positive_key = []
    negative_key = []
    for i,sequence in enumerate(fastas):
        if int(label[i]) == 1:
           positive.append(sequence[1])
           positive_key.append(sequence[0])
        else:
           negative.append(sequence[1])
           negative_key.append(sequence[0])

    nucleotides = ['A', 'C', 'G', 'T']
    trinucleotides = [n1 + n2 for n1 in nucleotides for n2 in nucleotides]
    order = {}
    for i in range(len(trinucleotides)):
        order[trinucleotides[i]] = i

    matrix_po = CalculateMatrix_Di(positive, order)
    matrix_ne = CalculateMatrix_Di(negative, order)

    positive_number = len(positive)
    negative_number = len(negative)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = [name]
        for j in range(len(sequence) - 1):
            if re.search('-', sequence[j: j+2]):
               code.append(0)
            else:
                p_num, n_num = positive_number, negative_number
                po_number = matrix_po[j][order[sequence[j: j+2]]]
                if i[0] in positive_key and po_number > 0:
                   po_number -= 1
                   p_num -= 1
                ne_number = matrix_ne[j][order[sequence[j: j+2]]]
                if i[0] in negative_key and ne_number > 0:
                   ne_number -= 1
                   n_num -= 1
                code.append(po_number/p_num - ne_number/n_num)
                    # print(sequence[j: j+3], order[sequence[j: j+3]], po_number, p_num, ne_number, n_num)
        encodings.append(code)
    return encodings,positive_key, negative_key

def CalculateMatrix(data, order):
    matrix = np.zeros((len(data[0]) - 2, 64))
    for i in range(len(data[0]) - 2): # position
        for j in range(len(data)):
            if re.search('-', data[j][i:i+3]):
                pass
            else:
                matrix[i][order[data[j][i:i+3]]] += 1
    return matrix

def PSTNP(input_data,label):
    fastas=readDNAFasta(input_data)
    for i in fastas:
        if re.search('[^ACGU-]', i[1]):
            print('Error: illegal character included in the fasta sequences, only the "ACGT[U]" are allowed by this encoding scheme.')
            return 0

    encodings = []
    header = ['#']
    for pos in range(len(fastas[0][1])-2):
        header.append('PSTNP.%d' %(pos+1))
    encodings.append(header)
    positive = []
    negative = []
    positive_key = []
    negative_key = []
    for i,sequence in enumerate(fastas):
        if int(label[i]) == 1:
           positive.append(sequence[1])
           positive_key.append(sequence[0])
        else:
           negative.append(sequence[1])
           negative_key.append(sequence[0])

    nucleotides = ['A', 'C', 'G', 'T']
    trinucleotides = [n1 + n2 + n3 for n1 in nucleotides for n2 in nucleotides for n3 in nucleotides]
    order = {}
    for i in range(len(trinucleotides)):
        order[trinucleotides[i]] = i

    matrix_po = CalculateMatrix(positive, order)
    matrix_ne = CalculateMatrix(negative, order)

    positive_number = len(positive)
    negative_number = len(negative)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = [name]
        for j in range(len(sequence) - 2):
            if re.search('-', sequence[j: j+3]):
               code.append(0)
            else:
                p_num, n_num = positive_number, negative_number
                po_number = matrix_po[j][order[sequence[j: j+3]]]
                ne_number = matrix_ne[j][order[sequence[j: j+3]]]
                code.append(po_number/p_num - ne_number/n_num)
        encodings.append(code)
    return encodings,positive_key, negative_key
# find the path
Script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))

# choose the method
option = sys.argv[1]

# the input sequence
fastas = sys.argv[2]

 
if(option == "1"):
    #kmer method
    vector=Kmer(fastas, k=3)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('Kmer_out.csv',header=False,index=False)
elif(option == "2"):
    #RCKmer method
    vector=RCKmer(fastas, k=2)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('RCKmer_out.csv',header=False,index=False)
elif(option == "3"):
    #Psednc method
    vector=Psednc(fastas,lamada=10)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('Psednc_out.csv',header=False,index=False)
elif(option == "4"):
    #Pseknc method
    vector=Pseknc(fastas,lamada=10)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('Pseknc_out.csv',header=False,index=False)
elif(option == "5"):
    #DAC method
    vector=DAC(fastas,lag=5)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DAC_out.csv',header=False,index=False)
elif(option == "6"):
    #DCC method
    vector=DCC(fastas,lag=1)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DCC_out.csv',header=False,index=False)
elif(option == "7"):
    #DACC method
    vector=DACC(fastas,lag=1)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DACC_out.csv',header=False,index=False)
elif(option == "8"):
    #TAC method
    vector=TAC(fastas,lag=1)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('TAC_out.csv',header=False,index=False)
elif(option == "9"):
    #TCC method
    vector=TCC(fastas,lag=1)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('TCC_out.csv',header=False,index=False)
elif(option == "10"):
    #TACC method
    vector=TACC(fastas,lag=1)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('TACC_out.csv',header=False,index=False)
elif(option == "11"):
    #NAC method
    vector=NAC(fastas)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('NAC_out.csv',header=False,index=False)
elif(option == "12"):
    #DNC method
    vector=DNC(fastas)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DNC_out.csv',header=False,index=False)
elif(option == "13"):
    #TNC method
    vector=TNC(fastas)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('TNC_out.csv',header=False,index=False)
elif(option == "14"):
    #zCurve method
    vector=zCurve(fastas)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('zCurve_out.csv',header=False,index=False)
elif(option == "15"):
    #monoMonoKGap method
    vector=MonoKGap_vector(fastas, g=2)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('MonoKGap_out.csv',header=False,index=False)
elif(option == "16"):
    #monoDiKGap method
    vector=MonoDiKGap_vector(fastas, g=2)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('MonoDiKGap_out.csv',header=False,index=False)
elif(option == "17"):
    #binary method
    vector=binary(fastas)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('binary_out.csv',header=False,index=False)
    
elif(option == "18"):
    #PSTNP method
    vector,A1,A2=PSDNP(fastas,label)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('PSDNP_out.csv',header=False,index=False)
    
elif(option == "19"):
    #PSDNP method
    vector,A1,A2=PSTNP(fastas,label)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('PSTNP_out.csv',header=False,index=False)
else:
    print("Invalid method number. Please check the method table!")



