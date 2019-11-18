import sys,re,os
from functools import reduce
from collections import Counter
import pandas as pd
import numpy as np
import itertools
import platform
from math import sqrt
from math import pow

ALPHABET='ACGU'

def readRNAFasta(file):
	with open(file) as f:
		records = f.read()

	if re.search('>', records) == None:
		print('The input RNA sequence must be fasta format.')
		sys.exit(1)
	records = records.split('>')[1:]
	myFasta = []
	for fasta in records:
		array = fasta.split('\n')
		name, sequence = array[0].split()[0], re.sub('[^ACGU-]', '-', ''.join(array[1:]).upper())
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
    ACGU_list=["".join(e) for e in itertools.product(alphabet, repeat=k)]
    return ACGU_list
        
def convert_dict(property_index):

    len_index_value = len(property_index[0])
    k = 0
    for i in range(1, 10):
        if len_index_value < 4**i:
            error_infor = 'error, the number of each index value is must be 4^k.'
            sys.stdout.write(error_infor)
            sys.exit(0)
        if len_index_value == 4**i:
            k = i
            break
    kmer_list = generate_list(k, ALPHABET)   
    len_kmer = len(kmer_list)
    phyche_index_dict = {}
    for kmer in kmer_list:
        phyche_index_dict[kmer] = []
    property_index = list(zip(*property_index))
    for i in range(len_kmer):
        phyche_index_dict[kmer_list[i]] = list(property_index[i])
    return phyche_index_dict

def standard(value_list):

    n = len(value_list)
    average_value = sum(value_list) * 1.0 / n
    std_value=sqrt(sum([pow(e - average_value, 2) for e in value_list]) * 1.0 / (n - 1))
    return std_value

def normalize_index(phyche_index, is_convert_dict=False):

    normalize_phyche = []
    for phyche_value in phyche_index:
        average_phyche = sum(phyche_value) * 1.0 / len(phyche_value)
        sd_phyche = standard(phyche_value)
        normalize_phyche.append([round((e - average_phyche) / sd_phyche, 2) for e in phyche_value])
    if is_convert_dict is True:
        return convert_dict(normalize_phyche)
    return normalize_phyche

def parallel_cor_function(nucleotide1, nucleotide2, phyche_index):
    temp_sum = 0.0
    phyche_index_values = list(phyche_index.values())
    len_phyche_index = len(phyche_index_values[0])
    for u in range(len_phyche_index):
        temp_sum += pow(float(phyche_index[nucleotide1][u]) - float(phyche_index[nucleotide2][u]), 2)

    parallel_value=temp_sum / len_phyche_index
    return parallel_value

def series_cor_function(nucleotide1, nucleotide2, big_lamada, phyche_value):
    
    series_value=float(phyche_value[nucleotide1][big_lamada]) * float(phyche_value[nucleotide2][big_lamada])    
    return series_value

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

def make_pseknc_vector(sequence_list, lamada, w, k, phyche_value, theta_type=1):

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
        if 1 == theta_type:
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

def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer

def RC(kmer):
    myDict = {
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'U': 'A'
    }
    return ''.join([myDict[nc] for nc in kmer[::-1]])


def generateRCKmer(kmerList):
    rckmerList = set()
    myDict = {
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'U': 'A'
    }
    for kmer in kmerList:
        rckmerList.add(sorted([kmer, ''.join([myDict[nc] for nc in kmer[::-1]])])[0])
    return sorted(rckmerList)

#########################################################################
def Kmer(input_data, k=2, normalize=True):
    
    fastas=readRNAFasta(input_data)
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

#vector=Kmer('S.txt',k=4)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('Kmer_out.csv',header=False,index=False)
#
def RCKmer(input_data, k=2, normalize=True):
    
    fastas=readRNAFasta(input_data)
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
#vector=RCKmer('S.txt',k=4)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('RCKmer_out.csv',header=False,index=False)
       
def get_property_dict(): 
    file_path =os.path.split(os.path.realpath(__file__))[0] + r'\data\phy.csv' if platform.system() == 'Windows' else  os.path.split(os.path.realpath(__file__))[0] + '/data/phy.csv'
    phy=pd.read_csv(file_path,header=-1,index_col=None)
    phyche=np.array(phy)
    property_dict=normalize_index(phyche, is_convert_dict=True)
    return property_dict

def Psednc(input_data,lamada=10, w=0.05, k = 2):

    phyche_value = get_property_dict()    
    fastas=readRNAFasta(input_data)
    vector = make_pseknc_vector(fastas, lamada, w, k, phyche_value, theta_type=1)
    return vector
#vector=Psednc('S.txt',lamada=10)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('Psednc_out.csv',header=False,index=False)

def DAC(input_data,k=2,lag=5):

    phyche_value = get_property_dict()     
    fastas=readRNAFasta(input_data)    
    vector,_=make_ac_vector(fastas, lag, phyche_value, k)    
    return vector
#vector=DAC('S.txt',lag=5)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('DAC_out.csv',header=False,index=False)

def DCC(input_data,k=2,lag=10):

    phyche_value = get_property_dict()     
    fastas=readRNAFasta(input_data)  
    vector,_=make_cc_vector(fastas, lag, phyche_value, k)  
    return vector

#vector=DCC('S.txt',lag=10)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('DCC_out.csv',header=False,index=False)

def DACC(input_data, k=2,lag=10):
    
    fastas=readRNAFasta(input_data) 
    phyche_value = get_property_dict() 
    vector1,ac_vector=make_ac_vector(fastas, lag, phyche_value, k)
    vector2,cc_vector=make_cc_vector(fastas, lag, phyche_value, k)    
    zipped = list(zip(vector1,cc_vector))    
    vector = [reduce(lambda x, y: x + y, e) for e in zipped]
    return vector
#vector=DACC('S.txt',lag=10)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('DACC_out.csv',header=False,index=False)

def NAC(input_data):
    fastas=readRNAFasta(input_data)
    encodings = []
    header = ['#']
    for i in ALPHABET:
        header.append(i)
    encodings.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = [name]
        for na in ALPHABET:
            code.append(count[na])
        encodings.append(code)
    return encodings
#vector=NAC('S.txt')
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('NAC_out.csv',header=False,index=False)

def DNC(input_data):
    fastas=readRNAFasta(input_data)
    encodings = []
    dinucleotides = [n1 + n2 for n1 in ALPHABET for n2 in ALPHABET]
    header = ['#'] + dinucleotides
    encodings.append(header)
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
        encodings.append(code)
    return encodings

#vector=DNC('S.txt')
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('DNC_out.csv',header=False,index=False)

def TNC(input_data):
    fastas=readRNAFasta(input_data)
    encodings = []
    triPeptides = [aa1 + aa2 + aa3 for aa1 in ALPHABET for aa2 in ALPHABET for aa3 in ALPHABET]
    header = ['#'] + triPeptides
    encodings.append(header)

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
        encodings.append(code)
    return encodings
#vector=TNC('S.txt')
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
    fastas=readRNAFasta(input_data)
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

#vector=zCurve_vector('S.txt')
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
    fastas=readRNAFasta(input_data)
    vector=[] 
    header=['#']
    for f in range((g-1)*16):
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

#vector=MonoKGap_vector('S.txt',g=2)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('MonoKGap_out.csv',header=False,index=False)   

def MonoDiKGap(x, g):  # 1___2    
    t=[]
    m = list(itertools.product(ALPHABET, repeat=3))
    L_sequence=(len(x)-g-2)
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
    fastas=readRNAFasta(input_data)
    vector=[] 
    header=['#']
    for f in range((g-1)*32):
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

#vector=MonoDiKGap_vector('S.txt',g=2)
#csv_data=pd.DataFrame(data=vector)
#csv_data.to_csv('MonoDiKGap.csv',header=False,index=False)  
    
def binary(input_data):
    fastas=readRNAFasta(input_data)
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
    
    fastas=readRNAFasta(input_data)
    for i in fastas:
        if re.search('[^ACGU-]', i[1]):
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

    nucleotides = ['A', 'C', 'G', 'U']
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
    fastas=readRNAFasta(input_data)
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

    nucleotides = ['A', 'C', 'G', 'U']
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
file = sys.argv[2]

if(option == "1"):
    #kmer method
    vector=Kmer(file, k=4)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('Kmer_out_4.csv',header=False,index=False)
elif(option == "2"):
    #RCKmer method
    vector=RCKmer(file, k=3)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('RCKmer_out.csv',header=False,index=False)
elif(option == "3"):
    #Psednc method
    vector=Psednc(file,lamada=10)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('Psednc_out.csv',header=False,index=False)
elif(option == "4"):
    #DAC method
    vector=DAC(file,lag=5)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DAC_out.csv',header=False,index=False)
elif(option == "5"):
    #DCC method
    vector=DCC(file,k=2,lag=1)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DCC_out.csv',header=False,index=False)
elif(option == "6"):
    #DACC method
    vector=DACC(file,k=2,lag=1)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DACC_out.csv',header=False,index=False)
elif(option == "7"):
    #NAC method
    vector=NAC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('NAC_out.csv',header=False,index=False)
elif(option == "8"):
    #DNC method
    vector=DNC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DNC_out.csv',header=False,index=False)
elif(option == "9"):
    #TNC method
    vector=TNC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('TNC_out.csv',header=False,index=False)
elif(option == "10"):
    #zCurve method
    vector=zCurve_vector(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('zCurve_out.csv',header=False,index=False)
elif(option == "11"):
    #monoMonoKGap method
    vector=MonoKGap_vector(file, g=2)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('MonoKGap_out.csv',header=False,index=False)
elif(option == "12"):
    #monoDiKGap method
    vector=MonoDiKGap_vector(file, g=2)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('MonoDiKGap_out.csv',header=False,index=False)
elif(option == "13"):
    #binary method
    vector=binary(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('binary_out.csv',header=False,index=False)
    
elif(option == "14"):
    #PSTNP method
    vector,A1,A2=PSDNP(file,label)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('PSDNP_out.csv',header=False,index=False)
    
elif(option == "15"):
    #PSDNP method
    vector,A1,A2=PSTNP(file,label)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('PSTNP_out.csv',header=False,index=False)
else:
    print("Invalid method number. Please check the method table!")



