import re, os, sys
from collections import Counter
import numpy as np
import platform
import math
import pandas as pd

ALPHABET='ACDEFGHIKLMNPQRSTVWY'
###############################
##read the protein sequences in fasta format
def readFasta(file):
    with open(file) as f:
         records=f.read()
    if re.search('>',records)==None:
       print('error in fasta format')
       sys.exit(1)
    records=records.split('>')[1:]
    myFasta=[]
    for fasta in records:
        array=fasta.split('\n')
        name, sequence=array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append([name,sequence])
    return myFasta
##################################
#read one sequence in pssm format 
def read_pssm_file(file):
    with open(file) as f:
         records=f.read()
    record=records.split('\n')[3:]
    matrix=[]
    for row in record:
        if len(row)==0:
            break
#        name=row[6]
        value=row[10:]
        value_=value.split()
        #value_list=[float(i) for i in value_]
        matrix.append(value_)
    pssm_matrix=np.array(matrix)
    PSSM_matrix=pssm_matrix.astype(float)
    return PSSM_matrix
################################
#read one sequence in spd format
def read_spd_file(file):
    with open(file) as f:
         records=f.read()
    record=records.split('\n')[1:]
    return record
##################################
#the minimize size of sequence length
def minSequenceLength(fastas):
	min_length = 10000
	for i in fastas:
		if min_length > len(i[1]):
		   min_length = len(i[1])
	return min_length
######################################
#sequence information
def AAC(input_data):
    #This is the amino acid compostion
    #AAC can generate 20-dim vector
    fastas=readFasta(input_data)
    vector=[]
    feature=['#']
    for i in range(20):
        feature.append('AAC.'+str(i))
    vector.append(feature)
    for i in fastas:
        name,sequence=i[0],re.sub('-','',i[1])
        count=Counter(sequence)
        sample=[name]
        for amino in count:
            count[amino]=count[amino]/len(sequence)
        for aa in ALPHABET:
            sample.append(count[aa])
        vector.append(sample)
    return vector

#AAC_out=AAC('protein.txt')
#csv_data=pd.DataFrame(data=AAC_out)
#csv_data.to_csv('AAC_out.csv',header=False,index=False)

def DC(input_data):
    #This is the dipeptide compostion 
    #DC can generate 400-dim vector
    fastas=readFasta(input_data)
    vector=[]
    feature=['#']
    for f in range(400):
        feature.append('DC.'+str(f))
    vector.append(feature)
    Dict={}
    for i in range((len(ALPHABET))):
        Dict[ALPHABET[i]]=i
    for i in fastas:
        name,sequence=i[0],re.sub('-','',i[1])
        sample=[name]
        each_fasta=[0]*400
        for j in range(len(sequence)-2+1):
            each_fasta[Dict[sequence[j]]*20+Dict[sequence[j+1]]]=each_fasta[Dict[sequence[j]]*20+Dict[sequence[j+1]]]+1
        if sum(each_fasta)!=0:
            each_fasta=[i/sum(each_fasta) for i in each_fasta]
        else:
            each_fasta=[0]*400
        sample=sample+each_fasta
        vector.append(sample)
    return vector

#DC_out=DC('protein.txt')
#csv_data=pd.DataFrame(data=DC_out)
#csv_data.to_csv('DC_out.csv',header=False,index=False)

def CKSAAP(input_data, gap=2):
    #This is the CKSAAP discriptor
    #CKSAAP can generate 400*(gap-1)-dim vector
    fastas=readFasta(input_data)
    if gap < 0 or minSequenceLength(fastas) < gap+2:
        print('error: the gap >=0 and sequence length should be > (gap+2)')
        return 0
    vector = []
    protein_pair=[]
    for aa1 in ALPHABET:
        for aa2 in ALPHABET:
            protein_pair.append(aa1 + aa2)
    header = ['#']              
    for f in range(400*gap):
        header.append('CKSAAP.' + str(f))
    vector.append(header)
    for i in fastas:
        name, sequence = i[0], i[1]
        code = [name]
        for g in range(gap+1):
            myDict = {}
            for pair in protein_pair:
                myDict[pair] = 0
        amount = 0
        for index1 in range(len(sequence)):
            index2 = index1 + gap + 1
            if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in ALPHABET and sequence[index2] in ALPHABET:
                myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                amount = amount + 1
        for pair in protein_pair:
            code.append(myDict[pair] / amount)
            vector.append(code)
    return vector

#CKSAAP_out=CKSAAP('protein.txt',gap=2)
#csv_data=pd.DataFrame(data=CKSAAP_out)
#csv_data.to_csv('CKSAAP_out.csv',header=False,index=False)

def GDC(input_data):
    #This is grouped dipeptide compostion
    #GDC can genarate 25-dim vector
    fastas=readFasta(input_data)
    group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
    }
    groupKey = group.keys()
    dipeptide = [g1+'.'+g2 for g1 in groupKey for g2 in groupKey]
    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key
    encodings = []
    header = ['#'] + dipeptide
    encodings.append(header)
    for i in fastas:
	     name, sequence = i[0], re.sub('-', '', i[1])
	     code = [name]
	     myDict = {}
	     for t in dipeptide:
		      myDict[t] = 0

	     sum = 0
	     for j in range(len(sequence) - 2 + 1):
	         myDict[index[sequence[j]]+'.'+index[sequence[j+1]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]] + 1        
	         sum=sum+1
	     if sum==0:
	        for t in dipeptide:
	            code.append(0)       
	     else:
	         for t in dipeptide:
	             code.append(myDict[t]/sum)
	     encodings.append(code)
    return encodings

#GDC_out=GDC('protein.txt')
#csv_data=pd.DataFrame(data=GDC_out)
#csv_data.to_csv('GDC_out.csv',header=False,index=False)

def GTC(input_data):
    #This is the grouped triad composition
    #GTC can genarate 125-dim vector
	fastas=readFasta(input_data)
	group = {
		'alphaticr': 'GAVLMI',
		'aromatic': 'FYW',
		'postivecharger': 'KRH',
		'negativecharger': 'DE',
		'uncharger': 'STCPNQ'
	}
	groupKey = group.keys()
	triple = [g1+'.'+g2+'.'+g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]
	index = {}
	for key in groupKey:
		for aa in group[key]:
			index[aa] = key
	encodings = []
	header = ['#'] + triple
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		myDict = {}
		for t in triple:
			myDict[t] = 0
		sum = 0
		for j in range(len(sequence) - 3 + 1):
			myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] = myDict[index[sequence[j]]+'.'+index[sequence[j+1]]+'.'+index[sequence[j+2]]] + 1
			sum = sum +1
		if sum == 0:
			for t in triple:
				code.append(0)
		else:
			for t in triple:
				code.append(myDict[t]/sum)
		encodings.append(code)
	return encodings

#GTC_out=GTC('protein.txt')
#csv_data=pd.DataFrame(data=GTC_out)
#csv_data.to_csv('GTC_out.csv',header=False,index=False)

def CalculateKSCTriad(sequence, gap, features, AADict):
    #The preprocessing of CTriad
	res = []
	for g in range(gap+1):
		myDict = {}
		for f in features:
			myDict[f] = 0

		for i in range(len(sequence)):
			if i+gap+1 < len(sequence) and i+2*gap+2<len(sequence):
				fea = AADict[sequence[i]] + '.' + AADict[sequence[i+gap+1]]+'.'+AADict[sequence[i+2*gap+2]]
				myDict[fea] = myDict[fea] + 1

		maxValue, minValue = max(myDict.values()), min(myDict.values())
		for f in features:
			res.append((myDict[f] - minValue) / maxValue)
	return res

def KSCTriad(input_data, gap = 0, **kw):
    #This is the k-space conjoint triad 
    #conjiont traid can generate 343-dim vecor
    #KSCTriad can generate 343*gap-dim vector 
	fastas=readFasta(input_data)
	AAGroup = {
		'g1': 'AGV',
		'g2': 'ILFP',
		'g3': 'YMTS',
		'g4': 'HNQW',
		'g5': 'RK',
		'g6': 'DE',
		'g7': 'C'
	}
	myGroups = sorted(AAGroup.keys())
	AADict = {}
	for g in myGroups:
		for aa in AAGroup[g]:
			AADict[aa] = g
	features = [f1 + '.'+ f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]
	encodings = []
	header = ['#']
	for g in range(gap+1):
		for f in features:
			header.append(f+'.gap'+str(g))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		if len(sequence) < 2*gap + 3:
			print('Error: for "KSCTriad" encoding, the input fasta sequences should be greater than (2*gap+3). \n\n')
			return 0
		code = code + CalculateKSCTriad(sequence, gap, features, AADict)
		encodings.append(code)
	return encodings

#CTriad_out=KSCTriad('protein.txt',g=0)
#csv_data=pd.DataFrame(data=CTriad_out)
#csv_data.to_csv('CTriad_out.csv',header=False,index=False)
#KSCTriad_out=KSCTriad('protein.txt',g=1)
#csv_data=pd.DataFrame(data=KSCTriad_out)
#csv_data.to_csv('KSCTriad_out.csv',header=False,index=False)

def get_CTD_group():
    #the goroups of amino acid sequence
    group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	   }
    group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	    }
    group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	   }
    return group1,group2,group3
def Count(seq1, seq2):
    #the counting process of CTD
	amount = 0
	for aa in seq1:
		amount = amount + seq2.count(aa)
	return amount

def CTDC(input_data):
    #This is the composition,transition and distribution 
    #CTDC can can produce 39-dim vector
	fastas=readFasta(input_data)
	group1, group2, group3=get_CTD_group()
	groups = [group1, group2, group3]
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	encodings = []
	header = ['#']
	for p in property:
		for g in range(1, len(groups) + 1):
			header.append(p + '.G' + str(g))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('X', '', i[1])
		code = [name]
		for p in property:
			c1 = Count(group1[p], sequence) / len(sequence)
			c2 = Count(group2[p], sequence) / len(sequence)
			c3 = 1 - c1 - c2
			code = code + [c1, c2, c3]
		encodings.append(code)
	return encodings

#CTDC_out=CTDC('protein.txt')
#csv_data=pd.DataFrame(data=CTDC_out)
#csv_data.to_csv('CTDC_out.csv',header=False,index=False)

def CTDT(input_data):
    #This is the compostion,transition and distribution
    #CTDD can generate 195-dim vector
	fastas=readFasta(input_data) 
	group1, group2, group3=get_CTD_group()
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')
	encodings = []
	header = ['#']
	for p in property:
		for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
			header.append(p + '.' + tr)
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('X', '', i[1])
		code = [name]
		aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
		for p in property:
			c1221, c1331, c2332 = 0, 0, 0
			for pair in aaPair:
				if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
					c1221 = c1221 + 1
					continue
				if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
					c1331 = c1331 + 1
					continue
				if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
					c2332 = c2332 + 1
			code = code + [c1221/len(aaPair), c1331/len(aaPair), c2332/len(aaPair)]
		encodings.append(code)
	return encodings

#CTDT_out=CTDT('protein.txt')
#csv_data=pd.DataFrame(data=CTDT_out)
#csv_data.to_csv('CTDT_out.csv',header=False,index=False)

def Count_D(aaSet, sequence):
    number=0
    code=[]
    select=[]
    for aa in sequence:
        if aa in aaSet:
            number=number+1
    cutoffNums=[1, math.floor(0.25*number),math.floor(0.5*number),math.floor(0.75*number),number]
    myCount=0
    for i in range(len(sequence)):
        if sequence[i] in aaSet:
            myCount=myCount+1
            if myCount in cutoffNums:
                code.append((i+1)/len(sequence))
    if len(code)<5:
        code=[]
        for i in range(len(sequence)):          
            if sequence[i] in aaSet:
               select.append(i)
        if len(select)<1:
            code=[0,0,0,0,0]  
        else:
            if 0 in cutoffNums:
               cutoffNums=np.array(cutoffNums)
               cutoffNums[cutoffNums==0]=1
            for j in range(5):
                label=select[cutoffNums[j]-1]
                code.append((label+1)/len(sequence))
                label=[]
    return code

def CTDD(input_data):
    #This is the composition transition and distribution 
    #CTDD can generate 39-dim vector
	fastas=readFasta(input_data)     
	group1, group2, group3=get_CTD_group()
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')
	encodings = []
	header = ['#']
	for p in property:
		for g in ('1', '2', '3'):
			for d in ['0', '25', '50', '75', '100']:
				header.append(p + '.' + g + '.residue' + d)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('X', '', i[1])
		code = [name]
		for p in property:
			code = code + Count_D(group1[p], sequence) + Count_D(group2[p], sequence) + Count_D(group3[p], sequence)
		encodings.append(code)
	return encodings

#CTDD_out=CTDD('protein.txt')
#csv_data=pd.DataFrame(data=CTDD_out)
#csv_data.to_csv('CTDD_out.csv',header=False,index=False)

def CalculateEBGW(sequence): 
    groups={'g1': 'GAVLIMPFW',
            'g2': 'QNSTYC',
            'g3': 'DE',
            'g4': 'HKR'}
    fenlei1={'g1','g2'}
    fenlei2={'g1','g3'}
    fenlei3={'g1','g4'}
    sum1=0
    sum_one,sum_two,sum_three=[],[],[]
    for g in fenlei1:
        sum1 = sum1+Count(groups[g],sequence)
    sum_one=(sum1)/len(sequence)
    sum2=0
    for g in fenlei2:
        sum2 = sum2+Count(groups[g],sequence)
    sum_two=(sum2)/len(sequence)   
    sum3=0
    for g in fenlei3:
        sum3 = sum3+Count(groups[g],sequence)
    sum_three=(sum3)/len(sequence) 
    code=[sum_one,sum_two,sum_three]
    return code

def EBGW(input_data,L):
    #This is the Encoding Based on Grouped Weight descriptor
    #EBGW can generate L*3 vector
    fastas=readFasta(input_data)    
    vector=[]  
    header=['#']            
    for f in range(L*3):
        header.append('EBGW'+str(f))
    vector.append(header)     
    for i in fastas:
        name,sequence=i[0],re.sub('X','',i[1])
        number=len(sequence)
        cutnumber=number/L
        cutoff=np.array(range(1,L+1))*cutnumber
        cutoffnum=[math.floor(cutoff[m]) for m in range(len(cutoff))]
        vector1=[]
        vector2=[name]
        for j in range(len(cutoffnum)):
            seq_new=sequence[:cutoffnum[j]]
            vector1=CalculateEBGW(seq_new)
            vector2=vector2+vector1
        vector.append(vector2)
    return vector

#EBGW_out=EBGW('protein.txt',L=5)
#csv_data=pd.DataFrame(data=EBGW_out)
#csv_data.to_csv('EBGW_out.csv',header=False,index=False)

####################################################################
#Physicochemical information
def AC(input_data, props=['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
						 'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201'],
				nlag = 30):
    # This is the Auto Covariance (AC) descriptor
    #AC can produce 8*nlag vector
	fastas=readFasta(input_data)
	if minSequenceLength(fastas) < nlag + 1:		
		print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
		return 0
	fileAAidx = os.path.split(os.path.realpath(__file__))[0] + r'\data\AAidx.txt' if platform.system() == 'Windows' else sys.path[0] + '/data/AAidx.txt'

	with open(fileAAidx) as f:
		records = f.readlines()[1:]
	myDict = {}
	for i in records:
		array = i.rstrip().split('\t')
		myDict[array[0]] = array[1:]
	AAidx = []
	AAidxName = []
	for i in props:
		if i in myDict:
			AAidx.append(myDict[i])
			AAidxName.append(i)
		else:
			print('"' + i + '" properties not exist.')
			return None
	AAidx1 = np.array([float(j) for i in AAidx for j in i])
	AAidx = AAidx1.reshape((len(AAidx), 20))
	propMean = np.mean(AAidx,axis=1)
	propStd = np.std(AAidx, axis=1)
	for i in range(len(AAidx)):
		for j in range(len(AAidx[i])):
			AAidx[i][j] = (AAidx[i][j] - propMean[i]) / propStd[i]
	index = {}
	for i in range(len(ALPHABET)):
		index[ALPHABET[i]] = i
	encodings = []
	header = ['#']
	for p in props:
		for n in range(1, nlag+1):
			header.append(p + '.lag' + str(n))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		N = len(sequence)
		for prop in range(len(props)):
			xmean = sum([AAidx[prop][index[aa]] for aa in sequence]) / N
			for n in range(1, nlag + 1):
				if len(sequence) > nlag:
					# if key is '-', then the value is 0
					jieguo = sum([(AAidx[prop][index.get(sequence[j], 0)] - xmean) * (AAidx[prop][index.get(sequence[j + n], 0)] - xmean) for j in range(len(sequence) - n)]) / (N - n)
					rn = jieguo
				else:
					rn = 'NA'
				code.append(rn)
		encodings.append(code)
	return encodings

#AC_out=AC('protein.txt')
#csv_data=pd.DataFrame(data=AC_out)
#csv_data.to_csv('AC_out.csv',header=False,index=False)

def NMBroto(input_data, props=['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
										 'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201'],
				nlag = 30):
    
	fastas=readFasta(input_data)
	if minSequenceLength(fastas) < nlag + 1:
		print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
		return 0
	fileAAidx = os.path.split(os.path.realpath(__file__))[0] + r'\data\AAidx.txt' if platform.system() == 'Windows' else sys.path[0] + '/data/AAidx.txt'
	with open(fileAAidx) as f:
		records = f.readlines()[1:]
	myDict = {}
	for i in records:
		array = i.rstrip().split('\t')
		myDict[array[0]] = array[1:]
	AAidx = []
	AAidxName = []
	for i in props:
		if i in myDict:
			AAidx.append(myDict[i])
			AAidxName.append(i)
		else:
			print('"' + i + '" properties not exist.')
			return None
	AAidx1 = np.array([float(j) for i in AAidx for j in i])
	AAidx = AAidx1.reshape((len(AAidx),20))
	pstd = np.std(AAidx, axis=1)
	pmean = np.average(AAidx, axis=1)
	for i in range(len(AAidx)):
		for j in range(len(AAidx[i])):
			AAidx[i][j] = (AAidx[i][j] - pmean[i]) / pstd[i]
	index = {}
	for i in range(len(ALPHABET)):
		index[ALPHABET[i]] = i
	encodings = []
	header = ['#']
	for p in props:
		for n in range(1, nlag + 1):
			header.append(p + '.lag' + str(n))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		N = len(sequence)
		for prop in range(len(props)):
			for n in range(1, nlag + 1):
				if len(sequence) > nlag:
					# if key is '-', then the value is 0
					rn = sum([AAidx[prop][index.get(sequence[j], 0)] * AAidx[prop][index.get(sequence[j + n], 0)] for j in range(len(sequence)-n)]) / (N - n)
				else:
					rn = 'NA'
				code.append(rn)
		encodings.append(code)
	return encodings

#NMBroto_out=NMBroto('protein.txt')
#csv_data=pd.DataFrame(data=NMBroto_out)
#csv_data.to_csv('NMBroto_out.csv',header=False,index=False)

def Moran(input_data, props=['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
						 'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201'],
				nlag = 30):
	fastas=readFasta(input_data)
	if minSequenceLength(fastas) < nlag + 1:		
		print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
		return 0
	fileAAidx = os.path.split(os.path.realpath(__file__))[0] + r'\data\AAidx.txt' if platform.system() == 'Windows' else sys.path[0] + '/data/AAidx.txt'

	with open(fileAAidx) as f:
		records = f.readlines()[1:]
	myDict = {}
	for i in records:
		array = i.rstrip().split('\t')
		myDict[array[0]] = array[1:]
	AAidx = []
	AAidxName = []
	for i in props:
		if i in myDict:
			AAidx.append(myDict[i])
			AAidxName.append(i)
		else:
			print('"' + i + '" properties not exist.')
			return None
	AAidx1 = np.array([float(j) for i in AAidx for j in i])
	AAidx = AAidx1.reshape((len(AAidx), 20))
	propMean = np.mean(AAidx,axis=1)
	propStd = np.std(AAidx, axis=1)
	for i in range(len(AAidx)):
		for j in range(len(AAidx[i])):
			AAidx[i][j] = (AAidx[i][j] - propMean[i]) / propStd[i]
	index = {}
	for i in range(len(ALPHABET)):
		index[ALPHABET[i]] = i
	encodings = []
	header = ['#']
	for p in props:
		for n in range(1, nlag+1):
			header.append(p + '.lag' + str(n))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		N = len(sequence)
		for prop in range(len(props)):
			xmean = sum([AAidx[prop][index[aa]] for aa in sequence]) / N
			for n in range(1, nlag + 1):
				if len(sequence) > nlag:
					# if key is '-', then the value is 0
					fenzi = sum([(AAidx[prop][index.get(sequence[j], 0)] - xmean) * (AAidx[prop][index.get(sequence[j + n], 0)] - xmean) for j in range(len(sequence) - n)]) / (N - n)
					fenmu = sum([(AAidx[prop][index.get(sequence[j], 0)] - xmean) ** 2 for j in range(len(sequence))]) / N
					rn = fenzi / fenmu
				else:
					rn = 'NA'
				code.append(rn)
		encodings.append(code)
	return encodings

#Moran_out=Moran('protein.txt')
#csv_data=pd.DataFrame(data=Moran_out)
#csv_data.to_csv('Moran_out.csv',header=False,index=False)

def Geary(input_data, props=['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
						 'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201'],
				nlag = 30):
    fastas=readFasta(input_data)	
    if minSequenceLength(fastas) < nlag + 1:
        print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
        return 0
    AAidx_file=os.path.split(os.path.realpath(__file__))[0] + r'\data\AAidx.txt' if platform.system() == 'Windows' else sys.path[0] + '/data/AAidx.txt' 
    with open(AAidx_file) as f: 
         records=f.readlines()[1:]
    myDict = {}
    for i in records:
        array=i.rstrip().split('\t')
        myDict[array[0]]=array[1:]
    AAidx = []
    AAidxName = []
    for i in props:
        for i in myDict:
            AAidx.append(myDict[i])
            AAidxName.append(i)
        else:
            print(i + 'properties not exist')
            return None
    AAdix1=np.array([float(j) for i in AAidx for j in i])
    AAidx=AAdix1.reshape(len(AAidx),20)
    propMean = np.mean(AAidx, axis=1)
    propStd=np.std(AAidx,axis=1)
    for i in range(len(AAidx)):
        for j in range(len(AAidx[i])):
            AAidx[i][j]=(AAidx[i][j]-propMean[i])/propStd[i]
    index={}
    for i in range(len(ALPHABET)):
        index[ALPHABET[i]]=i
    encodings=[]
    header=['#']
    for p in props:
        for n in range(1,nlag+1):
            header.append(p+str(n))
    encodings.append(header)
    for i in fastas:
        name,sequence=i[0],re.sub('-','',i[1])
        code=[name]
        N=len(sequence)
        for prop in range(len(props)):
            xmean=sum([AAidx[prop][index[aa]] for aa in sequence])/N
            for n in range(1,nlag+1):
                if len(sequence)>nlag:
                    rn=(N-1)/(2*(N-n)) * ((sum([(AAidx[prop][index.get(sequence[j], 0)] - AAidx[prop][index.get(sequence[j + n], 0)])**2 for j in range(len(sequence)-n)])) / (sum([(AAidx[prop][index.get(sequence[j], 0)] - xmean) ** 2 for j in range(len(sequence))])))
                else:
                    rn='NAN'
                code.append(rn)
        encodings.append(code)
    return encodings

#Geary_out=Geary('protein.txt')
#csv_data=pd.DataFrame(data=Geary_out)
#csv_data.to_csv('Geary_out.csv',header=False,index=False)

def QSOrder(input_data, nlag=30, w=0.1):
	fastas=readFasta(input_data)	
	if minSequenceLength(fastas) < nlag + 1:
		print('Error: all the sequence length should be larger than the nlag+1: ' + str(nlag + 1) + '\n\n')
		return 0
	dataFile = os.path.split(os.path.realpath(__file__))[0] + r'\data\Schneider-Wrede.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/Schneider-Wrede.txt'
	dataFile1 = os.path.split(os.path.realpath(__file__))[0] + r'\data\Grantham.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/Grantham.txt'
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	AA1 = 'ARNDCQEGHILKMFPSTWYV'
	DictAA = {}
	for i in range(len(AA)):
		DictAA[AA[i]] = i
	DictAA1 = {}
	for i in range(len(AA1)):
		DictAA1[AA1[i]] = i
	with open(dataFile) as f:
		records = f.readlines()[1:]
	AADistance = []
	for i in records:
		array = i.rstrip().split()[1:] if i.rstrip() != '' else None
		AADistance.append(array)
	AADistance = np.array(
		[float(AADistance[i][j]) for i in range(len(AADistance)) for j in range(len(AADistance[i]))]).reshape((20, 20))
	with open(dataFile1) as f:
		records = f.readlines()[1:]
	AADistance1 = []
	for i in records:
		array = i.rstrip().split()[1:] if i.rstrip() != '' else None
		AADistance1.append(array)
	AADistance1 = np.array(
		[float(AADistance1[i][j]) for i in range(len(AADistance1)) for j in range(len(AADistance1[i]))]).reshape(
		(20, 20))
	encodings = []
	header = ['#']
	for aa in AA1:
		header.append('Schneider.Xr.' + aa)
	for aa in AA1:
		header.append('Grantham.Xr.' + aa)
	for n in range(1, nlag + 1):
		header.append('Schneider.Xd.' + str(n))
	for n in range(1, nlag + 1):
		header.append('Grantham.Xd.' + str(n))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		arraySW = []
		arrayGM = []
		for n in range(1, nlag + 1):
			arraySW.append(
				sum([AADistance[DictAA[sequence[j]]][DictAA[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]))
			arrayGM.append(sum(
				[AADistance1[DictAA1[sequence[j]]][DictAA1[sequence[j + n]]] ** 2 for j in range(len(sequence) - n)]))
		myDict = {}
		for aa in AA1:
			myDict[aa] = sequence.count(aa)
		for aa in AA1:
			code.append(myDict[aa] / (1 + w * sum(arraySW)))
		for aa in AA1:
			code.append(myDict[aa] / (1 + w * sum(arrayGM)))
		for num in arraySW:
			code.append((w * num) / (1 + w * sum(arraySW)))
		for num in arrayGM:
			code.append((w * num) / (1 + w * sum(arrayGM)))
		encodings.append(code)
	return encodings

#QSOrder_out=QSOrder('protein.txt')
#csv_data=pd.DataFrame(data=QSOrder_out)
#csv_data.to_csv('QSOrder_out.csv',header=False,index=False)

def Rvalue(aa1, aa2, AADict, Matrix):
    R_value=sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)
    return R_value

def PAAC(input_data, lambdaValue=30, w=0.05):
	fastas=readFasta(input_data)
	if minSequenceLength(fastas) < lambdaValue + 1:
		print('Error: all the sequence length should be larger than the lambdaValue+1: ' + str(lambdaValue + 1) + '\n\n')
		return 0
	dataFile = os.path.split(os.path.realpath(__file__))[0] + r'\data\PAAC.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/PAAC.txt'
	with open(dataFile) as f:
		records = f.readlines()
	AA = ''.join(records[0].rstrip().split()[1:])
	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i
	AAProperty = []
	AAPropertyNames = []
	for i in range(1, len(records)):
		array = records[i].rstrip().split() if records[i].rstrip() != '' else None
		AAProperty.append([float(j) for j in array[1:]])
		AAPropertyNames.append(array[0])
	AAProperty1 = []
	for i in AAProperty:
		meanI = sum(i) / 20
		fenmu = math.sqrt(sum([(j-meanI)**2 for j in i])/20)
		AAProperty1.append([(j-meanI)/fenmu for j in i])
	encodings = []
	header = ['#']
	for aa in AA:
		header.append('PAAC.' + aa)
	for n in range(1, lambdaValue + 1):
		header.append('PAAC.lambda' + str(n))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		theta = []
		for n in range(1, lambdaValue + 1):
			theta.append(
				sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (
				len(sequence) - n))
		myDict = {}
		for aa in AA:
			myDict[aa] = sequence.count(aa)/100
		code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
		code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
		encodings.append(code)
	return encodings

#PAAC_out=PAAC('protein.txt')
#csv_data=pd.DataFrame(data=PAAC_out)
#csv_data.to_csv('PAAC_out.csv',header=False,index=False)

def APAAC(input_data, lambdaValue=30, w=0.05):
	fastas=readFasta(input_data)   
	if minSequenceLength(fastas) < lambdaValue + 1:
		print('Error: all the sequence length should be larger than the lambdaValue+1: ' + str(lambdaValue + 1) + '\n\n')
		return 0
	dataFile = os.path.split(os.path.realpath(__file__))[0] + r'\data\PAAC.txt' if platform.system() == 'Windows' else re.sub('codes$', '', os.path.split(os.path.realpath(__file__))[0]) + '/data/PAAC.txt'
	with open(dataFile) as f:
		records = f.readlines()
	AA = ''.join(records[0].rstrip().split()[1:])
	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i
	AAProperty = []
	AAPropertyNames = []
	for i in range(1, len(records) - 1):
		array = records[i].rstrip().split() if records[i].rstrip() != '' else None
		AAProperty.append([float(j) for j in array[1:]])
		AAPropertyNames.append(array[0])
	AAProperty1 = []
	for i in AAProperty:
		meanI = sum(i) / 20
		fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
		AAProperty1.append([(j - meanI) / fenmu for j in i])
	encodings = []
	header = ['#']
	for i in AA:
		header.append('Pc1.' + i)
	for j in range(1, lambdaValue + 1):
		for i in AAPropertyNames:
			header.append('Pc2.' + i + '.' + str(j))
	encodings.append(header)
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		theta = []
		for n in range(1, lambdaValue + 1):
			for j in range(len(AAProperty1)):
				theta.append(sum([AAProperty1[j][AADict[sequence[k]]] * AAProperty1[j][AADict[sequence[k + n]]] for k in
								  range(len(sequence) - n)]) / (len(sequence) - n))
		myDict = {}
		for aa in AA:
			myDict[aa] = sequence.count(aa)

		code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
		code = code + [w * value / (1 + w * sum(theta)) for value in theta]
		encodings.append(code)
	return encodings

#APAAC_out=APAAC('protein.txt')
#csv_data=pd.DataFrame(data=APAAC_out)
#csv_data.to_csv('APAAC_out.csv',header=False,index=False)

############################################################
#evolutionary information
def average(matrix, seqLen):
    # average the summary of rows
    matrix_array = np.array(matrix)
    matrix_array = np.divide(matrix_array, seqLen)
    matrix_array_shp = np.shape(matrix_array)
    matrix_average = [(np.reshape(matrix_array, (matrix_array_shp[0] * matrix_array_shp[1], )))]
    return matrix_average

def normalizePSSM(PSSM):    
    PSSM=np.array(PSSM)
    seq_cn=np.shape(PSSM)[0]
    PSSM_norm=[ [0.0] * 20 ] * seq_cn
    PSSM_norm=np.array(PSSM_norm)
    mean_matrix=np.mean(PSSM, axis=1)
    std_matrix=np.std(PSSM, axis=1)
    for i in range(seq_cn):
        for j in range(20):
            if std_matrix[i]==0.0:
                PSSM_norm[i][j] = PSSM[i][j]-mean_matrix[i]
            else:
                PSSM_norm[i][j]=(PSSM[i][j]-mean_matrix[i])/std_matrix[i]
    return PSSM_norm

def exponentPSSM(PSSM):
    PSSM=np.array(PSSM)
    seq_cn=np.shape(PSSM)[0]
    PSSM_exponent=[ [0.0] * 20 ] * seq_cn
    for i in range(seq_cn):
        for j in range(20):
            PSSM_exponent[i][j]=math.exp(PSSM[i][j])
    PSSM_exponent=np.array(PSSM_exponent)
    return  PSSM_exponent

def preHandleColumns(PSSM,STEP,ID):
    PSSM=PSSM.astype(float)
    matrix_final = [ [0.0] * 20 ] * 20
    matrix_final=np.array(matrix_final)
    seq_cn=np.shape(PSSM)[0]
    if ID==0:
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    matrix_final[i][j]+=(PSSM[k][i]*PSSM[k+STEP][j])
    elif ID==1:
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    matrix_final[i][j] += ((PSSM[k][i]-PSSM[k+STEP][j]) * (PSSM[k][i]-PSSM[k+STEP][j])/4.0)
    return matrix_final

def handleMixed(PSSM,ALPHA):    
    row1=[0.0] * 20
    row2=[0.0] * 20
    matrix_final = [ [0.0] * 40 ] * 1
    row1=np.array(row1)
    row2=np.array(row2)
    matrix_final=np.array(matrix_final)
    PSSM_norm=normalizePSSM(PSSM)#normalize
    seq_cn=np.shape(PSSM)[0]
    for i in range(seq_cn):
        row1=map(sum, zip(row1, PSSM_norm[i]))
    row1=np.array(list(row1))
    row1=np.divide(row1,seq_cn)
    for j in range(20):
        for i in range(seq_cn-ALPHA):
            row2[j]+=(PSSM_norm[i][j]-PSSM_norm[i+ALPHA][j])*(PSSM_norm[i][j]-PSSM_norm[i+ALPHA][j])
    row2=np.array(list(row2))
    row2=np.divide(row2,seq_cn-ALPHA)
    row=np.hstack((row1,row2))
    matrix_final[0]=row
    return matrix_final

def correlation(PSSM,ID,GROUP):
    PSSM=np.array(PSSM)
    seq_cn = np.shape(PSSM)[0]
    g=GROUP
    l=seq_cn
    if ID==0:
        matrix_final=pssm_ac_cal(PSSM,g,l)
    elif ID==1:
        matrix_final=pssm_cc_cal(PSSM,g,l)
    matrix_final =  average(matrix_final, l)
    return matrix_final

def pssm_ac_cal(PSSM, g, l):
    PSSM_AC = np.array([ [0.0] * 20 ] * g)
    for pg in range(g):
        l_g = l - pg - 1
        for pj in range(20):
            sum_jl = 0.0
            for i in range(l):
                sum_jl += PSSM[i][pj]
            sum_jl /= l
            pssm_acjg = 0.0
            for i in range(l_g):
                pssm_acjg +=  (PSSM[i][pj]-sum_jl) * (PSSM[i+pg+1][pj]-sum_jl)
            pssm_acjg /= l_g
            PSSM_AC[pg][pj] = pssm_acjg
    return PSSM_AC

def pssm_cc_cal(PSSM, g, l):
    PSSM_CC = np.array([ [0.0] * 380 ] * g)
    for pg in range(g):
        l_g = l - pg - 1
        for pj_1 in range(20):
            sum_jl_1 = 0.0
            for i in range(l):
                sum_jl_1 += PSSM[i][pj_1]
            sum_jl_1 /= l
            for pj_2 in range(20):
                if pj_2!=pj_1:
                    sum_jl_2 = 0.0
                    for i in range(l):
                        sum_jl_2 += PSSM[i][pj_2]
                    sum_jl_2 /= l
                    pssm_acjg = 0.0
                    for i in range(l_g):
                        pssm_acjg += (PSSM[i][pj_1]-sum_jl_1) * (PSSM[i+pg+1][pj_2]-sum_jl_2)
                    pssm_acjg /= l_g
                    if(pj_1<pj_2):
                        PSSM_CC[pg][19*pj_1+(pj_2-1)] = pssm_acjg
                    else:
                        PSSM_CC[pg][19*pj_1+(pj_2)] = pssm_acjg
    return PSSM_CC

def aac_pssm(input_matrix,exp=True):
    if exp==True:
       input_matrix=exponentPSSM(input_matrix)
    else:
        input_matrix=input_matrix
    seq_cn=float(np.shape(input_matrix)[0])
    aac_pssm_matrix=input_matrix.sum(axis=0)
    aac_pssm_vector=aac_pssm_matrix/seq_cn
    vec=[]
    result=[]
    header=[]
    for f in range(20):
        header.append('aac_pssm.'+str(f))
    result.append(header)   
    for v in aac_pssm_vector:
        vec.append(v)
    result.append(vec)
    return aac_pssm_vector,result

def dpc_pssm(input_matrix,exp=True):
    if exp==True:
       input_matrix=exponentPSSM(input_matrix)
    else:
        input_matrix=input_matrix
    STEP = 1
    ID = 0
    matrix_final = preHandleColumns(input_matrix, STEP, ID)
    seq_cn = float(np.shape(input_matrix)[0])
    dpc_pssm_vector = average(matrix_final, seq_cn-STEP)
    vector1=np.array(dpc_pssm_vector).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    vec=[]
    result=[]
    header=[]
    for f in range(400):
        header.append('dpc_pssm.'+str(f))
    result.append(header)
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return dpc_pssm_vector,result

def pse_pssm(input_matrix, ALPHA,exp=True):
    if exp==True:
       input_matrix=exponentPSSM(input_matrix)
    else:
        input_matrix=input_matrix
    pse_pssm_matrix=handleMixed(input_matrix,ALPHA)
    vector1=np.array(pse_pssm_matrix).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    vec=[]
    result=[]
    header=[]
    for f in range(40):
        header.append('pse_pssm.'+str(f))
    result.append(header)
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return pse_pssm_matrix,result

def ab_pssm(input_matrix,exp=True):
    if exp==True:
       input_matrix=exponentPSSM(input_matrix)
    else:
        input_matrix=input_matrix
    seq_cn=np.shape(input_matrix)[0]
    BLOCK=np.floor(seq_cn/20)
    BLOCK=int(BLOCK)
    matrix_final=[]
    for i in range(19):
        tmp=input_matrix[i*BLOCK:(i+1)*BLOCK]
        aac_matrix,_=aac_pssm(tmp,exp=False)
        matrix_final.append(aac_matrix)
    tmp=input_matrix[19*BLOCK:]
    aac_matrix_,_=aac_pssm(tmp,exp=False)
    matrix_final.append(aac_matrix_)
    ab_pssm_matrix=average(matrix_final,1.0)
    vector1=np.array(ab_pssm_matrix).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    vec=[]
    result=[]
    header=[]
    for f in range(400):
        header.append('ab_pssm.'+str(f))
    result.append(header)
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return ab_pssm_matrix,result

def ac_pssm(input_matrix, GROUP,exp=True):
    if exp==True:
       input_matrix=exponentPSSM(input_matrix)
    else:
        input_matrix=input_matrix
    ID=0
    pssm_ac_matrix=correlation(input_matrix, ID, GROUP)
    vector1=np.array(pssm_ac_matrix).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    vec=[]
    result=[]
    header=[]
    for f in range(20*GROUP):
        header.append('ac_pssm.'+str(f))
    result.append(header)
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return pssm_ac_matrix,result

def cc_pssm(input_matrix, GROUP,exp=True):
    if exp==True:
       input_matrix=exponentPSSM(input_matrix)
    else:
        input_matrix=input_matrix
    ID=1
    pssm_cc_matrix=correlation(input_matrix, ID, GROUP)
    vector1=np.array(pssm_cc_matrix).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    vec=[]
    result=[]
    header=[]
    for f in range(380*GROUP):
        header.append('cc_pssm.'+str(f))
    result.append(header)
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return pssm_cc_matrix,result

def acc_pssm(input_matrix, GROUP,exp=True):
    ID1=0
    ID2=1
    if exp==True:
       input_matrix=exponentPSSM(input_matrix)
    else:
        input_matrix=input_matrix
    pssm_ac_matrix=correlation(input_matrix, ID1, GROUP)
    pssm_cc_matrix=correlation(input_matrix, ID2, GROUP)
    pssm_acc_matrix=np.concatenate((pssm_ac_matrix,pssm_cc_matrix),axis=1)
    vector1=np.array(pssm_acc_matrix).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    vec=[]
    result=[]
    header=[]
    for f in range(400*GROUP):
        header.append('acc_pssm.'+str(f))
    result.append(header)
    for v in vector3:
        vec.append(v)
    result.append(vec)    
    return pssm_acc_matrix,result

def bi_pssm(input_matrix,exp=True):
    if exp==True:
       input_matrix=exponentPSSM(input_matrix)
    else:
        input_matrix=input_matrix
    PSSM=input_matrix
    PSSM=np.array(PSSM)
    header=[]
    for f in range(400):
        header.append('bi_pssm.'+str(f))
    result=[]
    result.append(header)
    vec=[]
    bipssm=[[0.0]*400]*(PSSM.shape[0]-1)
    p=0
    for i in range(20):
        for j in range(20):
            for h in range(PSSM.shape[0]-1):
                bipssm[h][p]=PSSM[h][i]*PSSM[h+1][j]
            p=p+1
    vector=np.sum(bipssm,axis=0)
    for v in vector:
        vec.append(v)
    result.append(vec)   
    return vector,result

############################################################
#Structural information
def SSC(file):
    spd_data=read_spd_file(file)
    C_count,E_count,H_count= 0,0,0
    header=[]
    for f in range(3):
        header.append('SSC.'+str(f))
    result=[]
    result.append(header)
    vec=[]
    L_sequence=0
    for z in range(len(spd_data)):
        seq=spd_data[z]
        if len(seq)==0:
           break
        L_sequence+=1
        matrix = seq.split()
        if matrix[2] == 'C':
           C_count += 1
        if matrix[2] == 'H':
           H_count += 1
        if matrix[2] == 'E':
           E_count += 1   
    C_count = C_count/L_sequence
    E_count = E_count/L_sequence
    H_count = H_count/L_sequence 
    vector=[C_count,E_count,H_count]
    for v in vector:
        vec.append(v)
    result.append(vec)
    return vector,result

#vector,result=[],[]
#vector,result=SSC(file)
    
def ASA(file):
    spd_data=read_spd_file(file)
    header=[]
    for f in range(1):
        header.append('ASA.'+str(f))
    result=[]
    result.append(header)
    ASA_sum =0
    L_sequence=0
    for z in range(len(spd_data)):
        seq=spd_data[z]
        L_sequence+=1
        if len(seq)==0:
           break
        matrix = seq.split()
        ASA_sum += float(matrix[3])
    vector = [ASA_sum/L_sequence]
    result.append(vector)
    return vector,result

#vector,result=[],[]
#vector,result=ASA(file)

def TAC(file):
    spd_data=read_spd_file(file)
    header=[]
    for f in range(8):
        header.append('TAC.'+str(f))
    result=[]
    result.append(header)
    vec=[]
    psi_sin_sum,psi_cos_sum,phi_sin_sum,phi_cos_sum=0,0,0,0
    theta_sin_sum,theta_cos_sum=0,0
    tau_sin_sum,tau_cos_sum=0,0
    L_sequence=0
    for z in range(len(spd_data)):
        seq=spd_data[z]
        if len(seq)==0:
           break
        L_sequence+=1
        matrix = seq.split()
        phi_sin_sum += np.sin(float(matrix[4]))
        phi_cos_sum += np.cos(float(matrix[4]))
        psi_sin_sum += np.sin(float(matrix[5]))
        psi_cos_sum += np.sin(float(matrix[5]))
        theta_sin_sum += np.sin(float(matrix[6]))
        theta_cos_sum += np.sin(float(matrix[6]))
        tau_sin_sum += np.sin(float(matrix[7]))
        tau_cos_sum += np.cos(float(matrix[7]))
    phi_sin = phi_sin_sum/L_sequence
    phi_cos = phi_cos_sum / L_sequence
    psi_sin = psi_sin_sum / L_sequence
    psi_cos = psi_cos_sum / L_sequence
    theta_sin = theta_sin_sum/ L_sequence
    theta_cos = theta_cos_sum /L_sequence
    tau_sin = tau_sin_sum/L_sequence
    tau_cos = tau_cos_sum / L_sequence  
    vector = [phi_sin, phi_cos , psi_sin ,psi_cos ,theta_sin, theta_cos,tau_sin,tau_cos]
    for v in vector:
        vec.append(v)
    result.append(vec)
    return vector,result

#vector,result=[],[]
#vector,result=TAC(file)

def TAAC(file): 
    spd_data=read_spd_file(file)
    header=[]
    for f in range(80):
        header.append('TAAC.'+str(f))
    result=[]
    result.append(header)
    vec=[]
    L_sequence=0
    i = 0
    feature = []
    for z in range(len(spd_data)):
        seq=spd_data[z]
        if len(seq)==0:
           break
        L_sequence+=1
        seq=seq.split()
        spd_a=seq[3:]   
        feature.append(spd_a)
    matrix = [[0 for x in range(8)] for y in range(10)]
    degree_matrix = [[0 for x in range(8)] for y in range(L_sequence-1)]
    degree_index = 0
    for x in range(0,L_sequence-1):
        degree_index = 0
        for y in range(1,5):
            degree_matrix[x][degree_index] = np.math.sin(float(feature[x][y]) * np.pi / 180 )
            degree_matrix[x][degree_index+1] = np.math.cos(float(feature[x][y]) * np.pi / 180 )
            degree_index += 2
    array = degree_matrix
    for k in range(0, 10):
        for j in range(0, 8):
            for i in range(0, L_sequence - 1 - k):
                matrix[k][j] += float(array[i][j]) * float(array[i + k][j])
                matrix[k][j] = matrix[k][j] / (L_sequence - 1)
    matrix=np.array(matrix)
    vector=np.reshape(matrix,(1,matrix.shape[0]*matrix.shape[1]))
    vector1=np.array(vector).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return vector,result

#vector,result=[],[]
#vector,result=TAAC(file)

def SPAC(file): 
    spd_data=read_spd_file(file)
    header=[]
    for f in range(30):
        header.append('SPAC.'+str(f))
    result=[]
    result.append(header)
    vec=[]
    L_sequence=0
    feature = []
    for z in range(len(spd_data)):
        seq=spd_data[z]
        if len(seq)==0:
           break
        L_sequence+=1
        seq=seq.split()
        spd_a=seq[3:]   
        feature.append(spd_a)
    probability_matrix = [[0 for x in range(3)] for y in range(10)]
    temp_feature = feature
    for k in range(0, 10):
        for j in range(5, 8):
            for i in range(0, L_sequence - 1 - k):
                probability_matrix[k][j - 5] += np.float(temp_feature[i][j]) * np.float(temp_feature[i + k][j])
                probability_matrix[k][j - 5] = probability_matrix[k][j - 5] / (L_sequence - 1)
    probability_matrix=np.array(probability_matrix)
    vector=np.reshape(probability_matrix,(1,probability_matrix.shape[0]*probability_matrix.shape[1]))
    vector1=np.array(vector).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return vector,result

#vector,result=[],[]
#vector,result=SPAC(file)

def TAB(file):
    spd_data=read_spd_file(file)
    header=[]
    for f in range(64):
        header.append('TAB.'+str(f))
    result=[]
    result.append(header)
    vec=[]
    L_sequence=0
    i = 0
    feature = []
    for z in range(len(spd_data)):
        seq=spd_data[z]
        if len(seq)==0:
           break
        L_sequence+=1
        seq=seq.split()
        spd_a=seq[3:]   
        feature.append(spd_a)      
    array = feature
    degree_matrix = [[0 for x in range(8)] for y in range(L_sequence-1)]
    matrix = [[0 for x in range(8)] for y in range(8)]
    degree_index = 0
    for x in range(0,L_sequence-1):
        degree_index = 0
        for y in range(1,5):
            degree_matrix[x][degree_index] = np.math.sin(float(array[x][y]) * np.pi / 180 )
            degree_matrix[x][degree_index+1] = np.math.cos(float(array[x][y]) * np.pi / 180  )
            degree_index += 2
    array = degree_matrix
    L_sequence -= 1
    for k in range(0, 8):
        for l in range(0, 8):
            for i in range(0,L_sequence - 1 ):
                matrix[k][l] += float(array[i][k]) * float(array[i + 1][l])
                matrix[k][l] = matrix[k][l] / (L_sequence - 1)
    matrix=np.array(matrix)
    vector=np.reshape(matrix,(1,matrix.shape[0]*matrix.shape[1]))
    vector1=np.array(vector).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return vector,result

#vector,result=[],[]
#vector,result=TAB(file)

def SPB(file):
    spd_data=read_spd_file(file)
    header=[]
    for f in range(9):
        header.append('SPB.'+str(f))
    result=[]
    result.append(header)
    vec=[]
    L_sequence=0
    i = 0
    feature = []
    for z in range(len(spd_data)):
        seq=spd_data[z]
        if len(seq)==0:
           break
        L_sequence+=1
        seq=seq.split()
        spd_a=seq[3:]   
        feature.append(spd_a)
    probability_matrix = [[0 for x in range(3)] for y in range(3)]      
    array = feature
    for k in range(0, 3):
        for l in range(5, 8):
            for i in range(0, L_sequence- 1 ):
                probability_matrix[k][l- 5] += np.float(array[i][k]) * np.float(array[i + 1][l -5 ])
                probability_matrix[k][l - 5] = probability_matrix[k][l - 5] / (L_sequence - 1)
    probability_matrix=np.array(probability_matrix)
    vector=np.reshape(probability_matrix,(1,probability_matrix.shape[0]*probability_matrix.shape[1]))
    vector1=np.array(vector).astype(float)
    vector2=vector1.T
    vector3=vector2[:,0]
    for v in vector3:
        vec.append(v)
    result.append(vec)
    return vector,result

#vector,result=[],[]
#vector,result=SPB(file)
    
# find the path
Script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
# choose the method
option = sys.argv[1]
# the input sequence
file = sys.argv[2]

if(option == "1"):

     vector=AAC(file)
     csv_data=pd.DataFrame(data=vector)
     csv_data.to_csv('AAC_out.csv',header=False,index=False)
elif(option == "2"):

    vector=DC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('DC_out.csv',header=False,index=False)
    
elif(option == "3"):

    vector=CKSAAP(file,gap=2)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('CKSAAP_out.csv',header=False,index=False)
    
elif(option == "4"):

    vector=GDC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('GDC_out.csv',header=False,index=False)
    
elif(option == "5"):

    vector=GTC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('GTC_out.csv',header=False,index=False)
    
elif(option == "6"):

    vector=KSCTriad(file,g=0)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('CTriad_out.csv',header=False,index=False)
    
elif(option == "7"):

    vector=KSCTriad(file,g=1)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('KSCTriad_out.csv',header=False,index=False)
    
elif(option == "8"):

    vector=CTDC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('CTDC_out.csv',header=False,index=False)
    
elif(option == "9"):

    vector=CTDT(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('CTDT_out.csv',header=False,index=False)
    
elif(option == "10"):

    vector=CTDD(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('CTDD_out.csv',header=False,index=False)
    
elif(option == "11"):

    vector=EBGW(file, L=5)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('EBGW_out.csv',header=False,index=False)
    
elif(option == "12"):

    vector=AC(file, nlag = 30)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('AC_out.csv',header=False,index=False)
    
elif(option == "13"):

    vector=NMBroto(file, nlag = 30)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('NMBroto_out.csv',header=False,index=False)
    
elif(option == "14"):

    vector=Moran(file,nlag = 30)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('Moran_out.csv',header=False,index=False)
    
elif(option == "15"):

    vector=Geary(file,nlag = 30)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('Geary_out.csv',header=False,index=False)
    
elif(option == "16"):
 
    vector=QSOrder(file,lambdaValue=30, w=0.05)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('QSOrder_out.csv',header=False,index=False)
    
elif(option == "17"):

    vector=PAAC(file,lambdaValue=30, w=0.05)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('PAAC_out.csv',header=False,index=False)
    
elif(option == "18"):

    vector=APAAC(file,lambdaValue=30, w=0.05)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('APAAC_out.csv',header=False,index=False)
    
elif(option == "19"):
  
    vector=aac_pssm(file,exp=True)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('aac_pssm_out.csv',header=False,index=False)
    
elif(option == "20"):

    vector=dpc_pssm(file,exp=True)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('dpc_pssm_out.csv',header=False,index=False)
    
elif(option == "21"):

    vector=bi_pssm(file,exp=True)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('bi_pssm_out.csv',header=False,index=False)
    
elif(option == "22"):

    vector=ac_pssm(file, GROUP=1,exp=True)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('ac_pssm_out.csv',header=False,index=False)
    
elif(option == "23"):

    vector=cc_pssm(file, GROUP=1,exp=True)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('cc_pssm_out.csv',header=False,index=False)
    
elif(option == "24"):

    vector=acc_pssm(file,GROUP=1,exp=True)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('acc_pssm_out.csv',header=False,index=False)
    
elif(option == "25"):

    vector=pse_pssm(file,ALPHA=1,exp=True)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('pse_pssm_out.csv',header=False,index=False)
    
elif(option == "26"):

    vector=ab_pssm(file,exp=True)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('ab_pssm_out.csv',header=False,index=False)
    
elif(option == "27"):

    vector=SSC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('SSC_out.csv',header=False,index=False)
    
elif(option == "28"):

    vector=ASA(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('ASA_out.csv',header=False,index=False)
    
elif(option == "29"):

    vector=TAC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('TAC_out.csv',header=False,index=False)
    
elif(option == "30"):

    vector=TAB(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('TAB_out.csv',header=False,index=False)
    
elif(option == "31"):

    vector=SPB(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('SPB_out.csv',header=False,index=False)
    
elif(option == "32"):

    vector=TAAC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('TAAC_out.csv',header=False,index=False)
    
elif(option == "33"):

    vector=SPAC(file)
    csv_data=pd.DataFrame(data=vector)
    csv_data.to_csv('SPAC_out.csv',header=False,index=False)
    
else:
    print("Invalid method number. Please check the method table!")



