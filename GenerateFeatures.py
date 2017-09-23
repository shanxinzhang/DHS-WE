# -*- coding: utf-8 -*-
#author: Shanxin Zhang

import numpy as np
from numpy import array
from itertools import combinations, combinations_with_replacement, permutations
from repDNA.nac import RevcKmer, Kmer
from repDNA.psenac import PCPseDNC,PCPseTNC,SCPseDNC,SCPseTNC
from repDNA.ac import DACC,TACC
from repDNA.util import normalize_index,get_data
import multiprocessing
import time
import sys


def GetSequences(f,alphabet):
    #return get_data(f,alphabet=alphabet)
    return get_data(f)

def GetKmerDict(alphabet,k):
    kmerlst=[]
    partkmers=list(combinations_with_replacement(alphabet, k))
    for element in partkmers:
        elelst=set(permutations(element, k))
        strlst=[''.join(ele) for ele in elelst]
        kmerlst+=strlst
    kmerlst=np.sort(kmerlst)
    kmerdict={kmerlst[i]:i for i in range(len(kmerlst))}
    return kmerdict


############################### Mismatch Profile ##############################
def GetMismatchProfile(alphabet, k, m,samples_file):
    kmerdict=GetKmerDict(alphabet, k)
    X=[]
    fp=open(samples_file,'r')
    posi_seq=GetSequences(fp,'ACGT')
    for sequence in array(posi_seq):
        vector=GetMismatchProfileVector(sequence, alphabet, kmerdict, k)
        X.append(vector)  
    X=array(X)
    np.savetxt('DHSs_MismatchProfile_'+str(k)+'_1.txt',X)   
    #return X

def GetMismatchProfileVector(sequence, alphabet, kmerdict, k):    
    vector=np.zeros((1,len(kmerdict)))
    n=len(sequence)
    for i in range(n-k+1):
        subsequence=sequence[i:i+k]
        position=kmerdict.get(subsequence)
        vector[0,position]+=1
        for j in range(k):
            substitution=subsequence
            for letter in set(alphabet)^set(subsequence[j]):
                substitution=list(substitution)
                substitution[j]=letter
                substitution=''.join(substitution)
                position=kmerdict.get(substitution)
                vector[0,position]+=1    
    return list(vector[0])

############################# Subsequence Profile ############################# 
def GetSubsequenceProfileByParallel(alphabet, k, delta,samples_file):
    fp=open(samples_file,'r')
    posi_seq=GetSequences(fp,'ACGT')
    instances=array(posi_seq)
    cpu_num=multiprocessing.cpu_count()   
    batches=ConstructPartitions(instances, cpu_num)
    pool = multiprocessing.Pool(cpu_num)
    results=[]
    for batch in batches:
        temp=pool.apply_async(GetSubsequenceProfile, (batch, alphabet, k, delta))
        results.append(temp)
    pool.close()
    pool.join()
    i=1
    for temp in results:
        temp_X=temp.get()
        if i==1:
            X=temp_X
        else:
            X=np.vstack((X,temp_X))
        i+=1
    X=array(X)
    np.savetxt('DHSs_SubsequenceProfile_'+str(k)+'_'+str(delta)+'.txt',X)
##    return X
    
def ConstructPartitions(instances, cpu_num):
    seqs_num=len(instances)
    batch_num=seqs_num//cpu_num
    batches=[]
    for i in range(cpu_num-1):
        batch=instances[i*batch_num:(i+1)*batch_num]
        batches.append(batch)
    batch=instances[(cpu_num-1)*batch_num:]
    batches.append(batch)
    return batches
    
def GetSubsequenceProfile(instances, alphabet, k, delta):
    kmerdict=GetKmerDict(alphabet, k)
    X=[]
    for sequence in instances:
        vector=GetSubsequenceProfileVector(sequence, kmerdict, k, delta)
        X.append(vector)
    X=array(X)    
    return X

def GetSubsequenceProfileVector(sequence, kmerdict, k, delta):      
    vector=np.zeros((1,len(kmerdict)))
    sequence=array(list(sequence))
    n=len(sequence)
    index_lst=list(combinations(range(n), k))
    for subseq_index in index_lst:
        subseq_index=list(subseq_index)
        subsequence=sequence[subseq_index]
        position=kmerdict.get(''.join(subsequence))     
        subseq_length=subseq_index[-1] - subseq_index[0] + 1
        subseq_score=1 if subseq_length==k else delta**subseq_length    
        vector[0,position]+=subseq_score
    return list(vector[0])
############################### Kmer  ########################################
def GetSpectrumProfile(k,samples_file):
    kmer = Kmer(k=k,normalize=True)
    vec=kmer.make_kmer_vec(open(samples_file))
    np.savetxt('DHSs_kmer_'+str(k)+'.txt',vec)

##    X = array(vec + neg_vec)
##    return X
########################### Reverse Compliment Kmer ###########################
def GetRevcKmer(k,samples_file,):
    rev_kmer = RevcKmer(k=k,normalize=True)
    vec = rev_kmer.make_revckmer_vec(open(samples_file))
    np.savetxt('DHSs_reckmer_'+str(k)+'.txt',vec) 
##    X = array(vec + neg_vec)
##    return X
############ Parallel Correlation Pseudo Dinucleotide Composition #############
def GetPCPseDNC(lamada,w, samples_file,):
    pc_psednc = PCPseDNC(lamada=lamada, w=w)
    phyche_index = [[0.04,0.06,0.04,0.05,0.04,0.04,0.04,0.04,0.05,0.05,0.04,0.06,0.03,0.05,0.04,0.04],
                    [0.08,0.07,0.06,0.10,0.06,0.06,0.06,0.06,0.07,0.07,0.06,0.07,0.07,0.07,0.06,0.08],
                    [0.07,0.06,0.06,0.07,0.05,0.06,0.05,0.05,0.06,0.06,0.06,0.06,0.05,0.06,0.05,0.07],
                    [6.69,6.80,3.47,9.61,2.00,2.99,2.71,3.47,4.27,4.21,2.99,6.80,1.85,4.27,2.00,6.6],
                    [6.24,2.91,2.80,4.66,2.88,2.67,3.02,2.80,3.58,2.66,2.67,2.91,4.11,3.58,2.88,6.24],
                    [21.34,21.98,17.48,24.79,14.51,14.25,14.66,17.48,18.41,17.31,14.25,21.98,14.24,18.41,14.51,21.34],
                    [1.05,2.01,3.60,0.61,5.60,4.68,6.02,3.60,2.44,1.70,4.68,2.01,3.50,2.44,5.60,1.05],
                    [-1.26,0.33,-1.66,0.00,0.14,-0.77,0.00,-1.66,1.44,0.00,-0.77,0.33,0.00,1.44,0.14,-1.26],
                    [35.02,31.53,32.29,30.72,35.43,33.54,33.67,32.29,35.67,34.07,33.54,31.53,36.94,35.67,35.43,35.02],
                    [-0.18,-0.59,-0.22,-0.68,0.48,-0.17,0.44,-0.22,-0.05,-0.19,-0.17,-0.59,0.04,-0.05,0.48,-0.18],
                    [0.01,-0.02,-0.02,0.00,0.01,0.03,0.00,-0.02,-0.01,0.00,0.03,-0.02,0.00,-0.01,0.01,0.01],
                    [3.25,3.24,3.32,3.21,3.37,3.36,3.29,3.32,3.30,3.27,3.36,3.24,3.39,3.30,3.37,3.25],
                    [-1.00,-1.44,-1.28,-0.88,-1.45,-1.84,-2.17,-1.28,-1.30,-2.24,-1.84,-1.44,-0.58,-1.30,-1.45,-1.00],
                    [-7.60,-8.40,-7.80,-7.20,-8.50,-8.00,-10.60,-7.80,-8.20,-9.80,-8.00,-8.40,-7.20,-8.20,-8.50,-7.60],
                    [-21.30,-22.40,-21.00,-20.40,-22.70,-19.90,-27.20,-21.00,-22.20,-24.40,-19.90,-22.40,-21.30,-22.20,-22.70,-21.30]]
    vec = pc_psednc.make_pcpsednc_vec(open(samples_file),extra_phyche_index=normalize_index(phyche_index,is_convert_dict=True))
    np.savetxt('DHSs_PC_PseDNC_'+str(lamada)+'_'+str(w)+'.txt',vec)    
##    X = array(vec + neg_vec)    
##    return X

############ Parallel Correlation Pseudo Trinucleotide Composition ############
def GetPCPseTNC(lamada,w,samples_file,):
    pc_psetnc = PCPseTNC(lamada=lamada, w=w)
    vec = pc_psetnc.make_pcpsetnc_vec(open(samples_file), all_property=True)
    np.savetxt('DHSs_PC_PseTNC_'+str(lamada)+'_'+str(w)+'.txt',vec) 
##    X = array(vec + neg_vec)
##    return X
    
############## Series Correlation Pseudo Dinucleotide Composition #############
def GetSCPseDNC(lamada, w, samples_file,):
    sc_psednc = SCPseDNC(lamada=lamada, w=w)
    phyche_index = [[0.04,0.06,0.04,0.05,0.04,0.04,0.04,0.04,0.05,0.05,0.04,0.06,0.03,0.05,0.04,0.04],
                    [0.08,0.07,0.06,0.10,0.06,0.06,0.06,0.06,0.07,0.07,0.06,0.07,0.07,0.07,0.06,0.08],
                    [0.07,0.06,0.06,0.07,0.05,0.06,0.05,0.05,0.06,0.06,0.06,0.06,0.05,0.06,0.05,0.07],
                    [6.69,6.80,3.47,9.61,2.00,2.99,2.71,3.47,4.27,4.21,2.99,6.80,1.85,4.27,2.00,6.6],
                    [6.24,2.91,2.80,4.66,2.88,2.67,3.02,2.80,3.58,2.66,2.67,2.91,4.11,3.58,2.88,6.24],
                    [21.34,21.98,17.48,24.79,14.51,14.25,14.66,17.48,18.41,17.31,14.25,21.98,14.24,18.41,14.51,21.34],
                    [1.05,2.01,3.60,0.61,5.60,4.68,6.02,3.60,2.44,1.70,4.68,2.01,3.50,2.44,5.60,1.05],
                    [-1.26,0.33,-1.66,0.00,0.14,-0.77,0.00,-1.66,1.44,0.00,-0.77,0.33,0.00,1.44,0.14,-1.26],
                    [35.02,31.53,32.29,30.72,35.43,33.54,33.67,32.29,35.67,34.07,33.54,31.53,36.94,35.67,35.43,35.02],
                    [-0.18,-0.59,-0.22,-0.68,0.48,-0.17,0.44,-0.22,-0.05,-0.19,-0.17,-0.59,0.04,-0.05,0.48,-0.18],
                    [0.01,-0.02,-0.02,0.00,0.01,0.03,0.00,-0.02,-0.01,0.00,0.03,-0.02,0.00,-0.01,0.01,0.01],
                    [3.25,3.24,3.32,3.21,3.37,3.36,3.29,3.32,3.30,3.27,3.36,3.24,3.39,3.30,3.37,3.25],
                    [-1.00,-1.44,-1.28,-0.88,-1.45,-1.84,-2.17,-1.28,-1.30,-2.24,-1.84,-1.44,-0.58,-1.30,-1.45,-1.00],
                    [-7.60,-8.40,-7.80,-7.20,-8.50,-8.00,-10.60,-7.80,-8.20,-9.80,-8.00,-8.40,-7.20,-8.20,-8.50,-7.60],
                    [-21.30,-22.40,-21.00,-20.40,-22.70,-19.90,-27.20,-21.00,-22.20,-24.40,-19.90,-22.40,-21.30,-22.20,-22.70,-21.30]]
    vec = sc_psednc.make_scpsednc_vec(open(samples_file), extra_phyche_index=normalize_index(phyche_index,is_convert_dict=True))
    np.savetxt('DHSs_SC_PseDNC_'+str(lamada)+'_'+str(w)+'.txt',vec)   
##    X = array(vec + neg_vec)
##    return X  
    
############## Series Correlation Pseudo Trinucleotide Composition ############
def GetSCPseTNC(lamada,w,samples_file,):
    sc_psetnc = SCPseTNC(lamada=lamada,w=w)
    vec = sc_psetnc.make_scpsetnc_vec(open(samples_file), all_property=True)
    np.savetxt('DHSs_SC_PseTNC_'+str(lamada)+'_'+str(w)+'.txt',vec)   
##    X = array(vec + neg_vec)
##    return X

#####################  Dinucleotide-based auto-cross covariance   ###############################################
def GetDACC(lag,samples_file,):
    dacc=DACC(lag)
    phyche_index = [[0.04,0.06,0.04,0.05,0.04,0.04,0.04,0.04,0.05,0.05,0.04,0.06,0.03,0.05,0.04,0.04],
                    [0.08,0.07,0.06,0.10,0.06,0.06,0.06,0.06,0.07,0.07,0.06,0.07,0.07,0.07,0.06,0.08],
                    [0.07,0.06,0.06,0.07,0.05,0.06,0.05,0.05,0.06,0.06,0.06,0.06,0.05,0.06,0.05,0.07],
                    [6.69,6.80,3.47,9.61,2.00,2.99,2.71,3.47,4.27,4.21,2.99,6.80,1.85,4.27,2.00,6.6],
                    [6.24,2.91,2.80,4.66,2.88,2.67,3.02,2.80,3.58,2.66,2.67,2.91,4.11,3.58,2.88,6.24],
                    [21.34,21.98,17.48,24.79,14.51,14.25,14.66,17.48,18.41,17.31,14.25,21.98,14.24,18.41,14.51,21.34],
                    [1.05,2.01,3.60,0.61,5.60,4.68,6.02,3.60,2.44,1.70,4.68,2.01,3.50,2.44,5.60,1.05],
                    [-1.26,0.33,-1.66,0.00,0.14,-0.77,0.00,-1.66,1.44,0.00,-0.77,0.33,0.00,1.44,0.14,-1.26],
                    [35.02,31.53,32.29,30.72,35.43,33.54,33.67,32.29,35.67,34.07,33.54,31.53,36.94,35.67,35.43,35.02],
                    [-0.18,-0.59,-0.22,-0.68,0.48,-0.17,0.44,-0.22,-0.05,-0.19,-0.17,-0.59,0.04,-0.05,0.48,-0.18],
                    [0.01,-0.02,-0.02,0.00,0.01,0.03,0.00,-0.02,-0.01,0.00,0.03,-0.02,0.00,-0.01,0.01,0.01],
                    [3.25,3.24,3.32,3.21,3.37,3.36,3.29,3.32,3.30,3.27,3.36,3.24,3.39,3.30,3.37,3.25],
                    [-1.00,-1.44,-1.28,-0.88,-1.45,-1.84,-2.17,-1.28,-1.30,-2.24,-1.84,-1.44,-0.58,-1.30,-1.45,-1.00],
                    [-7.60,-8.40,-7.80,-7.20,-8.50,-8.00,-10.60,-7.80,-8.20,-9.80,-8.00,-8.40,-7.20,-8.20,-8.50,-7.60],
                    [-21.30,-22.40,-21.00,-20.40,-22.70,-19.90,-27.20,-21.00,-22.20,-24.40,-19.90,-22.40,-21.30,-22.20,-22.70,-21.30]]
    vec=dacc.make_dacc_vec(open(samples_file), extra_phyche_index=normalize_index(phyche_index,is_convert_dict=True))
    np.savetxt('DHSs_dacc_'+str(lag)+'.txt',vec) 
##    X = array(vec + neg_vec)
##    return X  

###############################################################################

    
if __name__ == '__main__':

##    samples_file=sys.argv[1]
##    =sys.argv[2]
    
    samples_file=sys.argv[1]
    
    alphabet=['A','C','G','T']
   
    
    # Spectrum Profile for k=2,3,4,5,6
    for k in range(2,7):
        print('..........................................................................')
        print('Coding for feature:'+str(k)+'-Kmer, beginning')
        tic=time.clock()
        GetSpectrumProfile(k,samples_file)
        toc=time.clock()
        print('Coding time:%.3f minutes'%((toc-tic)/60))
        
    # Mismatch Profile for (k,m)=(3,1),(4,1),(5,1),(6,1)
    for (k,m) in [(3,1),(4,1),(5,1),(6,1)]:
        print('..........................................................................')
        print('Coding for feature:'+str((k,m))+'-Mismatch Profile, beginning')
        tic=time.clock()
        GetMismatchProfile(alphabet,k,m,samples_file)
        toc=time.clock()
        print('Coding time:%.3f minutes'%((toc-tic)/60))    

##    # Subsequence Profile for (k,delta)=(3,1),(4,1),(5,1),(6,1)
    for (k,delta) in [(2,0.5),(3,0.7)]:
            print('..........................................................................')
            print('Coding for feature:'+str((k,delta))+'-Subsequence Profile, beginning')
            print('The process may spend some time, please do not close the program')
            tic=time.clock()
            GetSubsequenceProfileByParallel(alphabet,k,delta,samples_file)
            toc=time.clock()
            print('Coding time:%.3f minutes'%((toc-tic)/60)) 
        
##    # Reverse Compliment Kmer for k=1,2,3,4,5,6
    for k in range(2,7):
        print('..........................................................................')
        print('Coding for feature:'+str(k)+'-RevcKmer, beginning')
        tic=time.clock()
        GetRevcKmer(k,samples_file)
        toc=time.clock()
        print('Coding time:%.3f minutes'%((toc-tic)/60))
        
##    # Parallel Correlation Pseudo Dinucleotide Composition   
    print('..........................................................................')
    print('Coding for feature:PCPseDNC, beginning')
    tic=time.clock()
    GetPCPseDNC(2,0.2,samples_file)
    toc=time.clock()
    print('Coding time:%.3f minutes'%((toc-tic)/60))
##
##    # Parallel Correlation Pseudo Trinucleotide Composition   
    print('..........................................................................')
    print('Coding for feature:PCPseTNC, beginning')
    tic=time.clock()
    GetPCPseTNC(6,0.1,samples_file)
    toc=time.clock()
    print('Coding time:%.3f minutes'%((toc-tic)/60))
##    
##    # Series Correlation Pseudo Dinucleotide Composition   
    print('..........................................................................')
    print('Coding for feature:SCPseDNC, beginning')
    tic=time.clock()
    GetSCPseDNC(1,0.1,samples_file)
    toc=time.clock()
    print('Coding time:%.3f minutes'%((toc-tic)/60))
##    
##    # Series Correlation Pseudo Trinucleotide Composition   
    print('..........................................................................')
    print('Coding for feature:SCPseTNC, beginning')
    tic=time.clock()
    GetSCPseTNC(6,0.1,samples_file)
    toc=time.clock()
    print('Coding time:%.3f minutes'%((toc-tic)/60))  
##
##    # Drinucleotide-based auto-cross covariance
    print('..........................................................................')
    print('Coding for feature:DACC, beginning')
    tic=time.clock()
    X=GetDACC(1,samples_file)
    toc=time.clock()
    print('Coding time:%.3f minutes'%((toc-tic)/60)) 



