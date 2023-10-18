# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 11:24:48 2023

@author: willi
"""

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align import PairwiseAligner
from Bio import SeqIO
import numpy as np
import itertools
import pandas as pd
import random

#Where the nexus file is saved
nexus_file = "C:/Users/willi/BioinformaticsProject/Data/ExampleData.nex"


class delta_plot():
    
    def __init__(self,nexus_file):
        #Alignment of the sequence data in the NEXUS file
        self.alignment = AlignIO.read(nexus_file, "nexus")  

    def distance_matrix(self):
        #Create an empty matrix
        distance_matrix = np.empty((4,4))
        #print(distance_matrix)
        
        #Series of loops to calculate the difference between sequences in an alignment.
        for i in range(0,len(list(self.alignment))):
            for k in range(0,len(list(self.alignment))):
                different_nucleotides=0 
                for j in range(0,self.alignment.get_alignment_length()): 
                   if self.alignment[i].seq[j] != self.alignment[k].seq[j]:
                       different_nucleotides=different_nucleotides+1
                       
                distance_matrix[i][k]=different_nucleotides
                    
        return distance_matrix
    
    def array(self):
        empty_array = []
        for i in range(0,4):
            empty_array.append(self.alignment[i])
        return empty_array
    
    def delta_value(self):
        distance_matrix = sequence.distance_matrix()
        distance_array = list(distance_matrix[np.triu_indices(4)])
        print(distance_array)

        d1 = distance_array[1]+distance_array[8]
        d2 = distance_array[5]+distance_array[3]
        d3 = distance_array[6]+distance_array[2]
        distance_list = [d1,d2,d3]
        distance_list.sort()

        delta_value = (distance_list[2]-distance_list[1])/(distance_list[2]-distance_list[0])
        return delta_value
    
    def random_sample(self):
        # This code generates random indexes, from 0 to 4 and pastes them in a list.
        # 5 samples are taken for each index
        list_sequences = self.array()
        print(list_sequences)
        sample=[]
        for i in range(0,len(list_sequences)):
            list1 = list(range(len(list_sequences)))
            list1.pop(i)
            for j in range(0,len(list_sequences)+1):
                sample.append(np.random.choice(list1,4))
                j=+1
        
        return sample
  
  
        
  
            
    
        
        

    
sequence = delta_plot(nexus_file)
print(sequence.distance_matrix())  
empty_list=[]
print(sequence.random_sample())



