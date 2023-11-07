# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 10:24:06 2023

@author: willi
"""

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
from matplotlib import pyplot as plt

#Where the nexus file is saved
nexus_file = "C:/Users/willi/BioinformaticsProject/Data/ExampleData.nex"


class delta_plot():
    
    def __init__(self,nexus_file):
        #Alignment of the sequence data in the NEXUS file
        file = open("C:/Users/willi/BioinformaticsProject/Data/practice_data.txt") 
        self.alignment=[]
        for line in file:
            line = line.split('#', 1)[0]
            line = line.rstrip()
            self.alignment.append(line)
        
        
        file.close()

    def distance_matrix(self):
        self.alignment = [list(string) for string in self.alignment]
        #Create an empty matrix
        distance_matrix = np.empty(((len(self.alignment)),len(self.alignment)))
        #print(distance_matrix)
        
        #Series of loops to calculate the difference between sequences in an alignment.
        for i in range(0,len(self.alignment)):
            for k in range(0,len(self.alignment)):
                different_nucleotides=0 
                #Work out length of each seqeunce you are comparing
                len_first = len(self.alignment[i])
                len_second= len(self.alignment[k])
                
                #Work out the minimum number of nucleotides in either sequence
                min_nucleotides = min(len_first,len_second)
                #Add nucleotides if the seqeunces differ in size
                if len_first>len_second:
                    different_nucleotides = different_nucleotides+ (len_first-len_second)
                elif len_second>len_first:
                    different_nucleotides = different_nucleotides + (len_second-len_first)
                
                #For loop to compare nucelotides in different sequences
                for j in range(0,min_nucleotides): 
                   if self.alignment[i][j] != self.alignment[k][j]:
                       different_nucleotides=different_nucleotides+1
                       
                #Add the number of different nucleotides to the distance matrix
                distance_matrix[i][k]=different_nucleotides                   
                
    
        return distance_matrix
    
    #Input any 4 indexed array and will output the distance matrix
    def distance_matrix2(self,array):
        distance_matrix2 = np.empty((4,4))
        array2=[]
        
        
        #Convert array values into integers
        for i in array:
            i =int(i)
            array2.append(i)
        
        #Compare the sequences, with one sequence at a time being constant.
        for i,j in enumerate(array2):
            for k in range(0,4):
                different_nucleotides=0
                for m in range(0,20): 
                   if self.alignment[j].seq[m] != self.alignment[array2[k]].seq[m]:
                       different_nucleotides=different_nucleotides+1
                
                distance_matrix2[i][k]=different_nucleotides
            
                
        return distance_matrix2    
            
            
        
    
    def array(self):
        empty_array = []
        for i in range(0,len(list(self.alignment))):
            empty_array.append(self.alignment[i])
        return empty_array
    
    def delta_value(self,matrix):
        distance_matrix = matrix
        distance_array = list(distance_matrix[np.triu_indices(4)])
        

        d1 = distance_array[1]+distance_array[8]
        d2 = distance_array[5]+distance_array[3]
        d3 = distance_array[6]+distance_array[2]
        distance_list = [d1,d2,d3]
        distance_list.sort()

        delta_value = (distance_list[2]-distance_list[1])/(distance_list[2]-distance_list[0])
        return delta_value
    
    def delta_value_2(self,array):
        
        delta_list = []
        for i in array:
            individual_matrix = self.different_length_matrix2(i)
            delta_value_matrix = self.delta_value(individual_matrix)
            delta_list.append(delta_value_matrix)
            
        return delta_list
            
            
        
    def random_sample(self):
        # This code generates random indexes, from 0 to 4 and pastes them in a list.
        # 5 samples are taken for each index
        list_sequences = self.array()
        
        #Make an empty samples list to store all the samples.
        sample=[]
        
        #For loop to run over from 0 to the number of taxa
        for i in range(0,len(list_sequences)):
            #List1 is a list of indexes 0 to 4
            list1 = list(range(len(list_sequences)))
            #Remove index
            list1.pop(i)
            #Another for loop to create the random choice-List of quartets
            # With 5 quartets for each taxa - can be changed.
            
            for j in range(0,500):
                
                random_choice = list(np.random.choice(list1,3))
                
                while True:
                    random_choice.sort()
                    if random_choice == list(set(random_choice)):
                        random_choice.append(i)
                        random_choice.sort()
                        break
                    else:
                        random_choice=list(np.random.choice(list1,3))
                        
                       
                
                    
                sample.append(random_choice)
                    
                j=+1
        
        return sample
  
    #This function returns a list of mean delta values for each taxa
    def delta_random_mean(self):
        #Call the random_smaple function to generate random lists of taxa
        sample= self.random_sample()
        #Call the delta_values funciton to produce a list of delta values for the samples
        list_delta_values = self.delta_value_2(sample)
        #Create an empty list to store mean delta values for each taxa
        list_mean_delta = []
        
        #Split the delta values into lists contianing each taxa
        array_delta_values = np.array(list_delta_values)
        
        split_list = np.split(array_delta_values,5)
        
       
        
        
        for i in split_list:
            mean =np.mean(i)
            list_mean_delta.append(mean)
        return list_mean_delta,array_delta_values
            
    def mean_delta_plot(self):
        mean_delta_array,list_delta = self.delta_random_mean()
        print(mean_delta_array)
        list_taxa = []
        for i in range(0,len(self.alignment)):
            list_taxa.append(i+1)
        
        mean_plot = plt.bar(list_taxa,mean_delta_array,width =0.2)
        plt.ylim(0, 0.3)
        plt.xlabel("Taxa")
        plt.ylabel("Mean delta value")
        plt.show()
        return mean_plot
    
    def overall_delta_plot(self):
        list_mean_delta, array_delta_values = self.delta_random_mean()
        plt.hist(array_delta_values,bins=10)
        plt.show()
        print(array_delta_values)
    
    def different_length_matrix2(self,array):
        distance_matrix2 = np.empty((4,4))
        array2=[]
        
        #Convert array values into integers
        for i in array:
            i =int(i)
            array2.append(i)
        
        #Compare the sequences, with one sequence at a time being constant.
        for i,j in enumerate(array2):
            for k in range(0,4):
                different_nucleotides=0
                if(len(self.alignment[j].seq)>len(self.alignment[array2[k]])):
                    different_nucleotides = different_nucleotides + (len(self.alignment[j].seq)-len(self.alignment[array2[k]]))
                elif (len(self.alignment[array2[k]]))>len(self.alignment[j].seq):
                    different_nucleotides = different_nucleotides + (len(self.alignment[array2[k]].seq)-len(self.alignment[j].seq))
                
                min_range = min(len(self.alignment[array2[k]]),len(self.alignment[j].seq))
                for m in range(0,min_range): 
                   if self.alignment[j].seq[m] != self.alignment[array2[k]].seq[m]:
                       different_nucleotides=different_nucleotides+1
                
                distance_matrix2[i][k]=different_nucleotides
            
                
        return distance_matrix2   
        
        
        
            
        
        
        

sequence = delta_plot(nexus_file)     
print(sequence.distance_matrix())
        

    

# 
# print(sequence.delta_random_mean())
# 
# #print(sequence.alignment[1].seq[1])
# print(sequence.overall_delta_plot())
# =============================================================================

