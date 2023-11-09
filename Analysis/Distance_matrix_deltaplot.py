# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 15:08:46 2023

@author: willi
"""

from Bio import AlignIO
import numpy as np
from matplotlib import pyplot as plt
import itertools

#Where the nexus file is saved
nexus_file = "C:/Users/willi/BioinformaticsProject/Data/HBV_distance_matrix.nex"


class delta_plot():
    
    def __init__(self,nexus_file):
        '''
        Initialise the class, this outputs the alignment of the sequences by 
        utilising the nexus file format.

        Parameters
        ----------
        nexus_file : nex
            This function recieves the nexus file as an input, which is a file 
            type commonly used in bioinformatics analysis.

        Returns
        -------
        None.

        '''
        
        self.nexus_file=nexus_file
        

    def distance_matrix(self):
        distance_matrix =[]
        file = open(nexus_file)
        for line in file:
            if 'matrix' in line:
                for line in file:
                    length_line = len(line)
                    line = line.split()
                    distance_matrix.append(line[2:length_line])

        file.close()
        distance_matrix2=np.zeros((24,24))
        for i,k in enumerate (distance_matrix):
            if len(k)==0:
                break
            for j in range(0,24):
                distance_value = k[j]
                distance_matrix2[i][j]=distance_value
                
                

            
            
        return distance_matrix2

            
            
    
    def delta_value(self,matrix):
        '''
        This function is used to calculate the delta value from a 4x4 matrix
        caluclated using the distance_matrix_2 matrix. This works by splitting
        the matrix into a list and indexing the list to calculate the delta value.

        Parameters
        ----------
        matrix : numpy array, 4x4
            This matrix is inputted 

        Returns
        -------
        delta_value : TYPE
            DESCRIPTION.

        '''
        distance_matrix = matrix
        distance_array = list(distance_matrix[np.triu_indices(4)])
        

        d1 = distance_array[1]+distance_array[8]
        d2 = distance_array[5]+distance_array[3]
        d3 = distance_array[6]+distance_array[2]
        distance_list = [d1,d2,d3]
        distance_list.sort()
        

        delta_value = (distance_list[2]-distance_list[1])/(distance_list[2]-distance_list[0])
        #print(f'Delta value is {distance_list}')
        return delta_value
    
  
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
                
                if(len(self.alignment[j])>len(self.alignment[array2[k]])):
                    different_nucleotides = different_nucleotides + (len(self.alignment[j])-len(self.alignment[array2[k]]))
                elif (len(self.alignment[array2[k]]))>len(self.alignment[j]):
                   different_nucleotides = different_nucleotides + (len(self.alignment[array2[k]])-len(self.alignment[j]))
                
                min_range = min(len(self.alignment[array2[k]]),len(self.alignment[j]))
                for m in range(0,min_range): 
                   if self.alignment[j][m] != self.alignment[array2[k]][m]:
                       different_nucleotides=different_nucleotides+1
                
                distance_matrix2[i][k]=different_nucleotides
            
        
        return distance_matrix2   
        
        
   
    def delta_plot1(self):
        
        delta_values_list = self.second_try()
        
        overall_delt = plt.hist(delta_values_list,bins=10)
        plt.show()
        return overall_delt
    
    def delta_value_calculator(self):
        delta_list=[]
        original_matrix = self.distance_matrix()
        s = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
        for quartet in itertools.combinations(s,4):
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            delt_value = self.delta_value(quartet_matrix)
            delta_list.append(delt_value)
        return delta_list
        
    
        

    
sequence = delta_plot(nexus_file)

#list_mean_delta,array_delta_values = sequence.delta_random_mean()

#print(sequence.mean_delta_plot())

file_path = "C:/Users/willi/Documents/YEAR2/Labsheets programming 1/Week6/Delta_plot/Data/output.txt"

#sequence.second_try()
print(sequence.delta_plot1())
