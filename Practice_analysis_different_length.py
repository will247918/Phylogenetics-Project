

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 11:24:48 2023

@author: willi
"""

from Bio import AlignIO
import numpy as np
from matplotlib import pyplot as plt
import itertools

#Where the nexus file is saved
nexus_file = "C:/Users/willi/BioinformaticsProject/Data/HBV_data.nex"


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
        #Alignment of the sequence data in the NEXUS file
        
        data = AlignIO.read(nexus_file, "nexus") 
        alignment_list=[]
        for x in data:
            x = x.seq.split('-', 1)[0]
            x = x.rstrip()
            alignment_list.append(x)
        self.alignment=alignment_list
        

    def distance_matrix(self):
        '''
        This function, forms the overall distance matrix from the alignment of 
        all sequences.

        Returns
        -------
        distance_matrix : numpy array
            The distance matrix calculates the difference in nucleotides between 
            the different sequences and ouputs them in a matrix form.
            

        '''
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
        '''
        This function, instead of caluclating the distance matrix for all the 
        sequences, calculates the distance matrix for 4 taxa. This is important
        when calculating the delta value between the taxa.

        Parameters
        ----------
        array : numpy array
            This numpy array contains random numbers which are used to identify
            the index of the sequence from its alignment.

        Returns
        -------
        distance_matrix2 : numpy array
            This matrix is a 4x4 matrix, showing the alignment between the 4 taxa
            seqeunces, which can be used to calculate the delta values.

        '''
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
        '''
        This function converts sequences and appends individual nucleotides 
        from the sequence to a list of lists

        Returns
        -------
        empty_array : list of lists
            DESCRIPTION.
            Empty array contains a lists of lists. With each list containing
            nucleotides from that sequence alignment

        '''
        #Create an empty list
        empty_array = []
        #Append the sequences to a list
        for i in range(0,len(list(self.alignment))):
            empty_array.append(self.alignment[i])
        return empty_array
    
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
    
    def delta_value_2(self,array):
        print(f'array is {array}')
        delta_list = []
        for i in array:
            
            individual_matrix = self.different_length_matrix2(i)
            delta_value_matrix = self.delta_value(individual_matrix)
            delta_list.append(delta_value_matrix)
            
        return delta_list
            
            
        
    def random_sample(self):
        '''
        This function is creating random samples each containing one constant
        taxa for the number defined. This allows us to calculate the mean
        delta value for each taxa.

        Returns
        -------
        sample : List of lists
            DESCRIPTION.
            List of lists containing random samples with a constant taxa.

        '''
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
            
            for j in range(0,20):
                
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
        print(sample)
        #Call the delta_values funciton to produce a list of delta values for the samples
        list_delta_values = self.delta_value_2(sample)
        #Create an empty list to store mean delta values for each taxa
        list_mean_delta = []
        
        #Split the delta values into lists contianing each taxa
        array_delta_values = np.array(list_delta_values)
        
        split_list = np.split(array_delta_values,len(self.alignment))
        
       
        
        
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
        plt.ylim(0, 0.5)
        plt.xlabel("Taxa")
        plt.ylabel("Mean delta value")
        plt.show()
        return mean_plot
    
    def overall_delta_plot(self):
        list_mean_delta, array_delta_values = self.delta_random_mean()
        overall_delt_plot = plt.hist(array_delta_values,bins=10)
        plt.show()
        return overall_delt_plot
    
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
        
        
    def all_samples(self):
        s = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
        n = 4        
        subset_array = list(itertools.combinations(s, n))
        
        return subset_array
        
    def delta_plot1(self):
        print(self.alignment[0])
        
        array4 = self.all_samples()
        
        array5=[]
        
        for i in range(0,len(array4)):
            array5.append(list(array4[i]))
        
        delta_values_list=self.delta_value_2(array5)
        
   
        
        print(delta_values_list)
        
           
            
        
        
        overall_delt = plt.hist(delta_values_list,bins=10)
        plt.show()
        return overall_delt
        
    def delta_plot_java(self,file_path):
        delta_value_list=[]
        list_delta_file = open(file_path,'r')
        lines = list_delta_file.readlines()
        
        for line in lines:
            format(line.strip('\n'))
            delta_value_list.append(float(line.strip()))
        
        return delta_value_list
            
    def java_histogram(self,file_path):
        list_delta_values = self.delta_plot_java(file_path)
        
        overall_delt_plot = plt.hist(list_delta_values,bins=10)
        plt.show()
        return overall_delt_plot
        

        

            
    
        

    
sequence = delta_plot(nexus_file)

#list_mean_delta,array_delta_values = sequence.delta_random_mean()

#print(sequence.mean_delta_plot())

file_path = "C:/Users/willi/Documents/YEAR2/Labsheets programming 1/Week6/Delta_plot/Data/output.txt"

print(sequence.java_histogram(file_path))

