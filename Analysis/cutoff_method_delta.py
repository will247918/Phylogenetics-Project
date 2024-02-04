# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 16:08:03 2024

@author: willi
"""



from Bio import AlignIO
import numpy as np
from matplotlib import pyplot as plt
import itertools
from matplotlib.ticker import PercentFormatter


#Where the nexus file is saved
nexus_file = "C:/Users/willi/BioinformaticsProject/Data/HBV_distance_matrix.nex"


class delta_plot():
    
    def __init__(self,nexus_file):
        '''
        Initialise the class, this just takes in and sets the nexus file.

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
        '''
        Extract the distance matrix from the NEXUS file, in the distances block
        to be used in the further analysis of calculating the delta plot.

        Returns
        -------
        distance_matrix2 : Numpy 2D array
            This is a 2D array representing the distance matrix of all the taxa
            taken from the NEXUS file.

        '''
        #Initialise empty list
        distance_matrix =[]
        #Open file
        file = open(nexus_file)
        #Iterate through the file
        for line in file:
            #Only read lines after the word matrix
            if 'matrix' in line:
                for line in file:
                    length_line = len(line)
                    #Split the line to only get the numbers of the matrix
                    line = line.split()
                    distance_matrix.append(line[2:length_line])

        file.close()
        #Add the numbers to a numpy array matrix
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
        
    def delta_value_2(self,matrix):
        '''
        This function takes in a 4x4 matrix, considering the quartet of 4 
        taxa and outputs, the delta value. It calculates the sides of the 
        quadrilateral and using this caluclates the ratio of:
        min(s,l)/max(s,l).

        Parameters
        ----------
        matrix : numpy array
            This is a distance matrix representing the distance between taxa.
            

        Returns
        -------
        float
            Represents the delta value.

        '''
        distance_matrix = matrix
        #Create distance array
        da = list(distance_matrix[np.triu_indices(4)])
        #List of dxy|uv
        l_1 =[da[1],da[8]]
        #List of dxu|yv
        l_2 =[da[2],da[6]]
        #List of dxv|yu
        l_3 =[da[3],da[5]]
        
        list_of_lists = [l_1,l_2,l_3]
        p = sorted(list_of_lists, key=lambda x: sum(x))

       
        
        
        
        list_1 =[da[3],da[5],da[1],da[8],da[2],da[6]]
        list_1.sort()
     
        x= (p[2][0]+p[2][1]-p[0][0]-p[0][1])/2
        
        y= (p[2][0]+p[2][1]-p[1][0]-p[1][1]) /2
        
        s = min(x,y)
        
        l = max(x,y)
        
        return s/l,s,l
  
    def different_length_matrix2(self,array):
        '''
        This function takes in a 4 membered array, containing 4 taxa. For example
        this could be taxa 4,7,13,20 and outputs the 4x4 distance matrix of 
        these taxa.

        Parameters
        ----------
        array : Numpy array
            This is array with 4 numbers representing the index of the taxa.

        Returns
        -------
        distance_matrix2 : TYPE
            DESCRIPTION.

        '''
        distance_matrix2 = np.empty((4,4))
        array2=[]
      
        #Convert array values into integers
        
        for i in array:
            i =int(i)
            array2.append(i)
        
        #Compare the sequences, with one sequence at a time being constant.
        #If there is a difference in a nucleotide you add one to the difference.
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
        '''
        

        Returns
        -------
        overall_delt : TYPE
            DESCRIPTION.

        '''
        
        delta_values_list = self.delta_value_calculator()
        

        overall_delt = plt.hist(delta_values_list,bins=10,density=True)
        plt.gca().yaxis.set_major_formatter(PercentFormatter(10))
        plt.ylabel("Proportion of total numbers")
        plt.xlabel("Delta values")
        plt.show()
        return overall_delt
    
    def delta_value_calculator(self):
        '''
        This function claculates the delta value for each quartet in the 
        distance matrix.

        Returns
        -------
        Array
            Returns an array of all the delta values for all quartets in the 
            distance matrix.

        '''
        #Define a list to store the delta values
        delta_list=[]
        
        original_matrix = self.distance_matrix()
        #Define number of taxa
        s = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
        #Use itertools to go through each quartet combination of taxa to find
        #the delta value.
        for quartet in itertools.combinations(s,4):
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            delt_value = self.delta_value_2(quartet_matrix)[0]
            delta_list.append(delt_value)
        print(np.mean(np.array(delta_list)))
        return np.array(delta_list)
        
    def average_length_quartet(self):
        '''
        This method incorporates the 'cut-off' method. If the delta value is 
        above a certain value, these quartets are isolated. Then if the area 
        of the quartet is also below a certain value the number of star like 
        quartets can be identified.

        Returns
        -------
        This outputs the number of quartets and th emean length, width and area 
        of the quartets.

        '''
        #Initiate 2 empty lists to store length values
        l1=[]
        l2=[]
        area_list = []
        original_matrix = self.distance_matrix()
        print(original_matrix)
        s = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
        for quartet in itertools.combinations(s,4):
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            if self.delta_value_2(quartet_matrix)[0] > 0.8:
                s = self.delta_value_2(quartet_matrix)[1]
                l = self.delta_value_2(quartet_matrix)[2]
                area=s*l
                if area<0.00001:
                    print([s,l,area])
                    #append the length values to the empty lists
                    l1.append(s)
                    l2.append(l)
                    area_list.append(area)
                
        #print(l2)
        #Convert the lists into arrays
        print(len(area_list))
        l1 = np.array(l1)
        l2=np.array(l2)
        area_list = np.array(area_list)

        return [l1.mean(),l2.mean(),area_list.mean()]

    
sequence = delta_plot(nexus_file)

#list_mean_delta,array_delta_values = sequence.delta_random_mean()

#print(sequence.mean_delta_plot())

file_path = "C:/Users/willi/Documents/YEAR2/Labsheets programming 1/Week6/Delta_plot/Data/output.txt"

#sequence.second_try()
print(sequence.delta_plot1())
print(sequence.average_length_quartet())