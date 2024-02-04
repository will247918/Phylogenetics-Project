# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 11:58:12 2024

@author: willi
"""


from Bio import AlignIO
import numpy as np
from matplotlib import pyplot as plt
import itertools
from matplotlib.ticker import PercentFormatter
import matplotlib.ticker as mtick


#Where the nexus file is saved
nexus_file = "C:/Users/willi/BioinformaticsProject/Data/AB.nex"


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
        number_taxa=0
        #Iterate through the file
        for line in file:
            #Only read lines after the word matrix
            if 'matrix' in line:
                for line in file:
                    number_taxa+=1
                    length_line = len(line)+2
                    #Split the line to only get the numbers of the matrix
                    line = line.split()
                    distance_matrix.append(line[2:length_line])
        print(number_taxa)
        file.close()
        #Add the numbers to a numpy array matrix
        distance_matrix2=np.zeros((number_taxa-2,number_taxa-2))
        for i,k in enumerate (distance_matrix):
            if len(k)==0:
                break
            for j in range(0,number_taxa-2):
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
  
    
    
    def delta_value_calculator(self):
        delta_list=[]
        
        original_matrix = self.distance_matrix()
        
        s = list(range(0,original_matrix.shape[0]))
        for quartet in itertools.combinations(s,4):
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            delt_value = self.delta_value(quartet_matrix)
            delta_list.append(delt_value)
            
        print(f'Mean delta value is {(np.mean(np.array(delta_list)))}')
        return np.array(delta_list)
    
    def delta_plot1(self):
        '''
        

        Returns
        -------
        overall_delt : TYPE
            DESCRIPTION.

        '''
        
        delta_values_list = self.delta_value_calculator()
        print(f'Number of delta values are :{len(delta_values_list)}')
        
        overall_delt = plt.hist(delta_values_list,bins=10)
        #plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(10))
        plt.ylabel("Proportion of total numbers")
        plt.xlabel("Delta values")
        plt.show()
        return overall_delt
        
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
        
        s = list(range(0,original_matrix.shape[0]))

        for quartet in itertools.combinations(s,4):
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            if self.delta_value_2(quartet_matrix)[0] > 0.8:
                s = self.delta_value_2(quartet_matrix)[1]
                l = self.delta_value_2(quartet_matrix)[2]
                area=s*l
                if area<0.001:
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
    def individual_taxa_plot(self):
        
        mean_delta_list=[]
        
        original_matrix = self.distance_matrix()
        
        s = list(range(0,original_matrix.shape[0]))
        for const in s:
            delta_list=[]
            # filter out the constant number from the list
            nums_filtered = [n for n in s if n != const]
        
            for quartet in itertools.combinations(nums_filtered,3):
                
                quartet_list = list(quartet)
                quartet_list.append(const)
                quartet = tuple(quartet_list)
                quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
              
                quartet_matrix = np.array(quartet_matrix)
                delt_value = self.delta_value_2(quartet_matrix)[0]
                
                delta_list.append(delt_value)
            delta_list2 = np.array(delta_list)
            delt_mean = np.mean(delta_list2)
                
            mean_delta_list.append(delt_mean)
        print(mean_delta_list)
        mean_delta_list.sort(reverse=True)
        bar_taxa = plt.bar(s,mean_delta_list)
        plt.show()
            
                
                
        
        
        
        
sequence = delta_plot(nexus_file)

#list_mean_delta,array_delta_values = sequence.delta_random_mean()

#print(sequence.mean_delta_plot())

file_path = "C:/Users/willi/Documents/YEAR2/Labsheets programming 1/Week6/Delta_plot/Data/output.txt"

#sequence.second_try()
#print(sequence.delta_plot1())
#print(sequence.average_length_quartet())
print(sequence.individual_taxa_plot())
#print(sequence.average_length_quartet())
'''
matrix = np.array([
    [0, 0.124063, 0.74625, 0.095312],
    [0.124063, 0, 0.744489, 0.048432],
    [0.74625, 0.744489, 0, 0.746042],
    [0.095312, 0.048432, 0.746042, 0]
])
print(sequence.delta_value(matrix))
'''



