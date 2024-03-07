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
import math

#Where the nexus file is saved
nexus_file = "C:/Users/willi/BioinformaticsProject/Data/Tomato_wilt_remove_recombinant.nex"


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
        #Initialise empty lists
        distance_matrix =[]
        names_taxa=[]
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
                    
                    names_taxa.append(line[1:2])
                    distance_matrix.append(line[2:length_line])
        
        file.close()
        #print(distance_matrix)
        #Add the numbers to a numpy array matrix
        distance_matrix2=np.zeros((21,21))
        for i,k in enumerate (distance_matrix):
            if len(k)==0:
                break
            for j in range(0,21):
                distance_value = k[j]
                distance_matrix2[i][j]=distance_value
        #Remove empty strings
        filtered_names = [name[0] for name in names_taxa if name]
        #Remove outer quotation marks
        names_taxa2 = [name[1:-1] for name in filtered_names]
        
        return distance_matrix2,names_taxa2

            
            
    
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
        
        if s ==0 or l==0:
            delta=0
        else:
            delta=s/l
        
        return delta,s,l
  
    
    
    def delta_value_calculator(self):
        '''
        This caluclates the delta values for each possible quartet in the 
        distance matrix.

        Returns
        -------
        array
            This array contains all the delta values from all possible quartets
            in the distance matrix.

        '''
        delta_list=[]
        
        #Extract the overall distance matrix
        original_matrix = self.distance_matrix()[0]
        
        #List of all taxa from 1 to n
        s = list(range(0,original_matrix.shape[0]))
        '''
        Use itertools to find every combination of 4 taxa (quartet). For each 
        of these quartets calculate the delta value. Then append the delta 
        values to a list.
        '''
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
        
        overall_delt = plt.hist(delta_values_list,bins=10,ec='White',density=True)
        plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(10))
        plt.ylim(0,9)
        plt.ylabel("Proportion of total number of quartets")
        plt.xlabel("Delta values")
        plt.title("Delta plot for Tomato spotted wilt virus taxa without \n recombinants ")
        plt.show()
        return overall_delt
        
    def cut_off_method(self):
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
        #Initiate 3 empty lists to store length values
        l1=[]
        l2=[]
        #Area list for smaller area
        area_list_small = []
        #Area list for normal area
        area_list=[]
        #Call the function to get the original distance matrix
        original_matrix = self.distance_matrix()[0]
        #Keep track of the number of quartets
        quartet_number =0 
        s = list(range(0,original_matrix.shape[0]))

        for quartet in itertools.combinations(s,4):
            quartet_number+=1
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            if self.delta_value_2(quartet_matrix)[0] >=0.8:
                s = self.delta_value_2(quartet_matrix)[1]
                l = self.delta_value_2(quartet_matrix)[2]
                area=s*l
                area_list.append(area)
                if area<=4.88978*10**-6:
                    
                    #append the length values to the empty lists
                    l1.append(s)
                    l2.append(l)
                    area_list_small.append(area)
      
        
        l1 = np.array(l1)
        l2=np.array(l2)
        area_list = np.array(area_list)
        proportion = len(area_list_small)/quartet_number
        return [l1.mean(),l2.mean(),area_list.mean(),proportion*100]
    def individual_taxa_plot(self):
        '''
        This function aims to plot the mean delta values for individual taxa in
        decending order. This works by calculating the mean delta value for 
        each taxa and then sorting by size. This is extended to sorting the 
        names as well.

        Returns
        -------
        bar_taxa : A bar chart showing the mean delta value for each taxa.
        sorted_names : A list of names of taxa, which correspond to the bars
        in the plot.
        

        '''
        #Empty list to store mean delta values
        mean_delta_list=[]
        
        #Overall distance matrix from the data
        original_matrix = self.distance_matrix()[0]
        print(original_matrix)
        #Extract the names of the taxa inputted from the NEXUS file.
        names_taxa = self.distance_matrix()[1]
       
        #List defining all the taxa (from 1 to n)
        
        s = list(range(0,original_matrix.shape[0]))

        '''
        This is a for loop to remove a taxa from the list s. Then calculate
        all possible combinations of 3 taxa from the remaining list, whilst 
        appending the taxa number back on. Therefore, for each taxa this will 
        calculate all possible quartets. Then the delta value is calculated 
        for all quartets containing a given taxa, hence the mean delta value 
        for each taxa can be calculated.
        '''
        for const in s:

            delta_list=[]
            # filter out the constant number from the list
            nums_filtered = [n for n in s if n != const]

            #Use itertools to iterate through all possible triplets
            for quartet in itertools.combinations(nums_filtered,3):
                
                quartet_list = list(quartet)
                #Append back the taxa number removed.
                quartet_list.append(const)
                quartet = tuple(quartet_list)
                #Extract the quartet from the overall distance matrix
                quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
                
                quartet_matrix = np.array(quartet_matrix)
                #Calculate the delta value
                delt_value = self.delta_value_2(quartet_matrix)[0]
                #Append to a list, representing the delta value list, where 
                #each quartet contains the same taxon.
                delta_list.append(delt_value)
                
            
            delta_list2 = np.array(delta_list)
            
            #Calculate the mean of delta value from the list for each taxon
            delt_mean = np.mean(delta_list2)
            
            mean_delta_list.append(delt_mean)
        
        #Extract the names from the NEXUS file
        names_delta = list(zip(mean_delta_list,names_taxa))
        #Sort both the mean delta value for each taxon and their corresponding 
        #names
        sort_names_delta = sorted(names_delta,key = lambda x:x[0],reverse=True)
        #Convert back to individual sorted lists 
        sorted_delta,sorted_names = zip(*sort_names_delta)
        #Plot the bar chart for each taxa
        bar_taxa = plt.bar(np.array(s)+1,sorted_delta)
        plt.ylim(0,0.5)
        plt.xlabel("Taxa")
        plt.ylabel("Mean delta value")
        plt.xticks(np.array(s)+1)
        plt.title("Mean delta values of the large RNA of Tomato spotted \n wilt virus without recombinant taxa")
        '''
        plt.annotate("Recombinant taxa",xy=(0.7,0.3),xytext=(0.5,0.45), arrowprops=dict(facecolor='red',  arrowstyle='->', linewidth=1),
             fontsize=12, color='blue')
        plt.annotate("Recombinant taxa",xy=(4,0.18),xytext=(0.5,0.45), arrowprops=dict(facecolor='red',  arrowstyle='->', linewidth=1),
             fontsize=12, color='blue')
        '''
        plt.show()
        print(sorted_delta[0]-sorted_delta[1])
        
        return bar_taxa,sorted_names
                
        
        
        
        
sequence = delta_plot(nexus_file)

#list_mean_delta,array_delta_values = sequence.delta_random_mean()

#print(sequence.mean_delta_plot())

file_path = "C:/Users/willi/Documents/YEAR2/Labsheets programming 1/Week6/Delta_plot/Data/output.txt"

#sequence.second_try()
#print(sequence.distance_matrix())
#print(sequence.average_length_quartet())
#print(sequence.individual_taxa_plot())
#print(sequence.cut_off_method())
print(sequence.delta_plot1())
'''
matrix = np.array([
    [0, 0.124063, 0.74625, 0.095312],
    [0.124063, 0, 0.744489, 0.048432],
    [0.74625, 0.744489, 0, 0.746042],
    [0.095312, 0.048432, 0.746042, 0]
])
print(sequence.delta_value(matrix))
'''



