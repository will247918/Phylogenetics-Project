# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 11:35:51 2023

@author: willi
"""


from Bio import AlignIO
import numpy as np
from matplotlib import pyplot as plt
import itertools
from matplotlib.ticker import PercentFormatter

#Where the nexus file is saved
nexus_file = "C:/Users/willi/BioinformaticsProject/Data/Insobium_outgroup.nex"


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
        #Add the numbers to a numpy array matrix
        distance_matrix2=np.zeros((63,63))
        for i,k in enumerate (distance_matrix):
            if len(k)==0:
                break
            for j in range(0,63):
                distance_value = k[j]
                distance_matrix2[i][j]=distance_value
        #Remove empty strings
        filtered_names = [name[0] for name in names_taxa if name]
        #Remove outer quotation marks
        names_taxa2 = [name[1:-1] for name in filtered_names]
        print(distance_matrix2)
        return distance_matrix2,names_taxa2


            
    def quartet_area(self,matrix):
        '''
        This function calculates the quartet area

        Parameters
        ----------
        matrix : numpy array
            Distance matrix of the quartet.

        Returns
        -------
        area : float
            The value of the area of the quartet

        '''
    
        distance_matrix = matrix
        #Create distance array
        da = list(distance_matrix[np.triu_indices(4)])
        
        # Distance for dxy,uv
        d1 = da[1]+da[8]
        #Distance for dxv,yu
        d2 = da[5]+da[3]
        #Distance for dxu,yv
        d3 = da[6]+da[2]
        distance_list = [d1,d2,d3]
        #Sort the distance to get the largest and the smallest values.
        distance_list.sort()
        
        #Multiply the different values of the edges of the quadrilateral in the 
        #quartet to output the area.
        area = ((distance_list[2]-distance_list[1])/2)*((distance_list[2]-distance_list[0])/2)
        
        return area      
    
    def a_value(self,matrix):
        '''
        This funciton calculates the a value in the quartet distance matrix.

        Parameters
        ----------
        matrix : TYPE
            DESCRIPTION.

        Returns
        -------
        a : float
            Calculated distance in the quartet diagram.

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
        
        a = (p[0][0]+p[1][0]-p[2][1])/2
        return a
    def b_value(self,matrix):
        '''
        Caluclates the distance value, b in the quartet.

        Parameters
        ----------
        matrix : numpy array
            Distance matrix of a quartet.

        Returns
        -------
        b : float
            Distance value b in the quartet.

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
        
        b = (p[1][0]+p[0][1]-p[2][0])/2 
        
        return b
    
    def c_value(self,matrix):
        '''
        Caluclates the distance value, c in the quartet.

        Parameters
        ----------
        matrix : numpy array
            Distance matrix of a quartet.

        Returns
        -------
        c : float
            Distance value c in the quartet.

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
        
        c = (p[1][1]+p[0][1]-p[2][1])/2
        return c
    def d_value(self,matrix):
        '''
        Caluclates the distance value, d in the quartet.

        Parameters
        ----------
        matrix : numpy array
            Distance matrix of a quartet.

        Returns
        -------
        d : float
            Distance value d in the quartet.

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
        
        d = (p[0][0]+p[1][1]-p[2][0])/2
        return d
        
    def s_l_value(self,matrix):
        '''
        This funciton calculates the value of s and l in the quartet diagram.

        Parameters
        ----------
        matrix : numpy array 
            This contains the distance matrix for the quartet.

        Returns
        -------
        s : float
            Distance value representing the width of the quadrilateral
        l : float
            Distance value representing the length of the quadrilateral.

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
     
        x= abs((p[2][0]+p[2][1]-p[0][0]-p[0][1])/2)
        
        y= abs((p[2][0]+p[2][1]-p[1][0]-p[1][1]) /2)
        
        s = min(x,y)
        
        l = max(x,y)
        
        return s,l
        
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
        list_1 = [distance_array[1],distance_array[8],distance_array[5],distance_array[3],distance_array[6],distance_array[2]]
        list_1.sort()
        x = list_1[4]+list_1[5]-list_1[0]-list_1[1]
       
        y= list_1[4]+list_1[5]-list_1[2]-list_1[3]
        
        s = min(x,y)
        
        l = max(x,y)
        
        
        distance_list = [d1,d2,d3]
        distance_list.sort()
        
        #print(distance_list[2]-distance_list[1])
        #print(distance_list[2]-distance_list[0])
        delta_value = (distance_list[2]-distance_list[1])/(distance_list[2]-distance_list[0])
        #print(f'Delta value is {distance_list}')
        
        return delta_value
    
    def delta_value_2(self,matrix):
        '''
        This is an alternative way to calculate the delta value, which uses the 
        distances in the delta plot.

        Parameters
        ----------
        matrix : numpy array
            This contains the distance matrix for a quartet of 4 taxa.

        Returns
        -------
        float
            Returns the delta value,for the quartet matrix inputted.

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
        
        #Create a list of lists of the different distance values and sort the
        #list by summing the individual values in the list
        list_of_lists = [l_1,l_2,l_3]
        p = sorted(list_of_lists, key=lambda x: sum(x))
        print(p)
        #Extract the different values of the list to calculate the distances
        #Of the values on the edge of the quadrilateral in the quartet.
        x= abs((p[2][0]+p[2][1]-p[0][0]-p[0][1])/2)
        
        y= abs((p[2][0]+p[2][1]-p[1][0]-p[1][1]) /2)
        
        #Calculate the minimum and maximum of the edges and divide to output
        #delta value
        s = min(x,y)
        
        l = max(x,y)
        if s ==0 or l==0:
            delta=0
        else:
            delta=s/l
        
        return delta
        
    
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
        for i,j in enumerate(array2):
            for k in range(0,4):
                different_nucleotides=0
                #Compare the sequence alignments of the different nucleotides
                #If the nucleotides are not the same add one to the different
                #nucleotides. At the end then append this value to the distance
                #matrix
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
        
    
    def list_area(self):
        '''
        This funciton calculates the area of the quardilateral in the quartet
        diagram and appends it to a list

        Returns
        -------
        a_list : list
            Area of each quartet.

        '''
        #Define lists
        area_list=[]
        a_list = []
        original_matrix = self.distance_matrix()[0]
        #Define a dictionary containing indexes of taxa and iterate through all
        #possible quartet combinations of them 
        s = list(range(0,original_matrix.shape[0]))
        for quartet in itertools.combinations(s,4):
            #Form a distance matrix from the quartet indexes
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            
            #Calculate the values of s and l and multiply to calculate the area
            #of the quadrilateral
            s,l = self.s_l_value(quartet_matrix)
            a_list.append(s*l)
       
        return a_list
    
    def perimeter_area_method(self):
        '''
        Method to compare the ratio of the perimeter to the area of every quartet,
        to form a new rescalled delta plot.

        Returns
        -------
        area_perimeter_list : list
            This is a list containing the area/perimeter values for all quartets

        '''
        #Initialise empty lists
        area_list=[]
        perimeter_list=[]
        area_perimeter_list=[]
        #Import distance matrix
        original_matrix = self.distance_matrix()[0]
        #Define a dictionary of all possible taxa values
        s = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
        #Iterate through all quartet combinations, using itertoools
        for quartet in itertools.combinations(s,4):
            #From the index of the quartet, produce a quartet_matrix, showing 
            #the distances.
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
           
            
            #Define the value of s and l of the quartet diagram, which represent
            #the sides of the quadrilateral
            s,l = self.s_l_value(quartet_matrix)
            #Calculate the area of the quadrilateral and append to the list
            area_list.append(s*l)
            
            #Calculate the perimeter and append to a list
            perimeter = 2*s+s*l
            perimeter_list.append(perimeter)
            
            #Calculate the ratio of area/perimeter for each quartet
            area_perimeter_list.append((s*l)/(perimeter))
            
        
        return area_perimeter_list
    
    def min_diagonal_perimeter(self):
        '''
        This function caluculates the ratio of the minimum diagonal distance 
        value in the quartet diagram to the perimeter of the quadrilateral.

        Returns
        -------
        list_delta : list
            List contianing the ratio of the minimum diagonal/perimeter.

        '''
        #Initialise final list
        list_delta=[]
        #Obtain the overall matrix
        original_matrix = self.distance_matrix()[0]
        #Define dictionary and use itertools to go through all index quartets
        s = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
        for quartet in itertools.combinations(s,4):
            #Calculate the distance matrix for each quartet
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            #Calculate the difefrent distance values in the quartet
            a,b,c,d = self.a_value(quartet_matrix),self.b_value(quartet_matrix),self.c_value(quartet_matrix),self.d_value(quartet_matrix)
            s,l = self.s_l_value(quartet_matrix)
            #Calculate the rescalled delta value by calulctaing the width +
            #teh length of the quadrilateral divided by the diagonal length 
            re_delta=min(((l+s)/(a+s+l+c)),((l+s)/(d+s+l+b)))
            list_delta.append(re_delta)
        
        return list_delta
    def average_diagonal_perimeter(self):
        '''
        Calculates the mean value of the sum of the length and width of the 
        quadrilateral divided by the different diagonal lengths.

        Returns
        -------
        list_delta : list
            List of the rescalled delta values for this example.

        '''
        list_delta=[]
        original_matrix = self.distance_matrix()[0]
        #Iterate through different combinations of the quartet indexes
        s = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
        for quartet in itertools.combinations(s,4):
            print(quartet)
            #Form the distance quartet matrix
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            print(quartet_matrix)
            #calculate the distance values
            a,b,c,d = self.a_value(quartet_matrix),self.b_value(quartet_matrix),self.c_value(quartet_matrix),self.d_value(quartet_matrix)
            s,l = self.s_l_value(quartet_matrix)
            array_values = np.array((l+s)/(a+s+l+c),(l+s)/(d+s+l+b))
            #Append the rescalled delta value to a list.
            re_delta=np.mean(array_values)
            print(re_delta)
            list_delta.append(re_delta)
        
        return list_delta
    
    def area_length_mean(self):
        '''
        Method that rescalles the delta value by using the 
        value of the area of the quartet divided by the 
        length between diagonals
        

        Returns
        -------
        list_delta : list
            List of all rescalled delta value using the method explained.

        '''
        list_delta=[]
     
        original_matrix = self.distance_matrix()[0]
        #Define a dictionary containing indexes of taxa and iterate through all
        #possible quartet combinations of them 
        s = list(range(0,original_matrix.shape[0]))
        for quartet in itertools.combinations(s,4):
            #Form a distance matrix from the quartet indexes
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            
            quartet_matrix = np.array(quartet_matrix)
            #Calculate value of a,b,c and d
            a,b,c,d = self.a_value(quartet_matrix),self.b_value(quartet_matrix),self.c_value(quartet_matrix),self.d_value(quartet_matrix)
            
            
            #Calculate the values of s and l and multiply to calculate the area
            #of the quadrilateral
            s,l = self.s_l_value(quartet_matrix)
            #Method for dividing my minimum diagonal
            #one_diagonal = a+s+l+c
            #other_diagonal = d+s+l+b
            #array_value = s*l/min(one_diagonal,other_diagonal)
            
            #Method for dividing my mean diagonal
            array_value = (s*l)/((a+s+l+c+d+s+l+b)/2)
            
            #Just the area
            #array_value = s*l
            list_delta.append(array_value)
            
            #Just the area
            
        max_delta = max(list_delta)
        min_delta = min(list_delta)
        norm_data=[]
        for i in list_delta:
            norm_data.append((i-min_delta)/(max_delta-min_delta))
            
        
        return norm_data,min_delta,max_delta
            
        '''
        for quartet in itertools.combinations(s,4):
            q_num +=q_num
            #Form the distance quartet matrix
            quartet_matrix = [[original_matrix[i][j] for j in quartet] for i in quartet]
            quartet_matrix = np.array(quartet_matrix)
            #calculate the distance values
            a,b,c,d = self.a_value(quartet_matrix),self.b_value(quartet_matrix),self.c_value(quartet_matrix),self.d_value(quartet_matrix)
            s,l = self.s_l_value(quartet_matrix)
            for i in self.list_area():
                array_value = i/((a+s+l+c+d+s+l+b)/2)
                list_delta.append(array_value)
            #Append the rescalled delta value to a list.
            print(f'q_num is {q_num}')
            return list_delta
    
    '''
    def delta_value_3(self,matrix):
        distance_matrix = matrix
        #Create distance array
        da = list(distance_matrix[np.triu_indices(4)])
        

        a,b,c,d = self.a_value(distance_matrix),self.b_value(distance_matrix),self.c_value(distance_matrix),self.d_value(distance_matrix)
       
       
        #Calculate the values of s and l and multiply to calculate the area
        #of the quadrilateral
        s,l = self.s_l_value(distance_matrix)
        #For area divided by minimum length
        #delta_value = (s*l)/min(a+s+l+c,d+s+l+b)
        
        #Area divided by average length
        #delta_value = (s*l)/(a+s+l+c+d+s+l+b)/2
        
        #Just the area
        delta_value = s*l
        
        
        
       
        
        return delta_value
    
    def delta_plot1(self):
        '''
        

        Returns
        -------
        overall_delt : TYPE
            DESCRIPTION.

        '''
        list_cutoff=[]
        #list_area = self.list_area()
        delta_values_list = self.area_length_mean()[0]
        
        
        
        
        for j in range(0,len(delta_values_list)):
            if delta_values_list[j]>0.8:
               list_cutoff.append(j) 
        print(f'Mean delta value is {np.mean(np.array(delta_values_list))}')
        overall_delt = plt.hist(delta_values_list,bins=10,density=True)
        plt.gca().yaxis.set_major_formatter(PercentFormatter(10))
        plt.ylim(0,7)
        
        plt.show()
        print(f'Number of quartets: {len(list_cutoff)}')
        return overall_delt
    
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
        #Extract the names of the taxa inputted from the NEXUS file.
        names_taxa = self.distance_matrix()[1]
       
        #List defining all the taxa (from 1 to n)
        s = list(range(0,original_matrix.shape[0]))
        min_delta = self.area_length_mean()[1]
        max_delta = self.area_length_mean()[2]
        
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
                delt_value = self.delta_value_3(quartet_matrix)
                delta_value_norm = (delt_value-min_delta)/(max_delta-min_delta)
                #Append to a list, representing the delta value list, where 
                #each quartet contains the same taxon.
                delta_list.append(delta_value_norm)
            
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
        bar_taxa = plt.bar(s,sorted_delta)
        plt.show()
        return bar_taxa,sorted_names
    
    def cut_off(self):
        '''
        This method incorporates the 'cut-off' method. If the delta value is 
        above a certain value, these quartets are isolated. Then if the area 
        of the quartet is also below a certain value the number of star like 
        quartets can be identified.

        Returns
        -------
        overall_delt : TYPE
            DESCRIPTION.

        '''
        list_cutoff=[]
        #list_area = self.list_area()
        delta_values_list = self.area_length_mean()[0]
      
        
        for j in range(0,len(delta_values_list)):
            if delta_values_list[j]>=0.8:
               list_cutoff.append(j) 
        print(f'Number of quartets: {len(list_cutoff)}')
      
        
    
sequence = delta_plot(nexus_file)

#list_mean_delta,array_delta_values = sequence.delta_random_mean()

#print(sequence.mean_delta_plot())

file_path = "C:/Users/willi/Documents/YEAR2/Labsheets programming 1/Week6/Delta_plot/Data/output.txt"

'''
# (19,21,22,23)
matrix = np.array([[0,0.745129,0.02357,0.731301],
                  [0.745129,0,0.739786,0.73350],
                  [0.02357,0.739786,0,0.732244],
                  [0.731301,0.733501,0.73224,0]])



print(sequence.delta_value(matrix))

a,b,c,d = sequence.a_value(matrix),sequence.b_value(matrix),sequence.c_value(matrix),sequence.d_value(matrix)
print(f'a is {a}, b is {b}, c is {c},d is {d}')
s,l = sequence.s_l_value(matrix)
print(f' s is {s}, l is {l}')
array_values = np.array((l+s)/(a+s+l+c),(l+s)/(d+s+l+b))
#Append the rescalled delta value to a list.
re_delta=np.mean(array_values)
print(sequence.delta_value(matrix))
print(re_delta)
'''
#print(sequence.delta_plot1())
#print(sequence.cut_off())
#print(sequence.individual_taxa_plot())
x = np.array([[ 0, 0.022888, 0.023449, 0.023112],
 [0.022888, 0, 0.005722, 0.005385],
 [0.023449, 0.005722, 0, 0.003478],
 [0.023112, 0.005385, 0.003478, 0]])
sequence.delta_value_2(x)



