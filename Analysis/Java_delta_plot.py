# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 19:02:51 2023

@author: willi
"""
from Bio import AlignIO
import numpy as np
from matplotlib import pyplot as plt
import itertools

class java_delta_plot():
     
    def __init__(self,file_path):
        
        self.file_path = file_path

    def delta_plot_java(self):
        delta_value_list=[]
        list_delta_file = open(self.file_path,'r')
        lines = list_delta_file.readlines()
        
        for line in lines:
            format(line.strip('\n'))
            delta_value_list.append(float(line.strip()))
        
        return delta_value_list
            
    def java_histogram(self):
        list_delta_values = self.delta_plot_java()
        
        overall_delt_plot = plt.hist(list_delta_values,bins=10)
        plt.show()
        return overall_delt_plot
        


file_path = "C:/Users/willi/Documents/YEAR2/Labsheets programming 1/Week6/Delta_plot/Data/output.txt"
sequence = java_delta_plot(file_path)
print(sequence.java_histogram())
   