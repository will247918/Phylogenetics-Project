# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 16:12:13 2024

@author: willi
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import collections
'''
#undirected graph
G=nx.Graph()

#Directed graph
G=nx.DiGraph()

#Multigraph undirected
G=nx.MultiGraph()
'''
'''
#Undirected Graph
G=nx.Graph()
#Define an edge between 1 and 2. Also implying there is
#node 1 and 2. Also define the weight, which can mean anything.
#Can also define 1 and 2 as strings.
G.add_edge(1,2,weight =0.9)

#Also can add a node
G.add_node("C")
G.add_edge(1,"C",weight=0.5)
G.add_edge(2,"C",weight =0.8)

nx.draw_spring(G,with_labels=True)
plt.show()
'''
'''
#Define edge list where 1 is connected to 2, 2 is connected to 3 etc..
edge_list = [(1,2),(2,3),(3,4),(3,5),(4,6),(6,7)]
G=nx.Graph()
G.add_edges_from(edge_list)

print(nx.adjacency_matrix(G))
'''



'''
If this represents A,B,C therefore A id connected to B
B is connected to C etc.. A one means there is a connection
'''
'''
G = nx.from_numpy_array(np.array([[0, 1, 0],
         [1, 0, 1],
         [0, 0, 0]]))

nx.circular_layout(G,with_labels=True)
plt.show()
'''
random_seed =1
def generate_weighted_graph(X):
    G=nx.Graph()
    G.add_nodes_from(X)
    for i in range(len(X)):
        for j in range(i+1,len(X)):
            weight = random.random()
            G.add_edge(X[i],X[j],weight=weight)
            
    return G
X = ['A','B','C','D','E','F','G']
G = generate_weighted_graph(X)
print(nx.shortest_path_length(G,X[0],X[1]))
nx.draw(G,with_labels=True)

def sort_values_by_multiple_keys_preserving_order(my_dict):
      """Sorts values in a dictionary by multiple keys in alphabetical order,
      preserving the original ordering for keys with the same value.
    
      Args:
        my_dict: A dictionary with values containing multiple keys.
    
      Returns:
        A new dictionary with sorted values.
      """
      sorted_dict = collections.OrderedDict()
      for key, value in my_dict.items():
            # Sort keys alphabetically with custom key function
            sorted_keys = sorted(value.keys(), key=lambda k: (value[k], k))
            # Create new OrderedDict with sorted keys and values
            sorted_value = collections.OrderedDict([(k, value[k]) for k in sorted_keys])
            sorted_dict[key] = sorted_value
      return sorted_dict

length = dict(nx.all_pairs_dijkstra_path_length(G))

distance_array = np.zeros((len(X),len(X)))
for i,j in enumerate(X):
    x = length.get(j)
    sorted_dict = dict(sorted(x.items()))
    print(sorted_dict)
    for k in range(0,len(X)):
        print(sorted_dict.get(X[k]))
        distance_array[i][k]=sorted_dict.get(X[k])
