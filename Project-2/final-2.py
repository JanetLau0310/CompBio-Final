#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os

# read edge list to networkx 
# the format of each line: (protein-1, protein-2, confidence of the interaction existence )
def createGraph(filename) :
    G = nx.Graph()
    for line in open(filename) :
        strlist = line.split('\t')
        n1 = strlist[0]
        n2 = strlist[1]
        weight = float(strlist[2])
        G.add_weighted_edges_from([(n1, n2, weight)]) #G.add_edges_from([(n1, n2)])
    return G

yeast_file = 'yeast.ppi'
G = createGraph(yeast_file)


# ### (b)Compute the degree distribution for the 5,001 nodes of the graph
# plot it as a bar graph with degree on the x-axis and counts on the y-axis, and submit this image as FPp2-degreedist.

import collections
degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
#print("Degree sequence", degree_sequence)
degreeCount = collections.Counter(degree_sequence)
deg, cnt = zip(*degreeCount.items())

fig, ax = plt.subplots()
plt.bar(deg, cnt, width=0.60)

#plt.hist(degree_values)
plt.xscale('log')
plt.title("Degree distribution")
plt.ylabel('Counts of the Nodes')
plt.xlabel("Degree")
plt.savefig('FPp2-degreedist.png')
plt.show()

# build the graph from the notes in class
# in order to test my result in the small size

G_test = nx.Graph()
G_test.add_nodes_from(['a','b','c','d','e','f','g'])
G_test.add_edges_from([('a','d'),('b','d'),('b','c'),('c','d'),('d','e'),('d','f'),('d','g'),('e','f'),('e','g')])

G1 = nx.Graph()
G1.add_nodes_from(['a','b','c','d','e','f','g'])
G1.add_edges_from([('a','b'),('a','e'),('b','c'),('b','e'),('b','f'),('c','f'),('c','e'),('d','f'),('e','f'),('g','f')])


def local_cluster(G, node):
    #print("degree ", G[node])
    x = sum([1 for u in G[node] for v in G[node] if not u == v and (u, v) in G.edges()])
    #print("x = ",x)
    if(len(G[node]) != 0 and len(G[node])!=1):
        y = (len(G[node]) * (len(G[node]) - 1))
        #print("y= ",y)
    else:
        return 0.0
    return  x/y

fp2c = 'FPp2c-clustcoef.txt'
nodes = list(G.nodes())

with open(fp2c,'a') as f:
    for node in nodes:
        f.write(node+'\t'+str(local_cluster(G,node))+'\n')

print("YGR296W: ",local_cluster(G,'YGR296W'))
print("YPL098C: " ,local_cluster(G,'YPL098C'))


# The result of local cluster coefficient is YGR296W = 0.1 and YPL098C = 0.8, so the latter is higher. <br>
# As the result shows, the difference between then is the Numerator, which means the number of triangles containing v. Consider that they have the same degree = 5, it means the interacting partners of YPL098C also interact with each other much more than YGR296W's partners.

def triangle(Graph):
    tri_test = list(nx.triangles(Graph).values())
    tmp = 0
    for i in tri_test:
        tmp += i
    return tmp/3.0

print(triangle(G))


# The number of triangles in G is 354514, I use the existing methods in networkX, which is the nx.triangles(), but this method will count the triangle of every nodes in Graph G, so to get the percise result, we need to divide the result by 3.

# ### (e)  A simple way to define the "closeness" is to use shortest path distance.
# In order to estimate the path length distribution for this graph, sample 1000 nodes at random from the graph, and compute the distribution of shortest path distances between these 1000 nodes.

from random import sample
random_nodes = sample(list(G.nodes()),1000)
new_edge = []
for i in range(0,1000):
    for j in range(i+1,1000):
        if G.get_edge_data(random_nodes[i],random_nodes[j]):
            new_edge.append((random_nodes[i],random_nodes[j]))

# build a new map, making it easier to deal with the following questions
G_tmp = nx.Graph()
G_tmp.add_nodes_from(random_nodes)
G_tmp.add_edges_from(new_edge)

tmp_nodes = list(G_tmp.nodes())
y = np.zeros((1000,1000),dtype=np.int)
for i in range(0,1000):
    for j in range(i+1,1000):
        if nx.has_path(G_tmp,tmp_nodes[i],tmp_nodes[j]):
            y[i][j] = nx.shortest_path_length(G_tmp, tmp_nodes[i],tmp_nodes[j])

def get_freq(arr):
    unique = np.unique(arr)
    res = [0]*len(unique)
    for i in range(len(arr)):
        for j in range(i+1,len(arr)):
            for k in range(len(unique)):
                if arr[i][j] == unique[k]:
                    res[k] +=1
    return res

xout = np.unique(y)
yout = get_freq(y)

'''
plt.title("Shortest Path Distribution")
plt.plot(xout,yout)
plt.xticks(xout)
plt.yscale('log')
plt.ylabel('Frequency')
plt.xlabel("Length of Shortest Path")
plt.savefig("FPp2e-spdists.png")
plt.show()
'''

tmp = list(nx.connected_components(G_tmp))
print(len(tmp))
diam = []
for c in nx.connected_components(G_tmp):
    diam.append(nx.diameter(G_tmp.subgraph(c)))
print(max(diam))
