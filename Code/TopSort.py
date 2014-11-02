#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  TopSort.py
#  
#  Copyright 2013 Raul Reina <raul@RMNet>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#
def make_link(G, node1, node2):
    if node1 not in G:
        G[node1] = {}
    (G[node1])[node2] = 1
    return G

def make_graph(connections):
    G = {}
    for (x,y) in connections:
        make_link(G,x,y)
    return G

def mark_component(G, node, marked, degrees):
    marked[node]=True
    for neighbor in G.get(node, []):
        if neighbor not in marked:
            degrees[neighbor]=1
            mark_component(G, neighbor, marked, degrees)
        else:
            degrees[neighbor]=degrees.get(neighbor, 0)+1

def get_indegree(G):
    degrees={}
    marked={}
    for node in G.keys():
        degrees[node]=0

    for node in G.keys():
        if node not in marked:
            mark_component(G, node, marked, degrees)
    return degrees

def topological_ordering(G):
    degrees=get_indegree(G)
    result=[]
    queue=[]
    for de in degrees:
        if degrees[de]==0:
            queue.append(de)

    while len(queue)>0:
        node=queue.pop(0)
        result.append(node)
        for neighbor in G.get(node, []):
            degrees[neighbor]=degrees[neighbor]-1
            if degrees[neighbor]==0:
                queue.append(neighbor)

    return result

def topSort(priorities):
	G = make_graph(priorities)
	return topological_ordering(G)
	
if __name__=='__main__':
    connections = [("ORD", "SEA"), ("ORD", "LAX"), ('ORD', 'DFW'), ('ORD', 'PIT'),
           ('SEA', 'LAX'), ('LAX', 'DFW'), ('ATL', 'PIT'), ('ATL', 'RDU'),
           ('RDU', 'PHL'), ('PIT', 'PHL'), ('PHL', 'PVD')]
    
    rules = xrange(4)
    connections = [(0,1), (0,3), (2,3)]
    G=make_graph(connections)
    TS = topological_ordering(G)
    TS.extend(list(set(rules) - set(TS)))
    print TS
