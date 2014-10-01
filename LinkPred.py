# Traverse a bipartite graph and compute the aa similarity of each node pair
# Usage: python aa_tranverse.py [input file] [output file]
# The MIT License (MIT)
# dataset used : edge list file inputdata.txt (dataset from drugbank of 727 drugs)
# Copyright (c) 2014 <Abhik Seal>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import sys
import math

# Check command line arguments
args = sys.argv
if len(args) != 3:
	print "Usage: python " + args[0] + " <input graph> <output file>"
	quit()
	
inputGraph = str(args[1])
outputFile = str(args[2])

leftNodes = {}
rightNodes = {}

# Dict for the weight of each pair of nodes (increment this)
weights = {}

# Read pairs of files into a dict
#pairs = []
input = open(inputGraph, "r")
for line in input:
	line = line.rstrip()
	[left, right] = line.split("\t")
	
	# Add to dictionaries
	if left in leftNodes:
		leftNodes[left].add(right)
	else:
		leftNodes[left] = set([right])
	
	if right in rightNodes:
		rightNodes[right].add(left)
	else:
		rightNodes[right] = set([left])

# Traverse the dictionaries for every path A, n1, n2, B for every 
# A and B not connected already
loops = len(leftNodes.keys())
count = 0
for startNode in leftNodes.keys():
	percent = (float(count) / loops) * 100
	print str(percent) + "%"
	count += 1
	
	# For all neighbors of start node (get n1)
	for n1 in leftNodes[startNode]:
		# For all neighbors of n1 (get n2)
		for n2 in rightNodes[n1]:
			# For all neighbors of n2 not neighbors
			# neighbors with startNode (get endNode)
			for endNode in leftNodes[n2]:
				if endNode in leftNodes[startNode]:
					continue
				# Calculate the running similarity between these nodes
				# Amount to increment = 1 / ln(deg(n1) + deg(n2))
				inc = 1 / math.log(len(rightNodes[n1]) + len(leftNodes[n2]))
				
				# Add the entry to the weights dictionary, making
				# sure it exists in the dictionary first
				if startNode in weights:
					if endNode in weights[startNode]:
						weights[startNode][endNode] += inc
					else:
						weights[startNode][endNode] = inc
				else:
					weights[startNode] = {endNode: inc}
					
# print weights to a file
# print weight first for easier sorting
output = open(outputFile, "wb")
for i in weights.keys():
	for j in weights[i].keys():
		output.write(str(weights[i][j]) + "\t" + str(i) + "\t" + str(j) + "\n")