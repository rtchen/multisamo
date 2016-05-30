import numpy as np
import argparse 
import operator
import itertools

objective = 0

def read_file(filename):
    fp = open(filename)
    line = fp.readline()
    numProteins = int(line.split()[1])
    line = fp.readline()
    numPos = int(line.split()[1])
    numAminos = numPos * numProteins
    line = fp.readline()
    aminoLen = [int(s) for s in line.split(',')]
    matches = np.loadtxt(filename,skiprows = 3,delimiter = ',')
    fp.close()
    return numProteins,numPos,aminoLen,matches

def map_protein(numProteins, numPos):
    numAminos = numProteins * numPos
    amino2protein = [0] * numAminos
    for i in range(numAminos):
        amino2protein[i] = i/numPos
    return amino2protein


def find_clique(matches,amino2protein,numProteins,numPos):
    numAminos = numProteins*numPos
    clique_map = {}
    for i in range(numAminos):
        for j in range(numAminos):
            if matches[i][j] == 0:
               continue
            for k in range(numAminos):
                if matches[i][k] ==1 and matches[j][k]==1:
                   key = tuple(sorted((amino2protein[i],amino2protein[j],amino2protein[k])))
                   if key not in clique_map:
                      clique_map[key] =1
                   else:
                      clique_map[key] = clique_map[key]+1
    return clique_map


def align_seeds(classification, matches, seeds,class_sets,numPos):
   global objective
   proteinA = seeds[0]
   for i in range(numPos):
      class_sets[i].append(proteinA*numPos+i)
      classification[proteinA][i] = i
   proteinB = seeds[1]
   unclassed = []
   for i in range(numPos):
       flag = 0
       for j in range(numPos):
           if matches[proteinB*numPos+i][proteinA*numPos+j] == 1:
              class_sets[j].append(proteinB*numPos+i)
              classification[proteinB][i] = j
              objective = objective + 1
              flag = 1
       if flag == 0:
           unclassed.append(i)
   for i in range(len(unclassed)):
        for j in range(numPos):
            if(len(class_sets[j])<2):
               class_sets[j].append(proteinB*i+unclassed[i])
               classification[proteinB][unclassed[i]] = j
               break

def add_layer(classification,matches,class_sets,numPos,protein,priority):
	global objective
	edgesNum = np.zeros((numPos,numPos))
	for i in range(numPos):
	    for j in range(numPos):
	        for k in range(len(class_sets[j])):
	            if matches[protein*numPos+i][class_sets[j][k]] == 1:
	               edgesNum[i][j] = edgesNum[i][j]+1
	maxEdgeSum = 0
	perms = list(itertools.permutations(range(0,numPos)))
	for i in range(len(perms)):
	    edgeSum = 0
	    for j in range(numPos):
	       classNo = perms[i][j]
	       edgeSum = edgeSum + edgesNum[j][classNo]
	    if edgeSum > maxEdgeSum:
	       maxEdgeSum = edgeSum
	       maxPerm = perms[i]
	objective = objective + maxEdgeSum
	for i in range(numPos):
	    class_sets[maxPerm[i]].append(protein*numPos+i)
	    classification[protein][i] = maxPerm[i]
	priority[protein] = -1
  

def initialize_priority(clique_map,aligned,numProteins,priority):
    proteinA = aligned[0]
    proteinB = aligned[1]
    for i in range(numProteins):
        if priority[i] == -1:
           continue
        key = tuple(sorted((proteinA,proteinB,i)))
        if key in clique_map:
           priority[i] = priority[i] + clique_map[key]


def update_priority(clique_map,aligned,numProteins,priority,proteinC):
    for proteinA in range(numProteins):
        if priority[proteinA] == -1:
           continue
        for proteinB in aligned:
           key = tuple(sorted((proteinA,proteinB,proteinC)))
           if key in clique_map:
              priority[proteinA] = priority[proteinA] + clique_map[key]    
    aligned.append(proteinC)

def multiple_align(clique_map,aligned,priority,numProteins,classification,matches,class_sets,numPos):
    while(len(aligned)!=numProteins):
        protein = priority.index(max(priority))
        add_layer(classification,matches,class_sets,numPos,protein,priority)
        update_priority(clique_map,aligned,numProteins,priority,protein)
        


if __name__== "__main__" :
   parser = argparse.ArgumentParser(description='protein multiple alignment')
   parser.add_argument('-i', '--dataFile',help='data file')
   parser.add_argument('-o','--outputFile',help='output alignment file')
   args = parser.parse_args()
   dataFile = args.dataFile
   numProteins,numPos,aminoLen,matches = read_file(dataFile)
   amino2protein = map_protein(numProteins, numPos)
   clique_map = find_clique(matches,amino2protein,numProteins,numPos)
   seeds = max(clique_map.iteritems(), key=operator.itemgetter(1))[0]
   classification = np.zeros((numProteins,numPos))
   class_sets = [[] for i in range(numPos)]
   priority = [0] * numProteins
   priority[seeds[0]] = -1
   priority[seeds[1]] = -1
   aligned = []
   align_seeds(classification,matches,seeds,class_sets,numPos)
   aligned.append(seeds[0])
   aligned.append(seeds[1])
   initialize_priority(clique_map,aligned,numProteins,priority)
   add_layer(classification,matches,class_sets,numPos,seeds[2],priority)
   update_priority(clique_map,aligned,numProteins,priority,seeds[2])
   multiple_align(clique_map,aligned,priority,numProteins,classification,matches,class_sets,numPos)
   print objective *2 # notice that in integer programming the objective of every connected edge is counted twice, so Im 
                      # multiplying 2 here just for comparison
   print classification












