import sys
import os
import numpy as np
import argparse
LINE_NUM = 7

def parse_args():
    args = parser.parse_args()

def parse_result():
    fp = open("temp.txt")
    for i, line in enumerate(fp):
       if i == LINE_NUM:
          matches = [int(s) for s in line.split()]
    fp.close()
    return matches

if __name__== "__main__" :
   parser = argparse.ArgumentParser(description='protein alignment in a cluster')
   parser.add_argument('-c', '--clusterFile',help='cluster file')
   parser.add_argument('-o','--outputFile',help='output distance file')
   args = parser.parse_args()
   clusterFile = args.clusterFile
   outputFile = args.outputFile
   pairs = [pairs.rstrip('\n') for pairs in open(clusterFile)]
   numProteins = len(pairs)
   numPairs = numProteins*numProteins
   aminoLen = [0]*numProteins
   matches = [None]*numPairs
   if numProteins < 2:
      sys.exit("not enough proteins")
   for i in range(numProteins):
      for j in range(i+1,numProteins):
         command = "samo0.exe -pocket %s %s -os temp.txt > screen.txt"%(pairs[i],pairs[j])
         os.system(command)
         matches[i*numProteins+j] = parse_result()
         if j==i+1:
            aminoLen[i]=len(matches[i*numProteins+j])
   command = "samo0.exe -pocket %s %s -os temp.txt > screen.txt"%(pairs[numProteins-1],pairs[numProteins-1])
   os.system(command)
   matches[numProteins*numProteins-1] = parse_result()
   aminoLen[numProteins-1] = len(matches[numProteins*numProteins-1])
   numPos = max(aminoLen)
   numAminos = numPos * numProteins;
   distance = np.zeros((numAminos,numAminos),dtype = int)
   for i in range(numProteins):
      for j in range(i+1,numProteins): 
        for k in range(len(matches[i*numProteins+j])):
          if matches[i*numProteins+j][k]!=-1:
             aminoj = i*numPos+k
             aminok = j*numPos+matches[i*numProteins+j][k]
             distance[aminoj][aminok] = 1
             distance[aminok][aminoj] = 1
   os.system("del temp.txt")
   os.system("del screen.txt")
   fp = open(outputFile,'w')
   fp.write("nummProteins %d\n" % (numProteins))
   fp.write("numPos %d\n" % (numPos))
   fp.write(",".join([str(n) for n in aminoLen]))
   for i in range(numAminos):
    fp.write("\n")
    for j in range(numAminos-1):
      fp.write("%d ," % distance[i][j])
    fp.write("%d\n" % distance[i][numAminos-1])
   fp.close()



