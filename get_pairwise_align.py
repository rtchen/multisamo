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
   parser.add_argument('-o','--outputFile',help='output distance file, suffix has to be dzn ')
   args = parser.parse_args()
   clusterFile = args.clusterFile
   outputFile = args.outputFile
   pairs = [pairs.rstrip('\n') for pairs in open(clusterFile)]
   numProteins = len(pairs)
   if numProteins < 2:
      sys.exit("not enough proteins")
   for i in range(numProteins):
      for j in range(i+1,numProteins):
         command = "samo0.exe -pocket %s %s -os temp.txt > screen.txt"%(pairs[i],pairs[j])
         os.system(command)
         fp = open("temp.txt")
         matches = parse_result()
         if i==0 and j==1:
            numAminos= len(matches)
            distance = np.zeros((numAminos*numProteins,numAminos*numProteins),dtype = int)
         for k in range(numAminos):
            if matches[k]!=-1:
               aminoj = i*numAminos+k
               aminok = j*numAminos+matches[k]
               distance[aminoj][aminok] = 1
               distance[aminok][aminoj] = 1
   fp.close()
   os.system("del temp.txt")
   os.system("del screen.txt")
   fp = open(outputFile,'w')
   fp.write("int: npos = %d;\n" % (numAminos))
   fp.write("int: nProteins= %d;\n" % (numProteins))
   fp.write("array[1..nAminos,1..nAminos] of int: adj=\n")
   fp.write("[")
   for i in range(numAminos*numProteins):
    fp.write("|")
    for j in range(numAminos*numProteins-1):
      fp.write("%d ," % distance[i][j])
    fp.write("%d\n" % distance[i][numAminos*numProteins-1])
   fp.write("|];\n")
   fp.close()



