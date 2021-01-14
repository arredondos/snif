#!/usr/bin/python
import sys, argparse

parser = argparse.ArgumentParser(description="Converts the output of MS into psmcfa (the input file of psmc)")
parser.add_argument("-s", "--bin_size", type=int, default=100, help="The equivalent of bin_size in psmc")
parser.add_argument("input_ms_results", help="The file produced by MS")

args = parser.parse_args()

BinSize = args.bin_size
fname = args.input_ms_results

f = open(fname,'r')#, encoding='utf-16')
f = f.read()
f = f.split('\n')

SeqLength = int(f[0].split(' -r ')[1].split(' ')[1])

nBins = int(SeqLength / BinSize) + 1
    
count=0

for i in range(len(f)-1):
    
    templine = f[i].split(' ')
    
    if templine[0] == 'segsites:':
        
        segsites = int(templine[1])
    
    if templine[0] == 'positions:':
        
        count = count + 1
        A=['T'] * nBins
        
        for j in range(1,segsites+1):
        
            pos = int ( round ( float(SeqLength) * float(templine[j]) ) / BinSize )
            A[pos] = 'K'
        
        #print ">"+str(count)
        sys.stdout.write(">"+str(count)+"\n")
        
        for i in range(len(A)):
            if i>0 and i%60==0:
                sys.stdout.write('\n')
            sys.stdout.write(A[i])

        sys.stdout.write('\n')