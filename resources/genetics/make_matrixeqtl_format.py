#!/usr/bin/python

import os
import sys

#
def Chunk(L, n):
    for i in xrange(0, len(L), n):
        yield L[i:i+n]
#
        
#
if __name__ == '__main__':
    #
    tp  = sys.argv[1]
    
    with open(str(tp) + ".raw", "r") as H:
        
        #
        # HEADER
        #
        
        Head = [ x.split("_")[0] for x in H.next()[:-1].split(" ")[6:] ] # just SNPs
        
        i = 0
        for x in Chunk(Head, 10000):
            i += 1
            
            OUT = open(str(tp) + "." + str(i) + ".tab", "w")                
            OUT.write("\t".join( [ "" ] + x ) + "\n")
        
        #
        # DATA
        #
        
        for record in H:
            ALN  = record[:-1].split(" ")[0]
            data = record[:-1].split(" ")[6:]
            
            i = 0
            for x in Chunk(data, 10000):
                i += 1
                
                OUT = open(str(tp) + "." + str(i) + ".tab", "a")
                OUT.write("\t".join( [ ALN ] + x ) + "\n")
    
