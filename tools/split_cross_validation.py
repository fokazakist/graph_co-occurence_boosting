#!/usr/bin/env python

import glob
import sys 


if len(sys.argv) != 3:
    print 'Usage: # python %s [split filename] [make file head]' % sys.argv[0]
    print 'Example: # python %s cpdb.gsp split/cpdb' % sys.argv[0]
    quit() 
filename = sys.argv[1]
head = sys.argv[2]

pos=0
neg=0

for i in range(10):
    pos=0
    neg=0
    f = open(filename)
    line = f.readline()

    trainname = head + "train" + str(i) + ".gsp"
    testname  = head + "test"  + str(i) + ".gsp"
    trf = open(trainname,'w')
    tef = open(testname,'w')
    
    while line:
        if line[:1]=="t":
            sline=line.split()
            print sline
            if sline[3]=='1':
                if pos % 10 == i:
                    tef.write(line)
                else:
                    trf.write(line)
                line=f.readline() 
                while line[:1]!="t" and line:
                    if pos % 10 == i:
                        tef.write(line)
                    else:
                        trf.write(line)
                    line=f.readline()
                pos = pos + 1
            else:
                if neg % 10 == i:
                    tef.write(line)
                else:
                    trf.write(line)
                line=f.readline()
                while line[:1]!="t" and line:
                    if neg % 10 == i:
                        tef.write(line)
                    else:
                        trf.write(line)
                    line=f.readline()
                neg = neg + 1
    f.close()
    trf.close()
    tef.close()
    
print "pos :",pos
print "neg :",neg


