#!/usr/bin/python

import sys
import math
import os.path

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def printUsage():
	print "Usage: python qos.py <original file> <nn file>"
	exit(1)
pass;


if(len(sys.argv) != 3):
	printUsage()

origFilename 	= sys.argv[1]
nnFilename		= sys.argv[2]

nnFilename_new = nnFilename
for k in range (1,1000): 
	nnFilename_new 		= str(k) + nnFilename
	origLines 		= open(origFilename).readlines()
	if (os.path.isfile(nnFilename_new) != True): continue
	nnLines			= open(nnFilename_new).readlines()


	e = 0.0
	absError = 0.0

	for i in range(len(origLines)):
		origLine 	= origLines[i].rstrip()
		nnLine 		= nnLines[i+1].rstrip()


		origReal 	= float(origLine.split("\t")[0])
		origImag 	= float(origLine.split("\t")[1])

		nnReal 		= float(nnLine.split("\t")[0])
 		nnImag 		= float(nnLine.split("\t")[1])

 		diffReal	= origReal - nnReal
 		diffImag	= origImag - nnImag

 		nominator   = math.sqrt(diffReal*diffReal + diffImag*diffImag)
 		denominator = math.sqrt(origReal*origReal + origImag*origImag)


 		if(denominator == 0):
 			e = 1.0
 		elif(math.isnan(nominator) or (math.isnan(denominator))):
 			e = 1.0
 		elif ((nominator / denominator > 1)):
 			e = 1.0
 		else:
 			e = nominator / denominator

 		absError += e
	pass;
	#print k
	print "*** Error: %1.8f " % (absError/float(len(origLines))) + str(k)
