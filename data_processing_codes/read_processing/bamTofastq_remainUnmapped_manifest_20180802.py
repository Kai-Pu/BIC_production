#!/usr/bin/env python2

######################################################
## input: manifest file, inputPath, outputPath
## usage: python bamTofastq_remainUnmapped_manifest_20180802.py ~/dir/manifest ~/inputDir/ ~/outFastqDir/
## output: uuid_unmapped.fastq
## aim: turn bam file into fastq format without mapped reads with samtools
## samtools version: 1.3.1
## author: Kai-Pu Chen
## date: 20190802
######################################################

import os
import sys

manifestPath = sys.argv[1]
inBamRootDir = sys.argv[2]
outFastqDir = sys.argv[3]

print "Manifest: " + manifestPath
print "inBamRootDir: " + inBamRootDir
print "outFastqDir: " + outFastqDir


fileList = []

try:
    with open(manifestPath, 'r') as f:
        for line in f:
            if not line.startswith("id"):
                line = line.strip().split('\t')
                #print line
                bamPath = inBamRootDir + '/' + line[0] + '.bam'
                fastqPath = outFastqDir + '/' + line[0] + '.fastq'
                fileList.append([bamPath, fastqPath])

except:
    print "Error in reading file!\nPlease check if InputPath and manifest files are correct!"


#print "Head 5 lines are:"
#for i in range(0,5):
#    print fileList[i][0] + '\t' + fileList[i][1]



for i in range(0, len(fileList)):
    cmd = "samtools fastq " + fileList[i][0] + " >" + fileList[i][1]
    stderrMsgCMD = str(i) + '......' + cmd + '......'
    print stderrMsgCMD
    #print >> sys.stderr, stderrMsgCMD

    cmdStateCode = os.system(cmd)
    stderrMsgSTATE_OK = fileList[i][1] + 'has been done!'
    stderrMsgSTATE_FAIL = "Error in transfering the " + fileList[i][1]
    if cmdStateCode == 0:
        print stderrMsgSTATE_OK
        #print >> sys.stderr, stderrMsgSTATE_OK
    else:
        print stderrMsgSTATE_FAIL
        #print >> sys.stderr, stderrMsgSTATE_FAIL
    print '\n\n'
