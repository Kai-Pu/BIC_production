#!/usr/bin/env python

import sys

TCGAProject = sys.argv[1]
dataDir = sys.argv[2]
fileNamePrefixListPath = sys.argv[3]
taxaLevel = sys.argv[4]
matchNum = 1

if dataDir[-1] != '/':
    dataDir = dataDir +'/'
    

fileNamePrefixList = []
with open(fileNamePrefixListPath, 'r') as f:
    for line in f:
        line = line.strip()
        fileNamePrefixList.append(line)
print(TCGAProject + " fileNum: "+ str(len(fileNamePrefixList)))


countDict = {}
for i, uuid in enumerate(fileNamePrefixList):
    print(i, uuid)
    filePath = dataDir + uuid + '/' + uuid + '_' + taxaLevel + '_matchNum_' + str(matchNum) + "_bacteria_seqFreq_count_selectedDB.txt"
    with open(filePath, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[0]=="bacteriaName":
                continue            
            bacteria = line[0]
            count = line[-2]
            if bacteria not in countDict:
                countDict[bacteria] = {}
            if i not in countDict[bacteria]:
                countDict[bacteria][i] = count

countList = [[taxaLevel] + fileNamePrefixList]
for k1,v1 in sorted(countDict.items()):
    tmp = [str(0)]*len(fileNamePrefixList)
    for k2,v2 in sorted(v1.items()):
        tmp[k2] = v2
    countList.append([k1] + tmp)



outputPath = dataDir + "TCGA_" + TCGAProject + "_miRNAseq_" + str(len(fileNamePrefixList)) + "samples_" + taxaLevel + "_matchNum_" + str(matchNum) + "_combineCounts_selectedDB.txt"
fp = open(outputPath, 'w')
for line in countList:
    fp.write('\t'.join(line) + '\n')
fp.close()
