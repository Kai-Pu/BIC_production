#!/usr/bin/env python

import sys


fileNamePrefix = sys.argv[1]
taxaLevel = sys.argv[2]

matchNum = 1
dataDir = "./"


if taxaLevel == "Species":
    taxaIndex = 7
    levelIndex = 8
elif taxaLevel == "Genus":
    taxaIndex = 6
    levelIndex =7
elif taxaLevel == "Family":
    taxaIndex = 5
    levelIndex =6
elif taxaLevel == "Order":
    taxaIndex = 4
    levelIndex =5
elif taxaLevel == "Class":
    taxaIndex = 3
    levelIndex =4
elif taxaLevel == "Phylum":
    taxaIndex = 2
    levelIndex =3




print("taxaLevel = " + taxaLevel)





    
fileNamePath_checkSeqMatchBacteriaNum = dataDir + fileNamePrefix + "_checkSeqMatchBacteriaNum.txt"
print(" ---------fileNamePath_checkSeqMatchBacteriaNum------------ ")
print(fileNamePath_checkSeqMatchBacteriaNum)

bacteriaSpeciesSeq = []    
with open(fileNamePath_checkSeqMatchBacteriaNum, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        if line[0] == "Sequence":
            continue

        if int(line[levelIndex]) == int(matchNum):
            bacteriaSpeciesSeq.append(line[0])
            

fileNamePath_anno_profile = dataDir + fileNamePrefix + "_Processed_anno.profile"
print(" ---------fileNamePath_anno_profile------------ ")
print(fileNamePath_anno_profile)

bacteriaList = []
annotPartInputList = []
with open(fileNamePath_anno_profile, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        if line[4]=="taxonomy" or line[4]=="n/a" or line[4].startswith("Archaea") or line[4].startswith("Bacteria;;;;;;;;"):
            continue
        if line[2] in bacteriaSpeciesSeq:
            taxList = line[4].split(';')
            for i in range(len(taxList)):
                if taxList[i] == '':
                    taxList[i] = "unknown"
            
            bacteria = taxList[taxaIndex]
            annotPartInputList.append((line[0:3]+[bacteria]))
            
            if bacteria not in bacteriaList:
                bacteriaList.append(bacteria)



bacateriaSeqFreqCountsDict = {}
for i in range(len(annotPartInputList)):
    database = annotPartInputList[i][0].split('.')[0]
    count = int(annotPartInputList[i][1].split('-')[1])
    seq = annotPartInputList[i][2]
    
    if annotPartInputList[i][3] not in bacateriaSeqFreqCountsDict:
        bacateriaSeqFreqCountsDict[annotPartInputList[i][3]] = [[seq], [str(count)], 1, count, [database]]
    
    else:
        if seq not in bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][0]:
            bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][0].append(seq)
            bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][1].append(str(count))
            bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][2] = bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][2] + 1
            bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][3] = bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][3] + count
            if database not in bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][4]:
                bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][4].append(database)
        else:
            if database not in bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][4]:
                bacateriaSeqFreqCountsDict[annotPartInputList[i][3]][4].append(database)
            
        

        
outputPath = dataDir + fileNamePrefix + '_' + taxaLevel + '_matchNum_' + str(matchNum) + "_bacteria_seqFreq_count_selectedDB" + ".txt"
fp = open(outputPath, 'w')
fp.write('\t'.join(["bacteriaName", "seq", "seq_individual_count", "seq_num", "count", "database"])+'\n')
for bacteria in sorted(bacateriaSeqFreqCountsDict.keys()):
    fp.write('\t'.join([bacteria] + [';'.join(bacateriaSeqFreqCountsDict[bacteria][0])] + [';'.join(bacateriaSeqFreqCountsDict[bacteria][1])] + [str(bacateriaSeqFreqCountsDict[bacteria][2]), str(bacateriaSeqFreqCountsDict[bacteria][3])] + [';'.join(sorted(bacateriaSeqFreqCountsDict[bacteria][4]))])+'\n')
fp.close()



                    

