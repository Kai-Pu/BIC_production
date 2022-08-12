#!/usr/bin/env python

import sys

fileNamePrefix = sys.argv[1]


fileNamePath = "./" + fileNamePrefix + "_Processed_anno.profile"
print(" ---------fileNamePath------------ ")
print(fileNamePath)

seqOTUdict = dict()
seqTaxFreqDict = dict()
with open(fileNamePath, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        database = line[0].split(".")[0]
        if line[4] == "taxonomy" or line[4] == "n/a" or line[4].startswith("Archaea") or line[4].startswith("Bacteria;;;;;;;;"):
            continue
#        if database is None:
#            print("data base is none")
#            print(database)
#            print(line)
        if line[2] not in seqOTUdict:
            # manipulate the otu table
            seqOTUdict[line[2]] = [line[1].split('-')[1], line[4], [database]]
            #print("First seq in dict: " + line[2])
            #print(database)
            #print(seqOTUdict[line[2]])
            # manipulate the sequence frequency
            taxName = line[4].split(';')
            seqTaxFreqDict[line[2]] = [[taxName[0]], [taxName[1]], [taxName[2]], [taxName[3]], [taxName[4]], [taxName[5]], [taxName[6]], [taxName[7]]]
        else:
            # manipulate the otu table
            newTax = line[4].split(';')
            oriTax = seqOTUdict[line[2]][1].split(';')
            taxLevelIndex = 0
            for i in range(8):
                if newTax[i] == oriTax[i]:
                    taxLevelIndex = i
                else:
                    break
            sameTax = []
            for i in range(8):
                if i<= taxLevelIndex:
                    sameTax.append(oriTax[i])
                else:
                    sameTax.append('')
            seqOTUdict[line[2]][1] = ';'.join(sameTax)

            if database not in seqOTUdict[line[2]][2]:
                seqOTUdict[line[2]][2].append(database)
#            try:
#                seqOTUdict[line[2]][2] = seqOTUdict[line[2]][2].append(database)
#            except:
#                print("Error line")
#                print(line)
#                print(database)
#                print(seqOTUdict[line[2]])
            # manipulate the sequence frequency
            for i in range(8):
                if newTax[i] not in seqTaxFreqDict[line[2]][i]:
                 seqTaxFreqDict[line[2]][i].append(newTax[i])
            
            
            
            
            
            
outputPath1 = "./" + fileNamePrefix + "_otuTable.txt"
fp = open(outputPath1, 'w')
fp.write("Sequence\tCount\tTaxonomy\tDatabase"+'\n')
for seq in sorted(seqOTUdict.keys()):
    fp.write('\t'.join([seq]+seqOTUdict[seq][0:2]+[";".join(sorted(seqOTUdict[seq][2]))])+'\n')
fp.close()


outputPath2 = "./" + fileNamePrefix + "_checkSeqMatchBacteriaNum.txt"
fp = open(outputPath2, 'w')
fp.write('\t'.join(["Sequence", "Superkingdom", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])+'\n')

for seq in sorted(seqTaxFreqDict.keys()):
    seqFreq = [seq]
    for i in range(8):
        seqFreq.append(str(len(seqTaxFreqDict[seq][i])))
    fp.write('\t'.join(seqFreq)+'\n')
fp.close()
    
