#!/usr/bin/env python

import sys


fileNamePrefix = sys.argv[1]

matchNum = 1
dataDir = "./"


taxaLevelList = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]




    
fileNamePath_otuTable = dataDir + fileNamePrefix + "_otuTable.txt"
print(" ---------fileNamePath_otuTable------------ ")
print(fileNamePath_otuTable)

bacteriaUniqueSeqBacteria = [{}, {}, {}, {}, {}, {}]
with open(fileNamePath_otuTable, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        
        
        if line[0] == "Sequence":
            continue
        
        seq = line[0]
        count = int(line[1])
        taxaList=line[2].split(";")[2:]
        database = line[3].split(";")
        
        for i in range(6):
            if taxaList[i] != '':
                if taxaList[i] not in bacteriaUniqueSeqBacteria[i]:
                    bacteriaUniqueSeqBacteria[i][taxaList[i]] = [[seq], [str(count)], 1, count, database]
                else:
                    if seq not in bacteriaUniqueSeqBacteria[i][taxaList[i]][0]:
                        bacteriaUniqueSeqBacteria[i][taxaList[i]][0].append(seq)
                        bacteriaUniqueSeqBacteria[i][taxaList[i]][1].append(str(count))
                        bacteriaUniqueSeqBacteria[i][taxaList[i]][2] = bacteriaUniqueSeqBacteria[i][taxaList[i]][2] + 1
                        bacteriaUniqueSeqBacteria[i][taxaList[i]][3] = bacteriaUniqueSeqBacteria[i][taxaList[i]][3] + count
                        for db in database:
                            if db not in bacteriaUniqueSeqBacteria[i][taxaList[i]][4]:
                                bacteriaUniqueSeqBacteria[i][taxaList[i]][4].append(db)
                    else:
                        for db in database:
                            if db not in bacteriaUniqueSeqBacteria[i][taxaList[i]][4]:
                                bacteriaUniqueSeqBacteria[i][taxaList[i]][4].append(db)

    
        

for i in range(len(taxaLevelList)):
    outputPath = dataDir + fileNamePrefix + '_' + taxaLevelList[i] + "_count_table" + ".txt"
    fp = open(outputPath, 'w')
    fp.write('\t'.join(["bacteriaName", "seq", "seq_individual_count", "seq_num", "count", "database"])+'\n')
    for bacteria in sorted(bacteriaUniqueSeqBacteria[i].keys()):
        fp.write('\t'.join([bacteria] + [';'.join(bacteriaUniqueSeqBacteria[i][bacteria][0])] + [';'.join(bacteriaUniqueSeqBacteria[i][bacteria][1])] + [str(bacteriaUniqueSeqBacteria[i][bacteria][2]), str(bacteriaUniqueSeqBacteria[i][bacteria][3])] + [';'.join(sorted(bacteriaUniqueSeqBacteria[i][bacteria][4]))])+'\n')
    fp.close()



                    

