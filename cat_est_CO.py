#Concatenate crossover intervals from est_CO_genoprob.py (ending in *posCO.txt)
##1st argument -- list of CO files
##2nd argument -- list of output file

import sys

with open(sys.argv[1]) as inf:
    flis = [line.rstrip() for line in inf]
with open(sys.argv[2], 'w') as outf:
    #write header line
    outf.writelines('rep\tCO_start\tCO_end\n')
    for f in flis:
        inf = open(f)
        #read off header line
        inf.readline()
        #get rep ID
        rep = '_'.join(f.split('_')[:4])
        #read remaining res to outfile
        for line in inf:
            line = line.rstrip().split('\t')
            outf.writelines(rep + '\t' + line[0] + '\t' + line[1] + '\n')
        inf.close()
