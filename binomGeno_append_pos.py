#script appending vcf positions to binomGenoprob output
##1st argument -- resfile
##2nd argument -- vcf

import sys

#get list containing positions
pos_ls = []
with open(sys.argv[2]) as vcf:
    for line in vcf:
        if not line.startswith('#'):
            pos_ls.append(line.rstrip().split('\t')[1])

#make new resfile with pos column
with open(sys.argv[1]) as resf, open(sys.argv[1][:-4] + '_wpos.txt', 'w') as outf:
    #read off resf headerline
    resf.readline()
    #write outf headerline
    outf.writelines('position\tprAA\tprAB\n')
    counter = 0
    for line in resf:
        line = line.rstrip().split(' ')
        outf.writelines(pos_ls[counter] + ' ' + line[1] + ' ' + line[2] + '\n')
        counter += 1
        
