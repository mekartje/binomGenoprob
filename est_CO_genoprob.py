##Estimate crossover positions from genotype probabilities as positions where Pr(het) goes from tau to 1-tau (Karl's suggestion)
##1st argument -- file with crossover probabilities
##2nd argument -- tau (0-1)

import sys

#function to look at dictionary values, find positions where
#hetProb: dictionary with positions (bp) as keys, pr(het) as values
#tau: probability cutoff
def getCO(hetProb, tau):
    #get sorted list of positions
    posLis = list(hetProb.keys())
    posLis.sort()
    #list of potential genotype switches (P[k] >= tau and P[k+1] < tau | P[k] <= 1 - tau and P[k+1] > 1 - tau)
    genoSwitch = []
    #list of k
    posk = posLis[:-1]
    #list of k + 1
    poskp1 = posLis[1:]
    for k in range(len(posk)):
        if hetProb[posk[k]] >= tau and hetProb[poskp1[k]] < tau:
            genoSwitch.append(posk[k])
        elif hetProb[posk[k]] <= 1 - tau and hetProb[poskp1[k]] > 1 - tau:
            genoSwitch.append(posk[k])
    #for all k satisfying either switch condition, find switch boundary satisfying p[k+i] <= 1 - tau if p[k] >= tau | p[k+i] >= tau if p[k] <= 1 - tau
    #store these intervals in posCO
    #posCO is a list of tuples [(left_border, right_border), ...]
    posCO = []
    for k in genoSwitch:
        #define i as next marker position -- i have chosen a strange way to do this!
        i = 1
        if hetProb[k] >= tau:
            #Explain strange indexing -- need to index hetProb with an actual position in bp, because these are its keys
                #to get the next position, use the list 'posLis'. First find the position for k (posLis.index(k)), and add i to it (posLis.index(k) + 1)
                #then take this entry from posLis (posLis[posLis.index(k) + 1]) to get the index to use for the next marker in hetProb
            while hetProb[posLis[posLis.index(k) + i]] > 1 - tau:
                i += 1
        elif hetProb[k] <= 1 - tau:
            while hetProb[posLis[posLis.index(k) + i]] < tau:
                i += 1
        posCO.append((k, posLis[posLis.index(k) + i]))
    return(posCO)
tau = float(sys.argv[2])
#initialize dictionary of P(het)
hetProb = {}
with open(sys.argv[1]) as inf:
    #read off header line
    inf.readline()
    #add prAB to hetProb with position as key
    for line in inf:
        line = line.rstrip().split(' ')
        hetProb[int(line[0])] = float(line[2])
posCO = getCO(hetProb = hetProb, tau = tau)

with open(sys.argv[1][:-4] + '_posCO.txt', 'w') as outf:
    #write header line
    outf.writelines('CO_start\tCO_end\n')
    #write CO boundaries
    for CO in posCO:
        outf.writelines(str(CO[0]) + '\t' + str(CO[1]) + '\n')
