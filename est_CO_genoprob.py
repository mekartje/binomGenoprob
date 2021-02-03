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

#some of the output from getCO will be redundant, making output like this:
#CO_start        CO_end
#90209440        101005788
#100305541       101005788

#function to refine CO intervals, keeping only the smallest. potentially come back to refine this later
#redundant intervals will only ever share an endpoint, b/c the loop in getCO moves from left to right
def refineCO(posCO):
    #cast posCO as dictionary with endpoint as key, start as value
    posCO_dict = {}
    for i in posCO:
        #if endpoing alread in posCO
        if i[1] in list(posCO_dict.keys()):
            #find which startpoint gives a smaller interval
            oldLength = i[1] - posCO_dict[i[1]]
            newLength = i[1] - i[0]
            #keep shorter interval
            if oldLength < newLength:
                continue
            elif newLength < oldLength:
                posCO_dict[i[1]] = i[0]
        else:
            posCO_dict[i[1]] = i[0]
    #rewrite posCO as list of tuples. inefficient, but this is the was things were set up originally
    posCO_new = []
    for k,v in posCO_dict.items():
        posCO_new.append((v, k))
    return(posCO_new)

#main
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
#refine CO interval list with refineCO
posCO = refineCO(posCO)
#write output
with open(sys.argv[1][:-4] + '_posCO.txt', 'w') as outf:
    #write header line
    outf.writelines('CO_start\tCO_end\n')
    #write CO boundaries
    for CO in posCO:
        outf.writelines(str(CO[0]) + '\t' + str(CO[1]) + '\n')
