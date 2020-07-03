#Estimating genotype probabilities from allele count data

args<-commandArgs(trailingOnly = TRUE)
#Functions for initial, transition and emission probabilities

#initiation probabilities
#for backcross, always return 0.5, from exp transmission ratio
init<-function(){return(log(0.5))}

#emission probabilities
emit<-function(n_ref_read, #integer -- observed count of reference allele
  n_tot_read, #total number of reads
  true_gen, #integer -- true genotype
  error_prob) { #double -- error probability
    #assume autosome
    if(true_gen == 'AA'){ #if true genotype homozygote, return prob (X=n_ref_read) where X~binom(n = n_tot, p = 1 - error_prob)
      em_pr<-dbinom(x = n_ref_read, size = n_tot_read, prob = 1 - error_prob, log = TRUE)
    }
    else if(true_gen == 'AB'){ #if true genotype heterozygote, return prob (X=n_ref_read) where X~binom(n = n_tot, p = 0.5)
      em_pr<-dbinom(x = n_ref_read, size = n_tot_read, prob = 0.5, log = TRUE)
    }
    return(em_pr)
}

#step probabilities
step<-function(gen_left, gen_right, rec_frac){
  #if no change in genotype, return prob of no recombination
  if(gen_left == gen_right){return(log(1 - rec_frac))}
  #if change in genotype, return prob of recombination
  else if(gen_left != gen_right){return(log(rec_frac))}
}

#function adding logs (log-transform probs to prevent underflow)
addlog<-function(a, b){
  if(is.nan(a) || is.nan(b)){return(NaN)}
  else{
    tol = 200
    if(b > (a + tol)){return(b)}
    else if(a > b + tol){return(a)}
    else{return(a + log1p(exp(b - a)))}
  }
}


#Functions for forward and backward equations

forwardEquations<-function(ref_read_ns, tot_read_ns, rec_frac, error_prob, poss_gen){
  n_pos<-length(ref_read_ns)
  n_gen<-length(poss_gen)
  #matrix containing alphas (ai(g) = P(O1, O2, O3, ..., Oi, Gi = g))
  alpha<-matrix(data = 0, nrow = n_gen, ncol = n_pos)
  #initialize alphas
  #for each possible genotype (2 possibilities in backcross) (i.e., iterate over each row in the matrix)
  for(i in 1:n_gen){
    #genotype of current iteration
    g<-poss_gen[i]
    #initialize 1st value in row i with initiation probability
    alpha[i,1]<-init()
    #add to this the emission probability of the observed reference read count given genotype of current iteration
    alpha[i,1]<-alpha[i,1] + emit(n_ref_read = ref_read_ns[1], n_tot_read = tot_read_ns[1], true_gen = g, error_prob = error_prob)
  }
  #for each position (pos; column) starting at position 2, and each possible genotype (ir)
  for(pos in 2:n_pos){
    for(ir in 1:n_gen){
      #initialize alpha for genotype ir, position pos
      alpha[ir, pos]<-alpha[1, pos - 1] + step(poss_gen[1], poss_gen[ir], rec_frac[pos - 1])
      #iterating over remaining genotypes, il
      for(il in 2:n_gen){
        #add to this the following sum:
        #alpha for genotype il at the position immediately left + the step prob from genotype il to genotype ir at the interval to the left of position pos
        alpha[ir, pos]<-addlog(alpha[ir, pos], alpha[il, pos - 1] + step(poss_gen[il], poss_gen[ir], rec_frac[pos - 1]))
      }
      #once both genotypes are iterated over, add the emission probability of the observed number of reference reads given true genotype ir
      alpha[ir, pos]<-alpha[ir, pos] + emit(n_ref_read = ref_read_ns[pos], n_tot_read = tot_read_ns[pos], true_gen = poss_gen[ir], error_prob = error_prob)
    }
  }
  return(alpha)
}

backwardEquations<-function(ref_read_ns, tot_read_ns, rec_frac, error_prob, poss_gen){
  n_pos<-length(ref_read_ns)
  n_gen<-length(poss_gen)
  #matrix containing betas (bi(g) = P(Oi+1, Oi+2, On | Gi = g)) where n = number of markers
  beta<-matrix(data = 0, nrow = n_gen, ncol = n_pos)
  for(pos in n_pos - 1:1){
    for(il in 1:n_gen){
      for(ir in 1:n_gen){
        to_add<-beta[ir, pos + 1] + step(poss_gen[il], poss_gen[ir], rec_frac[pos])
        to_add<-to_add + emit(n_ref_read = ref_read_ns[pos + 1], n_tot_read = tot_read_ns[pos + 1], true_gen = poss_gen[ir], error_prob)
        if(ir == 0){beta[ir, pos]<-to_add}
        else{beta[il, pos]<-addlog(beta[il, pos], to_add)}
      }
    }
  }
  return(beta)
}

#ref_read_ns and tot_read_ns are matricies with columns as individuals, rows as markers
#backcross individuals only, so poss_gen is the same for all individuals -- AA or AB
#copying rqtl2 hmm_calcgenoprob.cpp, except for structure of output array
calc_genoprob<-function(ref_read_ns, tot_read_ns, rec_frac, error_prob, poss_gen = c('AA', 'AB')){
  n_ind<-ncol(ref_read_ns)
  n_pos<-nrow(ref_read_ns)
  n_mar<-n_pos
  n_poss_gen<-length(poss_gen)
  n_gen<-2 #assuming autosome. 4 if X
  #initialize results array -- rows are gen, cols are pos, matrix for each indv
  genoprobs<-array(data = 0, dim = c(n_gen, n_pos, n_ind))
  #for each individual, get alphas and betas
  for(ind in 1:n_ind){
    alpha<-forwardEquations(ref_read_ns = ref_read_ns[,ind], tot_read_ns = tot_read_ns[,ind], rec_frac = rec_frac, error_prob = error_prob, poss_gen = poss_gen)
    beta<-backwardEquations(ref_read_ns = ref_read_ns[,ind], tot_read_ns = tot_read_ns[,ind], rec_frac = rec_frac, error_prob = error_prob, poss_gen = poss_gen)
    for(pos in 1:n_pos){
      g<-1 #starts at 0 for BC assuming autosomes only, see genotype enumeration in cross_bc.cpp
      sum_at_pos<-genoprobs[g,pos,ind]<-alpha[1,pos] + beta[1,pos]
      #loop over remaining genotypes -- change(i in 1:n_poss_gen) to (i in 2:n_poss_gen)
      for(i in 2:n_poss_gen){
        g<-i #just putting this here for consistency with cpp code, could just use 'i' to index genoprobs
        val<-genoprobs[g,pos,ind]<-alpha[i, pos] + beta[i, pos]
        sum_at_pos<-addlog(sum_at_pos, val)
      }
      for(i in 1:n_poss_gen){
        g<-i #same as above
        genoprobs[g,pos,ind]<-exp(genoprobs[g,pos,ind] - sum_at_pos)
      }
    }
  }
  return(genoprobs)
}

#function to pass data table to calc_genoprob. takes a VCF
#gwrr -- genome-wide recombination rate, in cM/Mb
#assumes that first count in AD vector is ref (i.e., ID of allele present in homozygotes).
pass_dat<-function(path, rec_frac_unif = TRUE, error_prob, gwrr){
    error_prob<-as.numeric(error_prob)
    gwrr<-as.numeric(gwrr)
    dat<-read.table(path, header = FALSE)
    ref_read_ns<-numeric(nrow(dat))
    tot_read_ns<-numeric(nrow(dat))
    #parse allele count data
    for(i in 1:nrow(dat)){
      print(as.character(dat[i]))
      allele_counts<-as.numeric(strsplit(strsplit(as.character(dat[[10]][i]), split = ':')[[1]][2], split = ',')[[1]])
      #again, need to adjust to assign parental origin of alleles. For now, treat 1/2 cases as 0/1.
      if(length(allele_counts) > 2){
        allele_counts<-allele_counts[2:3]
      }
      ref_read_ns[i]<-allele_counts[1]
      tot_read_ns[i]<-sum(allele_counts)
    }
    #obtain recombination fractions from genome-wide avg (espressed in cM/Mb) assuming uniform rate -- move this into above for loop
    rec_frac<-numeric(nrow(dat) - 1)
    for(pos in 1:length(rec_frac)-1){
      #distance (in bp) between marker at pos and marker at pos + 1
      bp_dist<-as.integer(dat[[2]][pos + 1]) - as.integer(dat[[2]][pos])
      #convert to megabases
      bp_dist<-bp_dist/1e6
      #convert to cM
      cM<-gwrr * bp_dist
      #convert to M
      rec_frac[pos]<-cM/100
    }
    calc_genoprob(ref_read_ns = matrix(ref_read_ns, ncol = 1), tot_read_ns = matrix(tot_read_ns, ncol = 1), rec_frac = rec_frac, error_prob = error_prob)
}

res<-pass_dat(path = args[1], error_prob = args[3], gwrr = args[4])

#coerce single individual data into dataframe
df<-t(as.data.frame(res[,,1]))
write.table(df, file = args[2])
