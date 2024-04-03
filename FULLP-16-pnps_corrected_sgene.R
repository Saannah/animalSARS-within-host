#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## set dir
dir = "7-pnps_corrected_sgene/"

### !!!this script makes the assumption that no site belongs to more than two ORFs.
#rm(list = ls())
library(phylotools)
library(stringr)

### load the reference sequence with phylotools
ref = read.fasta("ref_from_genbank.fasta")

library(seqinr)




#### load file names
vcfnames = read.delim("parsednames", header = FALSE)
vcfnames = vcfnames[, 1]
## load one file 
vcf_filename = vcfnames[as.numeric(args[1])]

## extract sra id from the vcf file name
SRA_id = unlist(str_split(vcf_filename, "\\."))[1]
print(SRA_id)

## check if the file is empty
if(file.size(vcf_filename) == 0){
  ### create dataframe for output
  out = data.frame(SRA_id = SRA_id,
                   pNpS_sgene = 1,
                   number_of_snps_filtered = 0,
                   correction = "yes")
  
  write.table(out, file = paste0(dir, SRA_id, "_pNpS_sgene_corrected.tsv"), sep = "\t", row.names = FALSE)
  print("file written")
  stop("file empty")
}








### define ORF boundaries:
start = c(21563)
end = c(25384)
name = c("S")


#### load the VCF
vcf = read.delim(vcf_filename)
### filter VCF based on strict criteria defined
vcf = vcf[which(vcf$VAF > 0.05), ]
vcf = vcf[which(vcf$DEPTH >= 50), ]


### remove indels
if(length(which(nchar(vcf$REF) > 1)) > 0){
  vcf = vcf[-which(nchar(vcf$REF) > 1), ]
}

if(length(which(nchar(vcf$ALT) > 1)) > 0){
  vcf = vcf[-which(nchar(vcf$ALT) > 1), ]
}



## convert everything to lowercase
vcf$REF = tolower(vcf$REF)
vcf$ALT = tolower(vcf$ALT)


#### define ORFs to consider, only coding regions (so the start and stop codpns are removed becuase they are not counted in the calculation of pnps)
ORFs = data.frame(start = start + 3,
                  end = end - 3,
                  name = name)
rownames(ORFs) = ORFs$name


ORFs_rowwise = rbind(start + 3,
                     end - 3)
colnames(ORFs_rowwise) = name




### create an array of the reference bases
ref_array = tolower(unlist(str_split(ref$seq.text, "")))



if(nrow(vcf) == 0){
  ### create dataframe for output
  out = data.frame(SRA_id = SRA_id,
                   pNpS_sgene = 1,
                   number_of_snps_filtered = nrow(vcf),
                   correction = "yes")
  write.table(out, file = paste0(dir, SRA_id, "_pNpS_sgene_corrected.tsv"), sep = "\t", row.names = FALSE)
  print("file written")
  stop("vcf had no eligible SNPS after filtering")
}

## find the index of mutations that are major allele
ind = which(vcf$VAF >= 0.5)
flipped_ref_array = ref_array
flipped_ref_array[vcf$POS[ind]] = tolower(vcf$ALT[ind])

## create a column for flipped SNPs
vcf$flipped_alt_allele = vcf$ALT
vcf$flipped_alt_allele[ind] = vcf$REF[ind]
vcf$flipped_ref_allele = vcf$REF
vcf$flipped_ref_allele[ind] = vcf$ALT[ind]


## create a column to count the number of synonymous and non-synonymous mutations
vcf$s = 0
vcf$n = 0


### create a function that finds the ORF of a site, and returns it. this only looks at coding sites so sites in start or stop codons are not considered
find_ORF_of_site = function(site){
  orfs_with_site = vector()
  for(i in 1:nrow(ORFs)){
    if(site >= ORFs$start[i]){
      if(site <= ORFs$end[i]){
        orfs_with_site = c(orfs_with_site, ORFs$name[i])
      }
    }
  }
  if(length(orfs_with_site) == 0){
    return("IG")
  }
  return(orfs_with_site)
}

find_codon_in_ref = function(site, ref_to_use){
  orfs_with_site = find_ORF_of_site(site)
  
  if("IG" %in% orfs_with_site){
    return(data.frame(site = vector(),
                      codons = vector(),
                      positions_in_each_codon = vector(),
                      ORFs = vector()))
  }
  
  codons_in_reftouse = vector()
  positions_in_each_codon = vector()
  
  for(orf in orfs_with_site){
    
    pos_in_codon = (site - (ORFs_rowwise[1, orf] - 1)) %% 3
    positions_in_each_codon = c(positions_in_each_codon, pos_in_codon)
    
    
    ## if the position is 1
    if(pos_in_codon == 1){
      codon = paste0(ref_to_use[site:(site + 2)], collapse = "")
      codons_in_reftouse = c(codons_in_reftouse, codon)
      next
    }
    
    ## if the position in the codon is 2
    if(pos_in_codon == 2){
      codon = paste0(ref_to_use[(site - 1):(site + 1)], collapse = "")
      codons_in_reftouse = c(codons_in_reftouse, codon)
      next
    }
    
    ## if the position in the codon is 3 (modulo is 0)
    if(pos_in_codon == 0){
      codon = paste0(ref_to_use[(site - 2):site], collapse = "")
      codons_in_reftouse = c(codons_in_reftouse, codon)
    }
  }
  
  positions_in_each_codon[which(positions_in_each_codon == 0)] = 3
  df = data.frame(site = rep(site, length(orfs_with_site)),
                  codons = codons_in_reftouse,
                  positions_in_each_codon = positions_in_each_codon,
                  ORFs = orfs_with_site)
  return(df)
  
}

### write a function that returns a dataframe with "s" and "n" columns, indicating whether a mutation is synonymous (s=1), non-synonymous (n=1) or neither (s=0, n=0)
### if from_codon is a start codons or stop codons then they are counted as 0

s_or_n = function(from_codon, to_codon){
  from_translated = translate(unlist(str_split(from_codon, "")))
  # if(from_translated == "*"){
  #   s_n_count = (c(0, 0))
  #   names(s_n_count) = c("s", "n")
  #   return(s_n_count)
  # }
  
  to_translated = translate(unlist(str_split(to_codon, "")))
  
  if(from_translated == to_translated){
    s_n_count = (c(1, 0))
    names(s_n_count) = c("s", "n")
    return(s_n_count)
  }
  
  if(from_translated != to_translated){
    s_n_count = (c(0, 1))
    names(s_n_count) = c("s", "n")
    return(s_n_count)
  }
  
}



#### count the number of synonymous and non-synonymous mutations

## for 1 row of vcf
find_s_or_n_mutation_for_1_row = function(vcf, rownum){
  s = 0
  n = 0
  codons = find_codon_in_ref(vcf$POS[rownum], flipped_ref_array)
  
  # if the mutation is not in a codon, then that mutation is not counted in the calculation and the sums in the VCF do not change
  if(nrow(codons) == 0){
    counts = data.frame(s = 0,
                        n = 0)
    
    return(counts)
  }
  
  
  # if the mutation is in 1 or more codons, go through each codon and determine if the product is synonymous or non-synonymous
  for(i in 1:nrow(codons)){
    new_alt = vcf$flipped_alt_allele[rownum]
    from_codon = unlist(str_split(codons$codons[i], ""))
    to_codon = from_codon
    to_codon[codons$positions_in_each_codon[i]] = new_alt
    
    ## if the product of both is the same, add 1 to synonymous tally
    if(translate(from_codon) == translate(to_codon)){
      s = s + 1
    }
    
    ## if the product of both is not the same, add 1 to non-synonymous count
    if(translate(from_codon) != translate(to_codon)){
      n = n + 1
    }
  }
  
  counts = data.frame(s = s,
                      n = n)
  
  return(counts)
  
}

vcf_mutation_tallies = do.call(rbind,
                               lapply(1:nrow(vcf), find_s_or_n_mutation_for_1_row, vcf = vcf))

number_synonymous_mutations = sum(vcf_mutation_tallies$s)
number_nonsynonymous_mutations = sum(vcf_mutation_tallies$n)

####!!!!!!!!!! when counting SITES only look between start and stop codon (remove positions from the ranges) AND dilter VCF/sites based on this)


#### aclculate the number of synonymoud and non-synonymous sites
#### for this, go through every ORF and count accordingly for that ORF
#### this automatically accounts for overlapping ORFs with frameshifts
### the function is run over ALL ORFs separately, so for a site in multiple oRFs, we only consider the codon in one orf at a time

all_bases = c("a", "t", "c", "g")
count_for_one_site = function(site, orf_name){
  s = 0
  n = 0
  codons = find_codon_in_ref(site, flipped_ref_array)
  from_codon = unlist(str_split(codons[orf_name, "codons"], ""))
  pos_in_codon = codons[orf_name, "positions_in_each_codon"]
  from_base = from_codon[pos_in_codon]
  to_bases = all_bases[-which(all_bases == from_base)]
  
  to_codon = from_codon
  for(alt_base in to_bases){
    to_codon[pos_in_codon] = alt_base
    
    if(translate(from_codon) == translate(to_codon)){
      s = s + 1
      next
    }
    
    if(translate(from_codon) != translate(to_codon)){
      n = n + 1
    }
  }
  
  s = s/3
  n = n/3
  
  return(data.frame(site = site,
                    orf_name = orf_name,
                    s_count = s,
                    n_count = n))
}




get_orf_counts = function(orf_name){
  orf_counts = do.call(rbind,
                       lapply(
                         ORFs[orf_name, "start"]:ORFs[orf_name, "end"], 
                         count_for_one_site, 
                         orf_name = orf_name)
  )
  
  return(data.frame(orf_name = orf_name,
                    s = sum(orf_counts$s_count),
                    n = sum(orf_counts$n_count)))
}



all_orf_counts = do.call(rbind,
                         lapply(name, get_orf_counts))

number_synonymous_sites = sum(all_orf_counts$s)
number_nonsynonymous_sites = sum(all_orf_counts$n)

if(number_synonymous_mutations == 0){
  pnps = ((number_nonsynonymous_mutations/number_nonsynonymous_sites) + 1)/((number_synonymous_mutations/number_synonymous_sites) + 1)
  correction = "yes"
} else{
  pnps = (number_nonsynonymous_mutations/number_nonsynonymous_sites)/(number_synonymous_mutations/number_synonymous_sites)
  correction = "no"
}

out = data.frame(SRA_id = SRA_id,
                 pNpS_sgene = pnps,
                 number_of_snps_filtered = nrow(vcf),
                 correction = correction)
write.table(out, file = paste0(dir, SRA_id, "_pNpS_sgene_corrected.tsv"), sep = "\t", row.names = FALSE)
print("file written")




