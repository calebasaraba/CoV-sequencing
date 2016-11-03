source("http://www.bioconductor.org/biocLite.R")
biocLite("msa")

library(msa)
library(Biostrings)
library(seqinr)



x1 <- readDNAStringSet("CoV_sequences.fasta")

ptm <- proc.time()
results <- rep(888, length(x1))
for(i in 1:20){
  seq1 <- x1[1]
  seq2 <- x1[i]
  pair <- append(seq1, seq2)
  clustal1 <- msaClustalW(pair, substitutionMatrix = "iub")
  cl2 <- msaConvert(clustal1, type="seqinr::alignment")
  dist<-dist.alignment(cl2, "identity")
  results[i] <- dist[1]
}

elapsed <- proc.time() - ptm

seq1 <- x1[1]

compare_one <- function(x,y){
  pair <- append(seq1, x)
  clustal1 <- msaClustalW(pair, substitutionMatrix = "iub")
  cl2 <- msaConvert(clustal1, type="seqinr::alignment")
  dist<-dist.alignment(cl2, "identity")
  dist[1]
}

pairwise_compare<- function(secondary, primary){
  pair <- append(primary, secondary)
  clustal_result <- msaClustalW(pair)
  converted_result <- msaConvert(clustal_result, type = "seqinr::alignment")
  distance_matrix <- dist.alignment(converted_result, "identity")
  dist[1]
}

pairwise_compare(x1[1], x1[2])

sequence_subset <- x1[1:20]

saveRDS(sequence_subset, "subset.rds")

lapply(unsequence_subset, pairwise_compare, primary = x1[1])
