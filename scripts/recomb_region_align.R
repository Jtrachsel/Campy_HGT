source('./scripts/03_HGT_ID.R')


plot_islands(RESULTS[[1]], clust_level = quat_cluster)
plot_islands(RESULTS[[2]], clust_level = quat_cluster)
plot_islands(RESULTS[[3]], clust_level = quat_cluster)
plot_islands(RESULTS[[4]], clust_level = quat_cluster)
plot_islands(RESULTS[[5]], clust_level = quat_cluster)
plot_islands(RESULTS[[6]], clust_level = quat_cluster) # I bet sampleF genome isnt complete
plot_islands(RESULTS[[7]], clust_level = quat_cluster)
plot_islands(RESULTS[[8]], clust_level = quat_cluster)
plot_islands(RESULTS[[9]], clust_level = quat_cluster)
plot_islands(RESULTS[[10]], clust_level = quat_cluster)
plot_islands(RESULTS[[11]], clust_level = quat_cluster) # Interesting... ALIGN THIS ONE
plot_islands(RESULTS[[12]], clust_level = quat_cluster) # Interesting... AND THIS ONE # ANNOTATE TetO
plot_islands(RESULTS[[13]], clust_level = quat_cluster) # Interesting... AND THIS ONE # ANNOTATE TetO location
plot_islands(RESULTS[[14]], clust_level = quat_cluster) # ERROR, empty tibble


LOOK <- RESULTS[[13]]


aln1 <- RESULTS[[11]] %>%
  group_by(seqid) %>% 
  summarise(chrom=seqid, 
            start=min(start), 
            end  =max(end)) %>%
  ungroup() %>%
  unique() %>% 
  select(chrom, start, end)

plot_islands(RESULTS[[12]] %>%filter(start<500000), clust_level = quat_cluster)
plot_islands(RESULTS[[12]] %>%filter(start>500000), clust_level = quat_cluster)

aln2A <- RESULTS[[12]] %>%filter(start<500000) %>% 
  group_by(seqid) %>% 
  summarise(chrom=seqid, 
            start=min(start), 
            end  =max(end)) %>%
  ungroup() %>%
  unique() %>% 
  select(chrom, start, end)

aln2B <- RESULTS[[12]] %>%filter(start>500000 & seqid !='13150_1') %>%
  group_by(seqid) %>% 
  summarise(chrom=seqid, 
            start=min(start), 
            end  =max(end)) %>%
  ungroup() %>%
  unique() %>% 
  select(chrom, start, end)

  


plot_islands(RESULTS[[13]] %>%filter(start<250000), clust_level = quat_cluster)
plot_islands(RESULTS[[13]] %>%filter(start>250000), clust_level = quat_cluster)


aln3A <- RESULTS[[13]] %>% filter(start < 250000) %>% 
  group_by(seqid) %>% 
  summarise(chrom=seqid, 
            start=min(start), 
            end  =max(end)) %>%
  ungroup() %>%
  unique() %>% 
  select(chrom, start, end)

aln3B <- RESULTS[[13]] %>% filter(start > 250000) %>% 
  group_by(seqid) %>% 
  summarise(chrom=seqid, 
            start=min(start), 
            end  =max(end)) %>%
  ungroup() %>%
  unique() %>% 
  select(chrom, start, end)

library(msa)
# BiocManager::install('msa')

vector_of_fasta_paths <- list.files(path = './outputs/annotations/', pattern = '.fna', recursive = TRUE, full.names = TRUE)


island_seqs_list <- sapply(vector_of_fasta_paths, readDNAStringSet)
names(island_seqs_list) <- NULL
island_seqs_set <- do.call(c, island_seqs_list)


library(BSgenome)

aln1_seqs <- getSeq(island_seqs_set, as(aln1, 'GRanges'))

aln2A_seqs <- getSeq(island_seqs_set, as(aln2A, 'GRanges'))
aln2B_seqs <- getSeq(island_seqs_set, as(aln2B, 'GRanges'))
aln3A_seqs <- getSeq(island_seqs_set, as(aln3A, 'GRanges'))
aln3B_seqs <- getSeq(island_seqs_set, as(aln3B, 'GRanges')) # rev comp

names(aln1_seqs) <- paste(aln1$chrom, paste(aln1$start, aln1$end, sep = '_'), sep = '_')
names(aln2A_seqs) <- paste(aln2A$chrom, paste(aln2A$start, aln2A$end, sep = '_'), sep = '_')
names(aln2B_seqs) <- paste(aln2B$chrom, paste(aln2B$start, aln2B$end, sep = '_'), sep = '_')
names(aln3A_seqs) <- paste(aln3A$chrom, paste(aln3A$start, aln3A$end, sep = '_'), sep = '_')
names(aln3B_seqs) <- paste(aln3B$chrom, paste(aln3B$start, aln3B$end, sep = '_'), sep = '_')




aln1_seqs[1,] <- reverseComplement(aln1_seqs[1,])
aln2B_seqs[1,] <- reverseComplement(aln2B_seqs[1,])
aln3B_seqs[1,] <- reverseComplement(aln3B_seqs[1,])





aln1_aln <- msa(aln1_seqs, method = 'Muscle')
aln2A_aln <- msa(aln2A_seqs, method = 'Muscle')
aln2B_aln <- msa(aln2B_seqs, method = 'Muscle')
aln3A_aln <- msa(aln3A_seqs, method = 'Muscle')
aln3B_aln <- msa(aln3B_seqs, method = 'Muscle')

aln1_aln
aln2A_aln
aln2B_aln # rev comp
aln3A_aln
aln3B_aln # rev comp


msaConvert()

unmasked(aln1_aln) %>% DNAMultipleAlignment() %>% detail()
unmasked(aln2A_aln)
unmasked(aln2B_aln)
unmasked(aln3A_aln)
unmasked(aln3B_aln)

PWaln1 <- pairwiseAlignment(aln1_seqs[1], aln1_seqs[2])
pid(PWaln1)

PWaln2A <- pairwiseAlignment(aln2A_seqs[1], aln2A_seqs[2])
pid(PWaln2A)
PWaln2B <- pairwiseAlignment(aln2B_seqs[1], aln2B_seqs[2])
pid(PWaln2B)

PWaln3A<- pairwiseAlignment(aln3A_seqs[1], aln3A_seqs[2])
pid(PWaln3A)
PWaln3B<- pairwiseAlignment(aln3B_seqs[1], aln3B_seqs[2])
pid(PWaln3B)



msa::msaConservationScore(aln1_aln)

get_islands <- function(island_info, genome){
  # using the island info dataframe, this function
  # extracts the sequence data for each island from the corresponding
  # fasta.
  locs <- island_info %>% select(seqid, Istart, Iend) %>%
    transmute(chrom=seqid, start=Istart, end=Iend)
  
  islands <- getSeq(genome, as(locs, "GRanges")) # this needs Biostrings and BSgenome loaded
  names(islands) <- island_info$island_ID
  
  return(islands)
}

