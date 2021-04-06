library(tidyverse)
library(msa)
#### rpsL point mutation comparisons.


# 6631_11601MD_11601MDx6631
gpa <- read_csv('./outputs/pan_genomes/6631_11601MD_11601MDx6631/gene_presence_absence.csv')

grep('rpsL',gpa$Gene)

look <-  gpa[grep('rpsL',gpa$Gene),]

THESE <- 
  look %>%
  gather(c(15:17), key = 'genome', value='locus_tag') %>% 
  select(Gene, Annotation, `Genome Fragment`, `Order within Fragment`, genome, locus_tag) %>% 
  na.omit() %>% 
  # mutate(genome_frag=ifelse(`Genome Fragment` == 77, 'chrom', 'plasmid')) %>% 
  mutate(ID=paste(genome, locus_tag, sep = '_'))

names(THESE$locus_tag) <- THESE$ID

name_swap <- THESE$ID
names(name_swap) <- THESE$locus_tag


gens <- THESE$genome %>% unique() %>% paste(collapse = '|')

ffns <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*ffn')
ffns <- ffns[!grepl('proteins', ffns)]
ffns <- ffns[grepl(gens, ffns)]
ffns <- ffns[-3]

nucseqs <- lapply(ffns, readDNAStringSet)
nucseqs <- do.call(c, nucseqs)
names(nucseqs) <- sub('(.*_[0-9]+) .*','\\1',names(nucseqs))

rpsL_nuc <- nucseqs[THESE$locus_tag]

names(rpsL_nuc) <- name_swap[names(rpsL_nuc)]



rpsL_nuc %>% 
  msa::msaClustalOmega() %>% 
  msa::msaPrettyPrint(file = 'rpsL_nuc.pdf')


# from  6631 and 13150 
# 6631_13150_13150x6631


gpa <- read_csv('./outputs/pan_genomes/6631_13150_13150x6631/gene_presence_absence.csv')

grep('rpsL',gpa$Gene)

look <-  gpa[grep('rpsL',gpa$Gene),]

# FAOCEDLB_00448
# BACMIICJ_00447
# KINOJHBP_00456


THESE <- 
  look %>%
  gather(c(15:17), key = 'genome', value='locus_tag') %>% 
  select(Gene, Annotation, `Genome Fragment`, `Order within Fragment`, genome, locus_tag) %>% 
  na.omit() %>% 
  # mutate(genome_frag=ifelse(`Genome Fragment` == 77, 'chrom', 'plasmid')) %>% 
  mutate(ID=paste(genome, locus_tag, sep = '_'))

names(THESE$locus_tag) <- THESE$ID

name_swap <- THESE$ID
names(name_swap) <- THESE$locus_tag


gens <- THESE$genome %>% unique() %>% paste(collapse = '|')

ffns <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*ffn')
ffns <- ffns[!grepl('proteins', ffns)]
ffns <- ffns[grepl(gens, ffns)]
ffns <- ffns[-c(1,4,5,6)]

nucseqs <- lapply(ffns, readDNAStringSet)
nucseqs <- do.call(c, nucseqs)
names(nucseqs) <- sub('(.*_[0-9]+) .*','\\1',names(nucseqs))

rpsL_nuc <- nucseqs[THESE$locus_tag]

names(rpsL_nuc) <- name_swap[names(rpsL_nuc)]


rpsL_nuc %>% 
  msa::msaClustalOmega() %>% 
  msa::msaPrettyPrint(file = '6631_13150_rpsL_nuc.pdf')


# 6631 / 11601MD / 6631x11601MD_NES

gpa <- read_csv('./outputs/pan_genomes/6631_11601MD_11601MDx6631/gene_presence_absence.csv')

grep('rpsL',gpa$Gene)

# NDOAKFMB_00470
# CKJNOHCG_00470
# KINOJHBP_00456

look <-  gpa[grep('rpsL',gpa$Gene),]

THESE <- 
  look %>%
  gather(c(15:17), key = 'genome', value='locus_tag') %>% 
  select(Gene, Annotation, `Genome Fragment`, `Order within Fragment`, genome, locus_tag) %>% 
  na.omit() %>% 
  # mutate(genome_frag=ifelse(`Genome Fragment` == 77, 'chrom', 'plasmid')) %>% 
  mutate(ID=paste(genome, locus_tag, sep = '_'))

names(THESE$locus_tag) <- THESE$ID

name_swap <- THESE$ID
names(name_swap) <- THESE$locus_tag


gens <- THESE$genome %>% unique() %>% paste(collapse = '|')


ffns <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*ffn')
ffns <- ffns[!grepl('proteins', ffns)]
ffns <- ffns[grepl(gens, ffns)]
ffns <- ffns[-c(3)]

nucseqs <- lapply(ffns, readDNAStringSet)
nucseqs <- do.call(c, nucseqs)
names(nucseqs) <- sub('(.*_[0-9]+) .*','\\1',names(nucseqs))

rpsL_nuc <- nucseqs[THESE$locus_tag]

names(rpsL_nuc) <- name_swap[names(rpsL_nuc)]



rpsL_nuc %>% 
  msa::msaClustalOmega() %>% 
  msa::msaPrettyPrint(file = '6631_11601MD_11601MDx6631_rpsL_nuc.pdf')





