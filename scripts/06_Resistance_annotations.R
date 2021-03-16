library(tidyverse)

gpa <- read_csv('./outputs/pan_genomes/ALL/gene_presence_absence.csv')

look <- gpa[grep('tetracycline ', gpa$Annotation),]



abricate <- read_tsv('./outputs/pan_genomes/ALL/pan_genome_abricate.ncbi') %>% 
  mutate(Gene= SEQUENCE)

vfdb <- read_tsv('./outputs/pan_genomes/ALL/pan_genome_abricate.vfdb') %>% 
  mutate(Gene=SEQUENCE)


abricate <- rbind(abricate, vfdb)

gpa_abric <- 
  gpa %>%
  pivot_longer(cols = c(15:24), names_to='genome', values_to='locus_tags') %>%
  left_join(abricate) %>% 
  select(Gene, Annotation,`Genome Fragment`, PRODUCT, RESISTANCE) %>% 
  filter(!is.na(PRODUCT)) %>% 
  unique()



gpa_abric %>% 
  select(-Gene, -`Genome Fragment`) %>%
  filter(!is.na(RESISTANCE)) %>% 
  transmute(prokka_annotation=Annotation, 
            abricate_annotation=PRODUCT, 
            RESISTANCE=RESISTANCE) %>%
  unique() %>% 
  write_tsv('./outputs/AMR_annotation_differences.tsv')



#### rpsL point mutation comparisons.

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


# faas <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*faa')
# faas <- faas[!grepl('proteins', faas)]
# faas <- faas[grepl(gens, faas)]


ffns <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*ffn')
ffns <- ffns[!grepl('proteins', ffns)]
ffns <- ffns[grepl(gens, ffns)]
ffns <- ffns[-3]

nucseqs <- lapply(ffns, readDNAStringSet)
nucseqs <- do.call(c, nucseqs)
names(nucseqs) <- sub('(.*_[0-9]+) .*','\\1',names(nucseqs))

rpsL_nuc <- nucseqs[THESE$locus_tag]

names(rpsL_nuc) <- name_swap[names(rpsL_nuc)]
library(msa)



rpsL_nuc %>% 
  msa::msaClustalOmega() %>% 
  msa::msaPrettyPrint(file = 'rpsL_nuc.pdf')


# from  6631 and 13150 



gpa <- read_csv('./outputs/pan_genomes/6631_13150_13150x6631/gene_presence_absence.csv')

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


# faas <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*faa')
# faas <- faas[!grepl('proteins', faas)]
# faas <- faas[grepl(gens, faas)]


ffns <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*ffn')
ffns <- ffns[!grepl('proteins', ffns)]
ffns <- ffns[grepl(gens, ffns)]
ffns <- ffns[-c(1,4,5,6)]

nucseqs <- lapply(ffns, readDNAStringSet)
nucseqs <- do.call(c, nucseqs)
names(nucseqs) <- sub('(.*_[0-9]+) .*','\\1',names(nucseqs))

rpsL_nuc <- nucseqs[THESE$locus_tag]

names(rpsL_nuc) <- name_swap[names(rpsL_nuc)]
library(msa)



rpsL_nuc %>% 
  msa::msaClustalOmega() %>% 
  msa::msaPrettyPrint(file = '6631_13150_rpsL_nuc.pdf')

