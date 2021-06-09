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



### check in-vivo tet(O) location


# look <- read_tsv('./input_genomes/6461x13150exp137_AMR.tsv', col_types = c('ccddcccccccccdddddcccc'))


