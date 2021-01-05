library(tidyverse)

gpa <- read_csv('./outputs/pan_genomes/ALL/gene_presence_absence.csv')

look <- gpa[grep('tetracycline ', gpa$Annotation),]



abricate <- read_tsv('./outputs/pan_genomes/ALL/pan_genome_abricate.tsv') %>% 
  mutate(locus_tags= SEQUENCE)

vfdb <- read_tsv('./outputs/pan_genomes/ALL/pan_genome_vfdb.tsv') %>% 
  mutate(locus_tags=SEQUENCE)


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




