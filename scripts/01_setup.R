library(tidyverse)
library(ggrepel)

#########


ani_tab <- read_tsv('./outputs/pyani_res/ANIm_percentage_identity.tab') %>% 
  pivot_longer(cols = -X1, names_to='genome', values_to='ANI') %>% 
  mutate(ani_dist=1-ANI) %>% select(-ANI)


ani_dist <- ani_tab %>% 
  pivot_wider(names_from = genome, values_from=ani_dist) %>% 
  column_to_rownames(var='X1') %>% as.dist()

ani_mds <- cmdscale(ani_dist) %>% as.data.frame() %>% 
  rownames_to_column(var='genome') 


ani_mds %>% ggplot(aes(x=V1, y=V2)) +
  geom_point()+
  geom_text_repel(aes(label=genome)) + 
  labs(x='MDS1', y='MDS2')


## TABLE VERIFYING DONOR RECIPIENT RESULT RELATIONSHIP


parental_genomes <- c('11601MD', '6631', '14229-5', '6461', '6067', '13150', 'JCC')


metadata <- ani_tab %>% 
  transmute(parental_genome = X1,
            result_genome = genome,
            ani_dist=ani_dist) %>% 
  filter(parental_genome %in% parental_genomes) %>% 
  filter(!(result_genome %in% parental_genomes)) %>% 
  group_by(result_genome) %>% 
  summarise(recipient_genome=parental_genome[which.min(ani_dist)]) %>% 
  mutate(
    donor_genome=case_when(
      result_genome == '11601MDx6631'        ~ '6631', 
      result_genome == '6461x13150'          ~ '6461', 
      result_genome == '6461x13150exp137'    ~ '6461', 
      result_genome == '6461x14229-5'        ~ '14229-5', 
      result_genome == '6461x6067'           ~ '6067', 
      grepl('Sample', result_genome)         ~ '14229-5', 
      result_genome == '13150x6631'          ~ '6631', 
      result_genome == '13150xJCC'           ~ 'JCC')) %>% 
  arrange(donor_genome, recipient_genome) %>% 
  mutate(pan_dir=paste(donor_genome,recipient_genome,result_genome,sep = '_'), 
         experiment=paste0('experiment_', seq_along(result_genome))) %>% 
  select(experiment, donor_genome, recipient_genome, result_genome, pan_dir)

metadata %>% write_tsv('./outputs/01_metadata.tsv')





