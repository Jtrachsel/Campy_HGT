# setwd('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes/all_6461s')
library(tidyverse)
library(cowplot)

## LOCUS TAG ORDERS ###


# gff_pa


cii <- read_csv('./outputs/pan_genomes/6461s/gifrop_out/clustered_island_info.csv')


gpa <- read_csv('./outputs/pan_genomes/6461s/gene_presence_absence.csv',
                col_types = 'cccdddcdcdcdddccccccccccc')



gpa <- gpa %>% arrange(`Genome Fragment`, `Order within Fragment`)


is.na(gpa$`6461`)

NOT6464 <- gpa[is.na(gpa$`6461`),]

loc_tag_class <- read_tsv('./outputs/result_locus_tag_classification.tsv')


abricate <- read_tsv('./outputs/pan_genomes/6461s/pan_genome_abricate.tsv') %>% 
  mutate(locus_tags= SEQUENCE)

vfdb <- read_tsv('./outputs/pan_genomes/6461s/pan_genome_vfdb.tsv') %>% 
  mutate(locus_tags=SEQUENCE)


abricate <- rbind(abricate, vfdb)

gpa_abric <- 
  gpa %>%
  pivot_longer(cols = c(15:24), names_to='genome', values_to='locus_tags') %>%
  left_join(abricate) %>% 
  select(Gene, Annotation,`Genome Fragment`, PRODUCT, RESISTANCE) %>% 
  filter(!is.na(PRODUCT)) %>% 
  unique()


gpa_long %>% 
  filter(!is.na(RESISTANCE)) %>% 
  group_by(`Order within Fragment`, RESISTANCE) %>% 
  tally()


gpa_long <-
  gpa %>%
  pivot_longer(cols = c(15:25), names_to='genome', values_to='locus_tags') %>% 
  left_join(loc_tag_class) %>% 
  left_join(gpa_abric)

gpa_long_foc <- 
  gpa_long %>% 
  filter(`Genome Fragment` == 1) %>% 
  filter(`Order within Fragment` > 580 & `Order within Fragment` < 620) 


  


gpa_long_foc$RESISTANCE <- factor(gpa_long_foc$RESISTANCE)


gpa_long_foc$AMR <- ifelse(is.na(gpa_long_foc$locus_tags),NA,gpa_long_foc$PRODUCT)
# gpa_long_foc$AMR <- sub("aminoglycoside O-phosphotransferase ", "", gpa_long_foc$PRODUCT)

unique(gpa_long_foc$PRODUCT)

p1 <- 
  gpa_long_foc %>% 
  filter(genome != '6461') %>% 
  # filter(`Genome Fragment` == 1) %>% 
  # filter(`Order within Fragment` > 1581 & `Order within Fragment` < 1675) %>% 
  ggplot(aes(x=`Order within Fragment`, y=genome, group=genome)) + 
  geom_vline(xintercept = 590, color='orange')+
  geom_vline(xintercept = 610, color='orange')+
  geom_point(aes(fill=classification), shape=22, size=3, na.rm = TRUE) +
  geom_point(data=filter(gpa_long_foc, !is.na(PRODUCT) &(`Order within Fragment` > 580 & `Order within Fragment` < 620 & genome != '6461')),aes(color=AMR), na.rm = TRUE,  show.legend = FALSE)+
  scale_fill_viridis_d() + 
  scale_color_brewer(palette = 'Set1') + 
  theme_cowplot() + 
  theme(legend.position = 'top')
p1
# need to abricate and merge into annot by locus tag?
annot <- 
  gpa_long_foc %>% 
  # filter(`Genome Fragment` == 1) %>% 
  filter(`Order within Fragment` > 590 & `Order within Fragment` < 610) %>% 
  select(Gene, Annotation, `Order within Fragment`) %>% unique() %>% 
  left_join(gpa_abric) %>% 
  mutate(Annotation2 = ifelse(is.na(PRODUCT), Annotation, PRODUCT))


p2 <- 
  gpa_long_foc %>% 
  filter(`Genome Fragment` == 1 & genome != '6461') %>% 
  filter(`Order within Fragment` > 590 & `Order within Fragment` < 610) %>% 
  ggplot(aes(x=`Order within Fragment`, y=genome, group=genome)) + 
  geom_point(aes(fill=classification), shape=22, size=4.5, show.legend = FALSE) +
  geom_point(data=
               filter(gpa_long_foc, !is.na(AMR) &
                        (`Order within Fragment` > 590 &
                           `Order within Fragment` < 610) &
                        genome != '6461'),
             aes(color=AMR), size=2,na.rm = TRUE , shape=17, show.legend = FALSE)+
  geom_text(data=annot, aes(x=`Order within Fragment`, y=-5.5, label=Annotation2), inherit.aes = FALSE, hjust='left', size=3.5) + 
  scale_y_discrete(expand=expansion(mult=c(.75,.09))) + coord_flip()+
  scale_fill_viridis_d() +
  # scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = -45, vjust = -.9)) + 
  xlab('Order')

p2



fig_1 <- ggdraw()+
  draw_plot(p1, 0,.4,1,.6)+
  draw_plot(p2, 0,0,1,.4)+
  draw_plot_label(x=c(0,0), y=c(1,.45), label = c('A', 'B'))
fig_1


ggsave(fig_1,
       filename = './outputs/figure1.jpeg',
       width = 260,
       height = 200,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')



### Other islands ###




# 
# df <- data.frame(
#   x = rep(c(2, 5, 7, 9, 12), 2),
#   y = rep(c(1, 2), each = 5),
#   z = factor(rep(1:5, each = 2)),
#   w = rep(diff(c(0, 4, 6, 8, 10, 14)), 2)
# )
# ggplot(df, aes(x, y)) +
#   geom_tile(aes(fill = z), colour = "grey50")
# ggplot(df, aes(x, y, width = w)) +
#   geom_tile(aes(fill = z), colour = "grey50")
# ggplot(df, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y, ymax = y + 1)) +
#   geom_rect(aes(fill = z), colour = "grey50")






gff_test <- 
  read_tsv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gifrop_out/my_islands/SampleA_29_1.gff')
gff_test %>% ggplot() + geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax=1), color='black') +
  geom_text(aes(x=(start+end)/2, y=.5, label=product), angle=90)
# 
gff_test1 <- read_tsv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gifrop_out/my_islands/14229-5_17_1.gff')
# gff_test2 <- read_tsv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gifrop_out/my_islands/6461_173_1.gff')
gff_test3 <- read_tsv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gifrop_out/my_islands/SampleA_29_1.gff') %>%  mutate(seqid=as.character(seqid))
# 
gff_test_all <- bind_rows(gff_test1,
                          #gff_test2,
                          gff_test3) %>%
  mutate(genome=sub('(.*)_[0-9]+', '\\1', seqid),
         gene_name=sub('.*','',attributes)) %>%
  group_by(genome) %>%
  mutate(yval=cur_group_id(),
         is_start=start-min(start),
         len=end-start,
         is_end=is_start+len)
# 
# gff_test_all$attributes
# 
#   ### getting there...
# #   
# gpa_test <- read_csv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gene_presence_absence.csv')%>%
#   pivot_longer(cols = c(15:17), names_to='genome', values_to='locus_tag') %>%
#   select(Gene, locus_tag, Annotation)
# # 
# # 
# gff_test_all <- gff_test_all %>% left_join(gpa_test, by = 'locus_tag')
# # 
# # 
# # 
# gff_test_all %>% ggplot()+
#   geom_rect(aes(xmin = is_start, xmax=is_end, ymin=yval, ymax=yval+1, fill=genome), color='black') +
#   geom_text(aes(x=(is_start+is_end)/2, y=yval+.5, label=Gene), angle=90)
# # 
# # colnames(gff_test_all)
# #
# # gpa_test %>% pivot_longer(cols = c(15:17), names_to='genome', values_to='locus_tag')
# 
# 
# 
# '14229-5_21_1'
# '6461_173_1'
# 'SampleA_72_1'
# 
#



# For all locus tags, assign to be present in donor, present in recipient, present in both
# do this from 3 way pangenomes
# then make big pangenome from only recipients and results
# use this to make fig showing tiles of genes that are near the interesting site
# get a few genes up and down stream of any insertions

tst <- RESULTS[[7]]
pii <- read_csv('./outputs/pan_genomes/14229-5_6461_SampleG/gifrop_out/pan_with_island_info.csv')


pii %>% filter(all)


look <- gpa_long[2450:2470,]

which(gpa_long$locus_tags == 'BEGKIOMD_00456')



gpa_long_foc <- 
  gpa_long %>% 
  filter(`Genome Fragment` == 1) %>% 
  filter(`Order within Fragment` > 203 & `Order within Fragment` < 245)



gpa_long_foc$RESISTANCE <- factor(gpa_long_foc$RESISTANCE)

gpa_long_foc$AMR <- sub("aminoglycoside O-phosphotransferase ", "", gpa_long_foc$PRODUCT)

unique(gpa_long_foc$PRODUCT)

p1 <- 
  gpa_long_foc %>% 
  filter(genome != '6461') %>% 
  # filter(`Genome Fragment` == 1) %>% 
  # filter(`Order within Fragment` > 1581 & `Order within Fragment` < 1675) %>% 
  ggplot(aes(x=`Order within Fragment`, y=genome, group=genome)) + 
  # geom_vline(xintercept = 1604, color='orange')+
  # geom_vline(xintercept = 1619, color='orange')+
  geom_point(aes(fill=classification), shape=22, size=3, na.rm = TRUE) +
  geom_text(aes(label=`Accessory Order with Fragment`))+
  geom_point(data=filter(gpa_long_foc, !is.na(RESISTANCE) &(`Order within Fragment` > 1581 & `Order within Fragment` < 1675 & genome != '6461')),aes(color=AMR), na.rm = TRUE,  show.legend = FALSE)+
  scale_fill_viridis_d() + 
  scale_color_brewer(palette = 'Set1') + 
  theme_cowplot() + 
  theme(legend.position = 'top')
p1
# need to abricate and merge into annot by locus tag?
annot <- 
  gpa_long_foc %>% 
  # filter(`Genome Fragment` == 1) %>% 
  filter(`Order within Fragment` > 1604 & `Order within Fragment` < 1619) %>% 
  select(Gene, Annotation, `Order within Fragment`) %>% unique() %>% 
  left_join(gpa_abric) %>% 
  mutate(Annotation2 = ifelse(is.na(PRODUCT), Annotation, PRODUCT))


p2 <- 
  gpa_long_foc %>% 
  filter(`Genome Fragment` == 1 & genome != '6461') %>% 
  filter(`Order within Fragment` > 1604 & `Order within Fragment` < 1619) %>% 
  ggplot(aes(x=`Order within Fragment`, y=genome, group=genome)) + 
  geom_point(aes(fill=classification), shape=22, size=4.5, show.legend = FALSE) +
  geom_point(data=
               filter(gpa_long_foc, !is.na(AMR) &
                        (`Order within Fragment` > 1604 &
                           `Order within Fragment` < 1619) &
                        genome != '6461'),
             aes(color=AMR), size=2,na.rm = TRUE , shape=17, show.legend = FALSE)+
  geom_text(data=annot, aes(x=`Order within Fragment`, y=-5.5, label=Annotation2), inherit.aes = FALSE, hjust='left', size=3.5) + 
  scale_y_discrete(expand=expansion(mult=c(.75,.09))) + coord_flip()+
  scale_fill_viridis_d() +
  # scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = -45, vjust = -.9)) + 
  xlab('Order')

p2


########

gpa_long %>% filter(`Genome Fragment` == 1) %>% 
  group_by(`Order within Fragment`, classification) %>% tally() %>% 
  filter(classification == 'donor_result') %>% 
  filter(n > 1)
