# setwd('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes/all_6461s')



cii <- read_csv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes/all_6461s/gifrop_out/clustered_island_info.csv')


gpa <- read_csv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes/all_6461s/gene_presence_absence.csv')


gpa <- gpa %>% arrange(`Genome Fragment`, `Order within Fragment`)


is.na(gpa$`6461`)

NOT6464 <- gpa[is.na(gpa$`6461`),]

loc_tag_class <- read_tsv('./outputs/result_locus_tag_classification.tsv')

gpa_long <-
  gpa %>%
  pivot_longer(cols = c(15:24), names_to='genome', values_to='locus_tags') %>% 
  left_join(loc_tag_class)

gpa_long %>% 
  filter(`Genome Fragment` == 1) %>% 
  filter(`Order within Fragment` > 435 & `Order within Fragment` < 525) %>% 
  ggplot(aes(x=`Order within Fragment`, y=genome, group=genome)) + 
  geom_point(aes(fill=classification), shape=22, size=3) + scale_fill_viridis_d()


gpa_long %>% 
  filter(`Genome Fragment` == 1) %>% 
  filter(`Order within Fragment` > 476 & `Order within Fragment` < 499) %>% 
  ggplot(aes(x=`Order within Fragment`, y=genome, group=genome)) + 
  geom_point(aes(fill=classification), shape=22, size=5) +
  scale_fill_viridis_d() + 
  geom_text(data=annot, aes(x=`Order within Fragment`, y=-6, label=Annotation), inherit.aes = FALSE, hjust='left', size=3) + 
  scale_y_discrete(expand=expansion(mult=c(.75,.03))) + coord_flip()
  
annot <- 
  gpa_long %>% 
  filter(`Genome Fragment` == 1) %>% 
  filter(`Order within Fragment` > 476 & `Order within Fragment` < 499) %>% 
  select(Gene, Annotation, `Order within Fragment`) %>% unique()





df <- data.frame(
  x = rep(c(2, 5, 7, 9, 12), 2),
  y = rep(c(1, 2), each = 5),
  z = factor(rep(1:5, each = 2)),
  w = rep(diff(c(0, 4, 6, 8, 10, 14)), 2)
)
ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = z), colour = "grey50")
ggplot(df, aes(x, y, width = w)) +
  geom_tile(aes(fill = z), colour = "grey50")
ggplot(df, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y, ymax = y + 1)) +
  geom_rect(aes(fill = z), colour = "grey50")






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
#   
gpa_test <- read_csv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gene_presence_absence.csv')%>%
  pivot_longer(cols = c(15:17), names_to='genome', values_to='locus_tag') %>%
  select(Gene, locus_tag, Annotation)
# 
# 
gff_test_all <- gff_test_all %>% left_join(gpa_test, by = 'locus_tag')
# 
# 
# 
gff_test_all %>% ggplot()+
  geom_rect(aes(xmin = is_start, xmax=is_end, ymin=yval, ymax=yval+1, fill=genome), color='black') +
  geom_text(aes(x=(is_start+is_end)/2, y=yval+.5, label=Gene), angle=90)
# 
# colnames(gff_test_all)
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



