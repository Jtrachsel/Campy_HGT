# setwd('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes/all_6461s')
library(tidyverse)
library(cowplot)

## LOCUS TAG ORDERS ###


# gff_pa


cii <- read_csv('./outputs/pan_genomes/6461s/gifrop_out/clustered_island_info.csv')


gpa <- read_csv('./outputs/pan_genomes/6461s/gene_presence_absence.csv',
                col_types = 'cccdddcdcdcdddcccccccccc')



gpa <- gpa %>% arrange(`Genome Fragment`, `Order within Fragment`)


# is.na(gpa$`6461`)

# NOT6464 <- gpa[is.na(gpa$`6461`),]

loc_tag_class <- read_tsv('./outputs/result_locus_tag_classification.tsv')


abricate <- read_tsv('./outputs/pan_genomes/6461s/pan_genome_abricate.ncbi') %>% 
  mutate(Gene= SEQUENCE)

vfdb <- read_tsv('./outputs/pan_genomes/6461s/pan_genome_abricate.vfdb') %>% 
  mutate(Gene=SEQUENCE)


abricate <- rbind(abricate, vfdb)

gpa_abric <- 
  gpa %>%
  pivot_longer(cols = c(15:24), names_to='genome', values_to='locus_tags') %>%
  left_join(abricate) %>% 
  select(Gene, Annotation,`Genome Fragment`, PRODUCT, RESISTANCE) %>% 
  filter(!is.na(PRODUCT)) %>% 
  unique()


insert_order <- function(dat){
  dat %>% 
    filter(!is.na(locus_tags)) %>%
    arrange(locus_tags) %>% 
    mutate(ORDER=seq_along(locus_tags))
}

order_rel_res <- function(dat){
  RES_ORDER <- 
    dat %>%
    filter(!is.na(RESISTANCE)) %>% 
    group_by(ORDER, RESISTANCE) %>% 
    tally() %>% 
    filter(grepl('AMIKACIN;GENTAMICIN;KANAMYCIN;TOBRAMYCIN', RESISTANCE)) %>% 
    pull(ORDER)
  
  dat %>% mutate(O_rel_res=ORDER-RES_ORDER)
  
}

gpa_long <-
  gpa %>%
  pivot_longer(cols = c(15:24), names_to='genome', values_to='locus_tags') %>% 
  left_join(loc_tag_class) %>% 
  left_join(gpa_abric) %>% 
  # filter(genome != '6461') %>% 
  group_by(genome) %>%
  nest() %>% 
  mutate(data2=map(.x = data, .f = insert_order)) %>% 
  select(genome, data2) %>% 
  mutate(data3=map(.x = data2, .f = order_rel_res)) %>% 
  select(genome, data3) %>% 
  unnest(data3) %>% 
  ungroup()

loook <- gpa_long %>% 
  group_by(Gene) %>% 
  summarise(AVORD=mean(ORDER)) %>% 
  arrange(AVORD) %>% 
  mutate(rankORD=rank(AVORD, ties.method = 'first')) %>% 
  select(Gene, rankORD)

gpa_long <- 
  gpa_long %>% left_join(loook)


resist_pan_coords <-
  gpa_long %>%
  # filter(!is.na(RESISTANCE)) %>%
  filter(grepl('AMIKACIN;GENTAMICIN;KANAMYCIN;TOBRAMYCIN', RESISTANCE)) %>%
  group_by(rankORD, RESISTANCE) %>%
  tally() %>%
  pull(rankORD)

pan_low <- signif(resist_pan_coords - 40, 3)
pan_high <- signif(resist_pan_coords + 30, 3)

zoom1_low <- pan_low + 30
zoom1_high <- pan_high - 20

# zoom2_low <- pan_low + 65
# zoom2_high <- pan_high - 65

# 
# gpa_long_foc <- 
#   gpa_long %>% 
#   filter(`Genome Fragment` == 1) %>% 
#   filter(`Order within Fragment` > pan_low & `Order within Fragment` < pan_high) 


gpa_long_foc <- 
  gpa_long %>% 
  filter(`Genome Fragment` == 1) %>% 
  filter(rankORD > pan_low & rankORD < pan_high) 


  


gpa_long_foc$RESISTANCE <- factor(gpa_long_foc$RESISTANCE)
# 
# 

# gpa_long_foc[grepl('AMR', gpa_long_foc$Annotation),]
gpa_long_foc$AMR <- ifelse(!is.na(gpa_long_foc$locus_tags) & !is.na(gpa_long_foc$RESISTANCE),gpa_long_foc$PRODUCT, NA)
# gpa_long_foc$AMR <- sub("aminoglycoside O-phosphotransferase ", "", gpa_long_foc$AMR)
unique(gpa_long_foc$AMR )

gpa_long_foc$vir <- ifelse(!is.na(gpa_long_foc$locus_tags) & !is.na(gpa_long_foc$PRODUCT)& is.na(gpa_long_foc$RESISTANCE) ,gpa_long_foc$PRODUCT,NA)

gpa_long_foc$vir <- sub('(\\(.*\\) [A-Za-z0-9\\/ \\-]+) \\[.*\\]','\\1',gpa_long_foc$vir)

unique(gpa_long_foc$vir)

# hclust order y axis by similarity #
check <- 
  gpa_long_foc %>% 
  mutate(present=ifelse(is.na(locus_tags), 0,1)) %>% 
  select(genome, present, Gene) %>% 
  spread(key=Gene, value=present, fill=0) %>% 
  column_to_rownames('genome') %>% 
  as.matrix() %>% 
  dist() %>% 
  hclust()

lab_orders <- check$labels[check$order]

gpa_long_foc <- 
  gpa_long_foc %>% 
  mutate(genome = factor(genome, levels =lab_orders))

LOOOOK <- gpa_long_foc %>% group_by(rankORD, Gene) %>% tally()

# look <- gpa_long_foc %>% spread(key=genome, value=locus_tags)

p1_data <- 
  gpa_long_foc %>% 
  filter(!is.na(locus_tags)) %>% 
  filter(genome != '6461') 
p1 <- 
  p1_data %>% 
  ggplot(aes(x=rankORD, y=genome, group=genome)) + 
  geom_vline(xintercept = zoom1_low, color='orange', size=2)+
  geom_vline(xintercept = zoom1_high, color='orange', size=2)+
  geom_point(aes(fill=classification), shape=22, size=4, na.rm = TRUE) +
  # geom_point(data=filter(gpa_long_foc, !is.na(PRODUCT) &(`Order within Fragment` > 580 & `Order within Fragment` < 620 & genome != '6461')),aes(color=AMR), na.rm = TRUE,  show.legend = FALSE)+
  geom_point(data=filter(p1_data, !is.na(AMR)),aes(size="APH(2'')-If"), color='black', fill='red', shape=21)+
  geom_point(data=filter(p1_data, !is.na(vir)),aes(alpha='virulence related'), color='black',fill='pink', shape=21)+
  scale_fill_viridis_d() + 
  scale_color_brewer(palette = 'Set1') + 
  scale_alpha_manual(name = "", values = 1) +
  scale_size_manual(name = "", values = 2) +
  theme_cowplot() + 
  theme(legend.position = 'top', 
        axis.title.x=element_text(size=11))+
  guides(fill=guide_legend(nrow = 2)) + 
  xlab('Order')
p1



# need to abricate and merge into annot by locus tag?
annot <- 
  gpa_long_foc %>% 
  # filter(`Genome Fragment` == 1) %>% 
  filter(rankORD >= zoom1_low & rankORD <= zoom1_high) %>% 
  select(Gene, Annotation, rankORD) %>% unique() %>% 
  left_join(gpa_abric) %>% 
  mutate(Annotation2 = ifelse(is.na(PRODUCT), Annotation, PRODUCT), 
         Annotation2=sub('\\[.*\\]','',Annotation2), 
         Annotation3=sub('\\((.*)\\).*','\\1',Annotation2), 
         Annotation4=sub('Putative ','',Annotation3), 
         Annotation5=sub(' protein','',Annotation4), 
         Annotation6=sub('hypothetical', 'hypothetical protein', Annotation5))




p2 <- 
  gpa_long_foc %>% 
  # filter(`Genome Fragment` == 1 & genome != '6461') %>% 
  filter(rankORD >= zoom1_low & rankORD <= zoom1_high) %>% 
  filter(!is.na(locus_tags)) %>% 
  ggplot(aes(x=rankORD, y=genome, group=genome)) + 
  geom_vline(xintercept = zoom1_low -1, color='orange', size=2)+
  geom_vline(xintercept = zoom1_high +1, color='orange', size=2)+
  geom_point(aes(fill=classification), shape=22, size=4.5, show.legend = FALSE) +
  geom_point(data=filter(gpa_long_foc, !is.na(AMR) & genome != '6461'),
             fill='red',color='black', size=2,na.rm = TRUE , shape=21, show.legend = FALSE)+
  geom_point(data=filter(gpa_long_foc, !is.na(vir) & genome != '6461' &
                           rankORD >=zoom1_low &
                           rankORD <= zoom1_high),
             size=2,na.rm = TRUE , shape=21, color='black',fill='pink', show.legend = FALSE)+
  geom_text(data=annot, aes(x=rankORD, y=-5.5, label=Annotation6), inherit.aes = FALSE, hjust='left', size=3.5) + 
  scale_y_discrete(expand=expansion(mult=c(.75,.09))) + coord_flip()+
  scale_fill_viridis_d() +
  # scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1') + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = -45, vjust = -.9, size=10), 
        axis.title.x = element_blank()) + 
  xlab('Order')

p2



fig_1 <- ggdraw()+
  draw_plot(p1, 0,.6,1,.4)+
  draw_plot(p2, 0,0,1,.6)+
  draw_plot_label(x=c(0,0), y=c(1,.6), label = c('A', 'B'))
fig_1


ggsave(fig_1,
       filename = './outputs/6461s.jpeg',
       width = 260,
       height = 210,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')

# lookshere <- ggplot_build(p1)
# 
# 
# unique(lookshere$data[[3]]$fill)

# Don_rec_res "#35B779FF"
# Don_res "#440154FF"
# Rec_res "#31688EFF"
# res_only "#FDE725FF"
### Two other transfer events islands ###

# Idea, 
# what functions or cell locations are enriched in:
# 1) Transferred genes,
# 2) Conserved genes adjacent to insertion
# 3) 






######## FUNCTION ######

# Next figure #

generate_transfer_plots <- function(base_dir,
                                    DONOR,
                                    RECIPIENT,
                                    RESULT,
                                    TARGET,
                                    POINT_MUT=F,
                                    adj1_low,
                                    adj1_high,
                                    low_shrink,
                                    high_shrink, 
                                    genomic_fragment=1, 
                                    TARGET_COL=Gene){
  # browser()
  # outputs/pan_genomes/6631_11601MD_11601MDx6631/ path to directory that contains pangenome calc by roary
  rtab_path <- paste0(base_dir, '/gene_presence_absence.Rtab')
  loc_tag_class <- read_tsv(rtab_path) %>% 
    column_to_rownames(var = 'Gene')
  # THIS NEEDS TO BE CHANGED, SPLIT BASE_DIR FOR PICKING WHICH COLUMN IS WHICH GENOME
  donor_recipient_result <- rownames(loc_tag_class[rowSums(loc_tag_class) == 3,]) # donor recipient result
  donor_recipient <- rownames(loc_tag_class[loc_tag_class[,DONOR] > 0 & loc_tag_class[,RECIPIENT] > 0 & loc_tag_class[,RESULT] == 0, ])# donor_recipient
  recipient_result <- rownames(loc_tag_class[loc_tag_class[,DONOR] == 0 & loc_tag_class[,RECIPIENT] > 0 & loc_tag_class[,RESULT] > 0, ])  # recipient_result 
  donor_result <- rownames(loc_tag_class[loc_tag_class[,DONOR] > 0 & loc_tag_class[,RECIPIENT] == 0 & loc_tag_class[,RESULT] > 0, ]) # donor_result
  
  donor <- rownames(loc_tag_class[loc_tag_class[,DONOR] > 0 & loc_tag_class[,RECIPIENT] == 0 & loc_tag_class[,RESULT] == 0, ])# donor only
  recipient <- rownames(loc_tag_class[loc_tag_class[,DONOR] == 0 & loc_tag_class[,RECIPIENT] > 0 & loc_tag_class[,RESULT] == 0, ])# recipient only
  result <- rownames(loc_tag_class[loc_tag_class[,DONOR] == 0 & loc_tag_class[,RECIPIENT] == 0 & loc_tag_class[,RESULT] > 0,] )# result only
  
  gpa_path <- paste0(base_dir,'/gene_presence_absence.csv')
  gpa <- read_csv(gpa_path,
                  col_types = 'cccdddcdcdcdddccc') %>%
    mutate(loc_tag_class=
             case_when(
               Gene %in% donor_recipient_result ~ 'donor_recipient_result', 
               Gene %in% donor_recipient        ~ 'donor_recipient', 
               Gene %in% recipient_result       ~ 'recipient_result', 
               Gene %in% donor_result           ~ 'donor_result', 
               Gene %in% donor                  ~ 'donor', 
               Gene %in% recipient              ~ 'recipient', 
               Gene %in% result                 ~ 'result'
             ))
  
  
  
  gpa <- gpa %>% arrange(`Genome Fragment`, `Order within Fragment`)
  
  
  # is.na(gpa$`6461`)
  
  # NOT6464 <- gpa[is.na(gpa$`6461`),]
  
  # loc_tag_class <- read_tsv('./outputs/result_locus_tag_classification.tsv')
  
  AMR_path <- paste0(base_dir, '/pan_genome_abricate.ncbi')
  abricate <- read_tsv(AMR_path) %>% 
    mutate(Gene= SEQUENCE)
  
  vfdb_path <- paste0(base_dir, '/pan_genome_abricate.vfdb')
  vfdb <- read_tsv(vfdb_path) %>% 
    mutate(Gene=SEQUENCE)
  
  
  abricate <- rbind(abricate, vfdb)
  
  gpa_abric <- 
    gpa %>%
    pivot_longer(cols = c(15:17), names_to='genome', values_to='locus_tags') %>%
    left_join(abricate) %>% 
    select(Gene, Annotation,`Genome Fragment`, PRODUCT, RESISTANCE) %>% 
    filter(!is.na(PRODUCT)) %>% 
    unique()
  
  
  
  gpa_long <-
    gpa %>%
    pivot_longer(cols = c(15:17), names_to='genome', values_to='locus_tags') %>% 
    # left_join(loc_tag_class) %>% 
    left_join(gpa_abric) 
    
  res_order <- gpa_long %>% filter(genome == all_of(RESULT)) %>%
    arrange(locus_tags) %>% 
    mutate(order_in_result_genome=seq_along(locus_tags)) %>% 
    select(Gene, order_in_result_genome)
  
  gpa_long <- gpa_long %>% left_join(res_order)
    
  # pick resist? return one result for each resist?
  # browser()
  if(POINT_MUT){
    resist_pan_coords <- 
      gpa_long %>% 
      filter(`Genome Fragment` == genomic_fragment) %>% 
      # filter(!is.na(RESISTANCE)) %>% 
      group_by(order_in_result_genome,`Genome Fragment`, {{TARGET_COL}}) %>% 
      tally() %>% 
      filter(grepl(TARGET, {{TARGET_COL}}, fixed = T)) %>% 
      pull(order_in_result_genome)
    
    pan_low <- signif(resist_pan_coords - adj1_low, 3)
    pan_high <- signif(resist_pan_coords + adj1_high, 3)
    
    
  }else{
    
    resist_pan_coords <- 
      gpa_long %>% 
      filter(`Genome Fragment` ==genomic_fragment) %>% 
      filter(!is.na(RESISTANCE)) %>% 
      group_by(order_in_result_genome,`Genome Fragment`, RESISTANCE) %>% 
      tally() %>% 
      filter(grepl(TARGET, RESISTANCE)) %>% 
      pull(order_in_result_genome)
    
    pan_low <- signif(resist_pan_coords - adj1_low, 3)
    pan_high <- signif(resist_pan_coords + adj1_high, 3)
    
  }
  
  
  # 
  # pan_low <- signif(resist_pan_coords - 70, 3)
  # pan_high <- signif(resist_pan_coords + 100, 3)
  # 
  zoom1_low <- pan_low + low_shrink
  zoom1_high <- pan_high - high_shrink
  
  # zoom2_low <- pan_low + 65
  # zoom2_high <- pan_high - 65
  
  
  gpa_long_foc <- 
    gpa_long %>% 
    filter(`Genome Fragment` == genomic_fragment) %>% 
    filter(order_in_result_genome > pan_low & order_in_result_genome < pan_high) 
  
  
  
  
  
  gpa_long_foc$RESISTANCE <- factor(gpa_long_foc$RESISTANCE)
  # 
  # 
  
  # gpa_long_foc[grepl('AMR', gpa_long_foc$Annotation),]
  gpa_long_foc$AMR <- ifelse(!is.na(gpa_long_foc$locus_tags) & !is.na(gpa_long_foc$RESISTANCE),gpa_long_foc$PRODUCT, NA)
  # gpa_long_foc$AMR <- sub("aminoglycoside O-phosphotransferase ", "", gpa_long_foc$AMR)
  unique(gpa_long_foc$AMR )
  
  gpa_long_foc$vir <- ifelse(!is.na(gpa_long_foc$locus_tags) & !is.na(gpa_long_foc$PRODUCT)& is.na(gpa_long_foc$RESISTANCE) ,gpa_long_foc$PRODUCT,NA)
  
  gpa_long_foc$vir <- sub('(\\(.*\\) [A-Za-z0-9\\/ \\-]+) \\[.*\\]','\\1',gpa_long_foc$vir)
  
  unique(gpa_long_foc$vir)
  
  # hclust order y axis by similarity #
  check <- 
    gpa_long_foc %>% 
    mutate(present=ifelse(is.na(locus_tags), 0,1)) %>% 
    select(genome, present, Gene) %>% 
    spread(key=Gene, value=present, fill=0) %>% 
    column_to_rownames('genome') %>% 
    as.matrix() %>% 
    dist() %>% 
    hclust()
  
  lab_orders <- check$labels[check$order]
  
  
  if(POINT_MUT){
    
    gpa_long_foc <- 
      gpa_long_foc %>% 
      filter(!is.na(locus_tags)) %>% 
      filter(genome == RESULT) %>% 
      mutate(genome = factor(genome, levels =lab_orders)) %>%
      mutate(AMR=ifelse(POINT_MUT & {{TARGET_COL}} == TARGET, TRUE, NA)) %>% 
      arrange((locus_tags))
    
  }else{
    gpa_long_foc <- 
      gpa_long_foc %>% 
      filter(!is.na(locus_tags)) %>% 
      filter(genome == RESULT) %>% 
      mutate(genome = factor(genome, levels =lab_orders)) %>%
      arrange((locus_tags))
    
  }
  
    
  
  
  
  
  p1_data <- 
    gpa_long_foc %>% 
    filter(!is.na(locus_tags)) %>% 
    filter(genome == RESULT) 
  # browser()
  p1 <- 
    p1_data %>%
    ggplot(aes(x=order_in_result_genome, y=genome, group=genome)) + 
    geom_vline(xintercept = zoom1_low, color='orange', size=2)+
    geom_vline(xintercept = zoom1_high, color='orange', size=2)+
    geom_point(aes(fill=loc_tag_class), shape=22, size=4, na.rm = TRUE) +
    # geom_point(data=filter(gpa_long_foc, !is.na(PRODUCT) &(`Order within Fragment` > 580 & `Order within Fragment` < 620 & genome != '6461')),aes(color=AMR), na.rm = TRUE,  show.legend = FALSE)+
    geom_point(data=filter(p1_data, !is.na(AMR)),  color='red', na.rm = TRUE)+
    geom_point(data=filter(p1_data, !is.na(vir)),  color='pink', na.rm = TRUE)+
    # scale_fill_viridis_d() + 
    scale_color_brewer(palette = 'Set1') + 
    theme_cowplot() + 
    theme(legend.position = 'top', 
          axis.title.x=element_text(size=11)) + 
    scale_fill_manual(values = c(recipient_result = "#35B779FF", 
                                 donor_recipient_result = "#440154FF", 
                                 donor_result = "#31688EFF", 
                                 result_only="#FDE725FF"))
  p1
  
  # Don_rec_res "#35B779FF"
  # Don_res "#440154FF"
  # Rec_res "#31688EFF"
  # res_only "#FDE725FF"
  
  # need to abricate and merge into annot by locus tag?
  
  annot <- 
    gpa_long_foc %>% 
    # filter(`Genome Fragment` == 1) %>% 
    filter(order_in_result_genome >= zoom1_low & order_in_result_genome <= zoom1_high) %>% 
    select(Gene, Annotation, order_in_result_genome) %>% unique() %>% 
    left_join(gpa_abric) %>% 
    mutate(Annotation2 = ifelse(is.na(PRODUCT), Annotation, PRODUCT), 
           Annotation2=sub('\\[.*\\]','',Annotation2), 
           Annotation3=sub('\\((.*)\\).*','\\1',Annotation2))
  
  
  p2_data <- 
    gpa_long_foc %>% 
    filter(genome == RESULT) %>%
    filter(order_in_result_genome >= zoom1_low & order_in_result_genome <= zoom1_high) %>% 
    filter(!is.na(locus_tags))
  
  
  p2 <- 
    p2_data %>% 
    ggplot(aes(x=order_in_result_genome, y=genome, group=genome)) + 
    geom_vline(xintercept = zoom1_low -1, color='orange', size=2)+
    geom_vline(xintercept = zoom1_high +1, color='orange', size=2)+
    geom_point(aes(fill=loc_tag_class), shape=22, size=4.5, show.legend = FALSE) +
    geom_point(data=filter(p2_data, !is.na(AMR)),fill='red', color='black', size=2,na.rm = TRUE , shape=21, show.legend = FALSE)+
    geom_point(data=filter(p2_data, !is.na(vir)),fill='pink', color='black', size=2,na.rm = TRUE , shape=21, show.legend = FALSE)+
    geom_text(data=annot, aes(x=order_in_result_genome, y=0, label=Annotation3), inherit.aes = FALSE, hjust='left', size=3.5) + 
    scale_y_discrete(expand=expansion(mult=c(.55,.9))) + coord_flip()+
    # scale_fill_viridis_d() +
    # scale_fill_brewer(palette = 'Set1')+
    scale_color_brewer(palette = 'Set1') + 
    theme_cowplot() + 
    theme(axis.text.x = element_text(angle = -45, vjust = -.7, size=10), 
          axis.title.x = element_blank()) + 
    xlab('Order')+
    scale_fill_manual(values = c(recipient_result = "#35B779FF", 
                                 donor_recipient_result = "#440154FF", 
                                 donor_result = "#31688EFF", 
                                 result_only="#FDE725FF"))
  
  
  
  
  ## return
  print(p1)
  
  print(p2)
  
  
  return(list(p1,p2))
  
}




x6461x13150 <- generate_transfer_plots(base_dir = 'outputs/pan_genomes/6461_13150_6461x13150', 
                        DONOR='6461', 
                        RECIPIENT = '13150', 
                        RESULT='6461x13150', 
                        adj1_low = 50, 
                        adj1_high = 50, 
                        low_shrink = 25,
                        high_shrink = 25, 
                        TARGET = 'TETRACYCLINE')



x6461x6067 <- generate_transfer_plots(base_dir ='outputs/pan_genomes/6067_6461_6461x6067/', 
                                DONOR = '6067', 
                                RECIPIENT=, '6461',
                                RESULT = '6461x6067',
                                TARGET = 'BETA-LACTAM', 
                                adj1_low = 50,
                                adj1_high= 50,
                                low_shrink=25, 
                                high_shrink=25)
####


x13150x6631 <- generate_transfer_plots(base_dir = 'outputs/pan_genomes/6631_13150_13150x6631',
                        DONOR = '6631',
                        RECIPIENT = '13150',
                        RESULT = '13150x6631', 
                        TARGET = 'rpsL', 
                        POINT_MUT = T,
                        adj1_low = 30,
                        adj1_high= 30,
                        low_shrink=15, 
                        high_shrink=15)



###
x6461x13150exp137 <- generate_transfer_plots(base_dir = 'outputs/pan_genomes/6461_13150_6461x13150exp137/',
                                       DONOR = '6461',
                                       RECIPIENT = '13150',
                                       RESULT = '6461x13150exp137', 
                                       TARGET = 'TETRACYCLINE', 
                                       POINT_MUT = T,
                                       adj1_low = 30,
                                       adj1_high= 30,
                                       low_shrink=15, 
                                       high_shrink=15,
                                       genomic_fragment = 1, 
                                       TARGET_COL = RESISTANCE)



x6461x13150exp137[[1]]

x6461x13150exp137[[2]]




x6461x13150_vivo <- ggdraw()+
  draw_plot(x6461x13150exp137[[1]], 0,.8,1,.2)+
  draw_plot(x6461x13150exp137[[2]], 0,0,1,.8)+
  draw_plot_label(x=c(0,0), y=c(1,.8), label = c('A', 'B'))
x6461x13150_vivo


ggsave(x6461x13150_vivo,
       filename = './outputs/x6461x13150_vivo.jpeg',
       width = 260,
       height = 290,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')




###



# 6461x13150



x6461x13150__2 <- generate_transfer_plots(base_dir = 'outputs/pan_genomes/6461_13150_6461x13150/',
                                             DONOR = '6461',
                                             RECIPIENT = '13150',
                                             RESULT = '6461x13150', 
                                             TARGET = 'TETRACYCLINE', 
                                             POINT_MUT = T,
                                             adj1_low = 30,
                                             adj1_high= 30,
                                             low_shrink=15, 
                                             high_shrink=15,
                                             genomic_fragment = 2, 
                                             TARGET_COL = RESISTANCE)




x6461x13150_plasmid <- ggdraw()+
  draw_plot(x6461x13150__2[[1]], 0,.8,1,.2)+
  draw_plot(x6461x13150__2[[2]], 0,0,1,.8)+
  draw_plot_label(x=c(0,0), y=c(1,.8), label = c('A', 'B'))
x6461x13150_plasmid


ggsave(x6461x13150_plasmid,
       filename = './outputs/x6461x13150_plasmid_tetO.jpeg',
       width = 260,
       height = 290,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')


x6461x13150__2[[1]]
x6461x13150__2[[2]]

#

####
LOOK <- x13150x6631[[2]]$data
LOOK$loc_tag_class
LOOK2 <- LOOK %>% filter(loc_tag_class != 'recipient_result')
gff <- gff_parse3('./outputs/pan_genomes/6631_13150_13150x6631/13150x6631.gff')

gff$locus_tag %in% LOOK2$locus_tags

res <- LOOK2 %>% left_join(gff, c('locus_tags'= 'locus_tag'))



fig_2 <- ggdraw()+
  draw_plot(x6461x13150[[1]], 0,.8,1,.2)+
  draw_plot(x6461x13150[[2]], 0,0,1,.8)+
  draw_plot_label(x=c(0,0), y=c(1,.8), label = c('A', 'B'))
fig_2


ggsave(fig_2,
       filename = './outputs/x6461x13150.jpeg',
       width = 260,
       height = 290,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')





fig_3 <- ggdraw()+
  draw_plot(x6461x6067[[1]], 0,.8,1,.2)+
  draw_plot(x6461x6067[[2]], 0,0,1,.8)+
  draw_plot_label(x=c(0,0), y=c(1,.8), label = c('A', 'B'))
fig_3


ggsave(fig_3,
       filename = './outputs/x6461x6067.jpeg',
       width = 260,
       height = 290,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')



fig_4 <- ggdraw()+
  draw_plot(x13150x6631[[1]], 0,.8,1,.2)+
  draw_plot(x13150x6631[[2]], 0,0,1,.8)+
  draw_plot_label(x=c(0,0), y=c(1,.8), label = c('A', 'B'))
fig_4


ggsave(fig_4,
       filename = './outputs/13150x6631.jpeg',
       width = 260,
       height = 290,
       device = 'jpeg',
       dpi = 300,
       units = 'mm')


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



# 
# 
# 
# gff_test <- 
#   read_tsv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gifrop_out/my_islands/SampleA_29_1.gff')
# gff_test %>% ggplot() + geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax=1), color='black') +
#   geom_text(aes(x=(start+end)/2, y=.5, label=product), angle=90)
# # 
# gff_test1 <- read_tsv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gifrop_out/my_islands/14229-5_17_1.gff')
# # gff_test2 <- read_tsv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gifrop_out/my_islands/6461_173_1.gff')
# gff_test3 <- read_tsv('/project/fsep_004/jtrachsel/Campy_HGT/outputs/pan_genomes_old/14229-5_6461_SampleA/pan/gifrop_out/my_islands/SampleA_29_1.gff') %>%  mutate(seqid=as.character(seqid))
# # 
# gff_test_all <- bind_rows(gff_test1,
#                           #gff_test2,
#                           gff_test3) %>%
#   mutate(genome=sub('(.*)_[0-9]+', '\\1', seqid),
#          gene_name=sub('.*','',attributes)) %>%
#   group_by(genome) %>%
#   mutate(yval=cur_group_id(),
#          is_start=start-min(start),
#          len=end-start,
#          is_end=is_start+len)
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

# tst <- RESULTS[[7]]
# pii <- read_csv('./outputs/pan_genomes/14229-5_6461_SampleG/gifrop_out/pan_with_island_info.csv')
# 
# 
# pii %>% filter(all)
# 
# 
# look <- gpa_long[2450:2470,]
# 
# which(gpa_long$locus_tags == 'BEGKIOMD_00456')
# 
# 
# 
# gpa_long_foc <- 
#   gpa_long %>% 
#   filter(`Genome Fragment` == 1) %>% 
#   filter(`Order within Fragment` > 203 & `Order within Fragment` < 245)
# 
# 
# 
# gpa_long_foc$RESISTANCE <- factor(gpa_long_foc$RESISTANCE)
# 
# gpa_long_foc$AMR <- sub("aminoglycoside O-phosphotransferase ", "", gpa_long_foc$PRODUCT)
# 
# unique(gpa_long_foc$PRODUCT)
# 
# p1 <- 
#   gpa_long_foc %>% 
#   filter(genome != '6461') %>% 
#   # filter(`Genome Fragment` == 1) %>% 
#   # filter(`Order within Fragment` > 1581 & `Order within Fragment` < 1675) %>% 
#   ggplot(aes(x=`Order within Fragment`, y=genome, group=genome)) + 
#   # geom_vline(xintercept = 1604, color='orange')+
#   # geom_vline(xintercept = 1619, color='orange')+
#   geom_point(aes(fill=classification), shape=22, size=3, na.rm = TRUE) +
#   geom_text(aes(label=`Accessory Order with Fragment`))+
#   geom_point(data=filter(gpa_long_foc, !is.na(RESISTANCE) &(`Order within Fragment` > 1581 & `Order within Fragment` < 1675 & genome != '6461')),aes(color=AMR), na.rm = TRUE,  show.legend = FALSE)+
#   scale_fill_viridis_d() + 
#   scale_color_brewer(palette = 'Set1') + 
#   theme_cowplot() + 
#   theme(legend.position = 'top')
# p1
# # need to abricate and merge into annot by locus tag?
# annot <- 
#   gpa_long_foc %>% 
#   # filter(`Genome Fragment` == 1) %>% 
#   filter(`Order within Fragment` > 1604 & `Order within Fragment` < 1619) %>% 
#   select(Gene, Annotation, `Order within Fragment`) %>% unique() %>% 
#   left_join(gpa_abric) %>% 
#   mutate(Annotation2 = ifelse(is.na(PRODUCT), Annotation, PRODUCT))
# 
# 
# p2 <- 
#   gpa_long_foc %>% 
#   filter(`Genome Fragment` == 1 & genome != '6461') %>% 
#   filter(`Order within Fragment` > 1604 & `Order within Fragment` < 1619) %>% 
#   ggplot(aes(x=`Order within Fragment`, y=genome, group=genome)) + 
#   geom_point(aes(fill=classification), shape=22, size=4.5, show.legend = FALSE) +
#   geom_point(data=
#                filter(gpa_long_foc, !is.na(AMR) &
#                         (`Order within Fragment` > 1604 &
#                            `Order within Fragment` < 1619) &
#                         genome != '6461'),
#              aes(color=AMR), size=2,na.rm = TRUE , shape=17, show.legend = FALSE)+
#   geom_text(data=annot, aes(x=`Order within Fragment`, y=-5.5, label=Annotation2), inherit.aes = FALSE, hjust='left', size=3.5) + 
#   scale_y_discrete(expand=expansion(mult=c(.75,.09))) + coord_flip()+
#   scale_fill_viridis_d() +
#   # scale_fill_brewer(palette = 'Set1')+
#   scale_color_brewer(palette = 'Set1') + 
#   theme_cowplot() + 
#   theme(axis.text.x = element_text(angle = -45, vjust = -.9)) + 
#   xlab('Order')
# 
# p2
# 
# 
# ########
# 
# gpa_long %>% filter(`Genome Fragment` == 1) %>% 
#   group_by(`Order within Fragment`, classification) %>% tally() %>% 
#   filter(classification == 'donor_result') %>% 
#   filter(n > 1)
