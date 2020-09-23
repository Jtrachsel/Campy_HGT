library(ggrepel)
library(ggmsa)
library(tidyverse)
library(Biostrings)


##########

transferred_gene_table <- 
  function(base_dir, donor, recipient, result){
    
    pii_path <- paste0(base_dir, '/gifrop_out/pan_with_island_info.csv')
    Rtab_path <- paste0(base_dir, '/gene_presence_absence.Rtab')
    
    pii <- read_csv(pii_path)
    
    Rtab <- read_tsv(Rtab_path)%>% 
      select(Gene, recipient, donor, result) %>% 
      column_to_rownames(var = 'Gene')
    
    THESE <- rownames(Rtab[(Rtab[,1] == 0) & (Rtab[,3] > 0) & (Rtab[,2] > 0),])
    
    res <- pii %>% filter(Gene %in% THESE)
    
    return(res)
    
  }


HGT_ID <- 
  function(base_dir, donor, recipient, result, clust_level){
    # this selects gifrop identified genomic islands (by secondary cluster)
    # that are present in donor and the result.
    # (they may also be present in the recipient)
    clust_level <- enquo(clust_level)
    
    cii_path <- paste(base_dir, 'gifrop_out/clustered_island_info.csv', sep = '')
    # gpai_path <- paste(base_dir, 'gifrop_out/pan_with_island_info.csv', sep = '')
    
    cii <- read_csv(cii_path, col_types = 'ccccccccddddcddlccccccc')
    # gpai <- read_csv(gpai_path, col_types = 'cddcccccdddcdcdcdddccc')
    
    # potentials <- gpai %>% filter(`No. isolates` == 2)
    
    donor_islands <- cii %>% filter(genome_name == donor) %>%
      pull(!!clust_level)
    
    recipient_islands <- cii %>% filter(genome_name == recipient) %>%
      pull(!!clust_level)
    
    result_islands <- cii %>% filter(genome_name == result) %>%
      pull(!!clust_level)
    
    # donor_islands %in% result_islands
    potentials <- donor_islands %in% result_islands
    potentials2  <- !(donor_islands %in% recipient_islands)
    # these <- donor_islands[potentials & potentials2]
    these <- donor_islands[potentials]
    # browser()
    cii <- cii %>%
      mutate(HGT_role=
               case_when(
                 genome_name == donor     ~ 'donor', 
                 genome_name == recipient ~ 'recipient', 
                 genome_name == result    ~ 'result'
               ))
    # browser()
    return(cii %>% filter(!!clust_level %in% these))
    
}


plot_islands <- function(HGT_ID_result){
  HGT_ID_result %>% 
    select(genome_name, start, end, quat_cluster, HGT_role, seqid_len) %>% 
    arrange(quat_cluster) %>% 
    pivot_wider(names_from = HGT_role, values_from=c(start, end, genome_name, seqid_len)) %>% 
    ggplot() + 
    geom_segment(aes(x=1,
                     xend=seqid_len_donor,
                     y=genome_name_donor, 
                     yend=genome_name_donor), color='grey') + 
    geom_segment(aes(x=1,
                     xend=seqid_len_result,
                     y=genome_name_result, 
                     yend=genome_name_result), color='grey') + 
    geom_segment(aes(x=start_donor, xend=end_donor,
                     y=genome_name_donor, yend=genome_name_donor, 
                     color=quat_cluster), size=5) + 
    geom_segment(aes(x=start_result, xend=end_result,
                     y=genome_name_result, yend=genome_name_result, 
                     color=quat_cluster), size=5) + 
    geom_segment(aes(x=start_result, xend=end_donor,
                     y=genome_name_result, yend=genome_name_donor, 
                     color=quat_cluster)) +
    ylab('Genome Name') + xlab('Genomic coordinates')
}

read_islands <-
  function(vector_island_ids, base_path){
    # base path is a path to a 'my_islands' directory output by gifrop
  vector_of_fasta_paths <- paste(base_path, vector_island_ids, '.fasta', sep = '')
  island_seqs_list <- sapply(vector_of_fasta_paths, readDNAStringSet)
  names(island_seqs_list) <- NULL
  island_seqs_set <- do.call(c, island_seqs_list)
  return(island_seqs_set)
  }


HGT_aln <-
  function(HGT_ID_tab,       # a tibble output by HGT_ID function
           islands_path,     # a path to a 'my_islands' directory output by gifrop
           aln_method='ClustalOmega',
           out_folder=NULL,  # if specified will write alignments to this folder
           donor=NULL,       # character, name of donor genome
           result,           # character, name of result genome
           recipient=NULL,   # character, name of recipient genome
           FILT = NULL){     # character vector of genomes to limit alignments to

    if (!is.null(FILT)){
      HGT_ID_tab <- HGT_ID_tab %>% filter(genome_name %in% FILT)
    }


    HGT_ID_tab_aln <- HGT_ID_tab %>%
      select(secondary_cluster, island_ID) %>%
      group_by(secondary_cluster) %>%
      summarise(IDs=list(island_ID),
                num_to_aln=length(island_ID)) %>%
      filter(num_to_aln > 1) %>%
      mutate(seqs=map(.x=IDs,.f=read_islands,base_path=islands_path),
             aln=map(.x=seqs, .f=msa, method='ClustalOmega'),
             aln_unmask=map(.x=aln, .f=unmasked))

    names(HGT_ID_tab_aln$aln_unmask) <- paste(result, '_', HGT_ID_tab_aln$secondary_cluster,'.aln', sep = '')


    if (is.null(out_folder)){
      sapply(names(HGT_ID_tab_aln$aln_unmask),
                       function (x) writeXStringSet(HGT_ID_tab_aln$aln_unmask[[x]], filepath =  x))
    } else {
      sapply(names(HGT_ID_tab_aln$aln_unmask),
                       function (x) writeXStringSet(HGT_ID_tab_aln$aln_unmask[[x]], filepath =  paste(out_folder, x, sep = '')))

    }

    return(HGT_ID_tab_aln)

  }



#########


ani_tab <- read_tsv('/home/Julian.Trachsel/Vanina/outputs/pyani_res/ANIm_percentage_identity.tab') %>% 
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


parental_genomes <- c('11601MD', '6631', '14229-5', '6461', '6067', '13150')



  ani_tab %>% 
  transmute(parental_genome = X1,
            result_genome = genome,
            ani_dist=ani_dist) %>% 
  filter(parental_genome %in% parental_genomes) %>% 
  filter(!(result_genome %in% parental_genomes)) %>% 
  group_by(result_genome) %>% 
  summarise(recipient_genome=parental_genome[which.min(ani_dist)]) %>% 
  mutate(donor_genome=case_when(
    result_genome == '11601MDx6631'        ~ '6631', 
    result_genome == '6461x13150'          ~ '6461', 
    result_genome == '6461x13150_exp137'   ~ '6461', 
    result_genome == '6461x14229-5'        ~ '14229-5', 
    result_genome == '6461x6067'           ~ '6067', 
    grepl('Sample', result_genome)         ~ '14229-5')) %>%
  mutate(pan_dir=paste(donor_genome,
                       recipient_genome,
                       result_genome,
                       sep = '_',
                       collapse = '_'))


metadata
# now I need to make a directory for each result genome
###



### END ANI MDS ###





###

#this one is weird....
# no AMR or anything...

# 11601MDx6631          11601MD = recipient
base_dir <- '/home/Julian.Trachsel/Vanina/Sept2020/11601MD_6631_11601MDx6631/pan99/'

split_11601MDx6631 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/11601MD_6631_11601MDx6631/pan/', 
         donor = '6631', 
         recipient = '11601MD', 
         result = '11601MDx6631', 
         clust_level = secondary_cluster)

plot_islands(split_11601MDx6631)

look <- transferred_gene_table(base_dir = base_dir, 
                       donor = '6631', 
                       recipient = '11601MD', 
                       result = '11601MDx6631')
######################################################################
######################################################################
# THis one seems like there were two recombination events. SEE FIGURE
# 6461x13150            13150 = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan/'

split_6461x13150 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan99/', 
         donor = '6461', 
         recipient = '13150', 
         result = '6461x13150', 
         clust_level = tertiary_cluster)


plot_islands(split_6461x13150)

# 6461x13150_exp137     13150 = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150_exp137/pan/'

split_6461x13150_exp137 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150_exp137/pan/', 
         donor = '6461', 
         recipient = '13150', 
         result = '6461x13150_exp137')


plot_islands(split_6461x13150_exp137)


#########

# 6461x14229-5          6461  = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_14229-5_6461x14229-5/pan/'

split_6461x14229_5 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_14229-5_6461x14229-5/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = '6461x14229-5')



plot_islands(split_6461x14229_5)

#############

# 
# 
# split_6461x13150 <- 
#   HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan/', 
#        donor = '6461', 
#        recipient = '13150', 
#        result = '6461x13150')
# 
# plot_islands(split_6461x13150)
# 

########


# 6461x6067             6461  = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_6067_6461x6067/pan/'

split_6461x6067 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_6067_6461x6067/pan/', 
         donor = '6067', 
         recipient = '6461', 
         result = '6461x6067')


plot_islands(split_6461x6067)






########### pan 99  ##########

split_11601MDx6631 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/11601MD_6631_11601MDx6631/pan99/', 
         donor = '6631', 
         recipient = '11601MD', 
         result = '11601MDx6631')

plot_islands(split_11601MDx6631)



######################################################################
######################################################################
# THis one seems like there were two recombination events. SEE FIGURE
# 6461x13150            13150 = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan/'

split_6461x13150 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan99/', 
         donor = '6461', 
         recipient = '13150', 
         result = '6461x13150')


plot_islands(split_6461x13150)

# 6461x13150_exp137     13150 = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150_exp137/pan/'

split_6461x13150_exp137 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150_exp137/pan99/', 
         donor = '6461', 
         recipient = '13150', 
         result = '6461x13150_exp137')


plot_islands(split_6461x13150_exp137)


#########

# 6461x14229-5          6461  = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_14229-5_6461x14229-5/pan/'

split_6461x14229_5 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_14229-5_6461x14229-5/pan99/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = '6461x14229-5')



plot_islands(split_6461x14229_5)

#############

# 
# 
# split_6461x13150 <- 
#   HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan/', 
#        donor = '6461', 
#        recipient = '13150', 
#        result = '6461x13150')
# 
# plot_islands(split_6461x13150)
# 

########


# 6461x6067             6461  = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_6067_6461x6067/pan/'

split_6461x6067 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_6067_6461x6067/pan99/', 
         donor = '6067', 
         recipient = '6461', 
         result = '6461x6067')


plot_islands(split_6461x6067)



#
########### pan 100  ##########

split_11601MDx6631 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/11601MD_6631_11601MDx6631/pan100/', 
         donor = '6631', 
         recipient = '11601MD', 
         result = '11601MDx6631')

plot_islands(split_11601MDx6631)


######################################################################
######################################################################
# THis one seems like there were two recombination events. SEE FIGURE
# 6461x13150            13150 = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan/'

split_6461x13150 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan100/', 
         donor = '6461', 
         recipient = '13150', 
         result = '6461x13150')


plot_islands(split_6461x13150)

# 6461x13150_exp137     13150 = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150_exp137/pan/'

split_6461x13150_exp137 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150_exp137/pan100/', 
         donor = '6461', 
         recipient = '13150', 
         result = '6461x13150_exp137')


plot_islands(split_6461x13150_exp137)


#########

# 6461x14229-5          6461  = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_14229-5_6461x14229-5/pan/'

split_6461x14229_5 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_14229-5_6461x14229-5/pan100/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = '6461x14229-5')



plot_islands(split_6461x14229_5)

#############

# 
# 
# split_6461x13150 <- 
#   HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan/', 
#        donor = '6461', 
#        recipient = '13150', 
#        result = '6461x13150')
# 
# plot_islands(split_6461x13150)
# 

########


# 6461x6067             6461  = recipient
'/home/Julian.Trachsel/Vanina/Sept2020/6461_6067_6461x6067/pan/'

split_6461x6067 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/Sept2020/6461_6067_6461x6067/pan100/', 
         donor = '6067', 
         recipient = '6461', 
         result = '6461x6067')


plot_islands(split_6461x6067)

#

### Summary ###
# I used a pan-genome framework to identify genes that were transferred from the donor
# strain to the host strain.  In this framework I calculate a pan-genome from 
# 3 genomes, the donor, the recipient, and the result. Any gene that occurs in all three
# is a 'core' gene, any gene that occurs in only the donor and the result strain is considered
# to be a transferred gene.
# This analysis identified several instances of multiple genes being transferred in close proximity.
# but these non-core genes are interrupted by core genes.
# I believe this is evidence of a homologous recombination event.
# I think these events occur where the Donor genome has close homology 
# with a segment of the recipient genome, but has some novel genes inserted into this 
# region of homology.  
# The result strain gains these genes through homologous recombination


# roary Pangenomes determine what genes are the 'same gene' based on protein sequence identity
# the default is 95% protein ID
# I ran this same analysis at 95, 99, and 100% id.
# 100% id may be most sensitive but may suffer from sequencing errors introducing false pos
# 95% id will miss the transfer of very similar, but different, versions of the same gene. (allele transfer)


# # CHANGED TO WITH FLANKING
# SampleA_aln <- 
#   HGT_aln(HGT_ID_tab = SampleA,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_f5k/gifrop_out/my_islands/with_flanking/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleA', 
#           out_folder = 'ALIGNMENTS/')
# 
# 
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleA,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleA', 
#           out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#           FILT = c('14229-5', 'SampleA'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleA,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleA', 
#           out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#           FILT = c('6461', 'SampleA'))
# ###
# 
# SampleB <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleB/pan/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleB') %>% 
#   write_tsv('SampleB_HGT_results.tsv')
# 
# SampleB_aln <- 
#   HGT_aln(HGT_ID_tab = SampleB,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleB/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleB', 
#           out_folder = 'ALIGNMENTS/')
# 
# 
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleB,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleB/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleB', 
#         out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#         FILT = c('14229-5', 'SampleB'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleB,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleB/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleB', 
#         out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#         FILT = c('6461', 'SampleB'))
# 
# 
# ########
# 
# 
# 
# SampleC <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleC/pan/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleC')%>% 
#   write_tsv('SampleC_HGT_results.tsv')
# 
# SampleC_aln <- 
#   HGT_aln(HGT_ID_tab = SampleC,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleC/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleC', 
#           out_folder = 'ALIGNMENTS/')
# 
# 
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleC,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleC/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleC', 
#         out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#         FILT = c('14229-5', 'SampleC'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleC,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleC/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleC', 
#         out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#         FILT = c('6461', 'SampleC'))
# 
# 
# 
# ######
# ######
# # SAMPLE E #
# SampleE <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleE/pan/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleE')%>% 
#   write_tsv('SampleE_HGT_results.tsv')
# 
# SampleE_aln <- 
#   HGT_aln(HGT_ID_tab = SampleE,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleE/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleE', 
#           out_folder = 'ALIGNMENTS/')
# 
# 
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleE,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleE/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleE', 
#         out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#         FILT = c('14229-5', 'SampleE'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleE,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleE/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleE', 
#         out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#         FILT = c('6461', 'SampleE'))
# 
# 
# 
# 
# 
# #
# 
# 
# 
# 
# 
# 
# #
# 
# SampleF <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleF/pan/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleF') %>%
#   write_tsv('SampleF_HGT_results.tsv')
# 
# SampleF_aln <- 
#   HGT_aln(HGT_ID_tab = SampleF,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleF/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleF', 
#           out_folder = 'ALIGNMENTS/')
# 
# 
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleF,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleF/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleF', 
#         out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#         FILT = c('14229-5', 'SampleF'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleF,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleF/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleF', 
#         out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#         FILT = c('6461', 'SampleF'))
# 
# 
# #
# 
# SampleG <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleG/pan/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleG')%>% 
#   write_tsv('SampleG_HGT_results.tsv')
# 
# SampleG_aln <- 
#   HGT_aln(HGT_ID_tab = SampleG,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleG/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleG', 
#           out_folder = 'ALIGNMENTS/')
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleG,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleG/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleG', 
#         out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#         FILT = c('14229-5', 'SampleG'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleG,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleG/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleG', 
#         out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#         FILT = c('6461', 'SampleG'))
# 
# #
# 
# SampleH <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleH/pan/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleH')%>% 
#   write_tsv('SampleH_HGT_results.tsv')
# 
# SampleH_aln <- 
#   HGT_aln(HGT_ID_tab = SampleH,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleH/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleH', 
#           out_folder = 'ALIGNMENTS/')
# 
# 
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleH,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleH/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleH', 
#         out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#         FILT = c('14229-5', 'SampleH'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleH,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleH/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleH', 
#         out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#         FILT = c('6461', 'SampleH'))
# 
# 
# SampleI <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleI/pan/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleI')%>% 
#   write_tsv('SampleI_HGT_results.tsv')
# 
# SampleI_aln <- 
#   HGT_aln(HGT_ID_tab = SampleI,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleI/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleI', 
#           out_folder = 'ALIGNMENTS/')
# 
# 
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleI,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleI/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleI', 
#         out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#         FILT = c('14229-5', 'SampleI'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleI,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleI/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleI', 
#         out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#         FILT = c('6461', 'SampleI'))
# 
# 
# #
# 
# 
# SampleJ <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleJ/pan/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleJ') %>% 
#   write_tsv('SampleJ_HGT_results.tsv')
# 
# 
# 
# SampleJ_aln <- 
#   HGT_aln(HGT_ID_tab = SampleJ,
#           islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleJ/pan/gifrop_flank/gifrop_out/my_islands/',
#           donor = '14229-5', 
#           recipient = '6461', 
#           result = 'SampleJ', 
#           out_folder = 'ALIGNMENTS/')
# 
# 
# 
# ## this one only aligns donor and result ##
# 
# HGT_aln(HGT_ID_tab = SampleJ,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleJ/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleJ', 
#         out_folder = 'DONOR_RESULT_ALIGNMENTS/',
#         FILT = c('14229-5', 'SampleJ'))
# 
# 
# ## this one only aligns recipient and result ##
# 
# HGT_aln(HGT_ID_tab = SampleJ,
#         islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleJ/pan/gifrop_flank/gifrop_out/my_islands/',
#         donor = '14229-5', 
#         recipient = '6461', 
#         result = 'SampleJ', 
#         out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
#         FILT = c('6461', 'SampleJ'))
# 
# 
# 
# 
# #############
# 
# 
# ### TODO ###
# # Align other non-AMR regions to check for potential recombination
# 
# 
# 
# 
# ### NEED TO DETERMINE DONOR RECIPIENT HERE ###
# # THESE NEED ATTENTION
# SampleM2 <- 
#   HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/June29_2020/shorts/6461_6067_SampleM/pan/gifrop_flank', 
#          donor = '6067', 
#          recipient = '6461', 
#          result = 'SampleM')%>% 
#   write_tsv('SampleM_HGT_results.tsv')
# # this didn't work...
# SampleM_aln <- SampleM %>% 
#   select(secondary_cluster, island_ID) %>% 
#   group_by(secondary_cluster) %>% 
#   summarise(IDs=list(island_ID)) %>% 
#   mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_6067_SampleM/pan/gifrop_flank/gifrop_out/my_islands/'), 
#          aln=map(.x=seqs, .f=msa, method='ClustalOmega'))
# 
# ### SampleM redo ###
# 
# 
# 
# SampleM <-
#   HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/July16/SampleM_redo/pan/',
#          donor = '6067',
#          recipient = '6461',
#          result = 'SampleM')%>%
#   write_tsv('SampleM_HGT_results_REDO.tsv')
# 
# ###
# 
# gpa <- read_csv('/home/Julian.Trachsel/Vanina/July16/SampleM_redo/pan/gene_presence_absence.csv')
# 
# 
# ###
# 
# SampleN <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/shorts/6461_13150_SampleN/pan/gifrop_flank/', 
#          donor = '6461', 
#          recipient = '13150', 
#          result = 'SampleN')%>% 
#   write_tsv('SampleN_HGT_results.tsv')
# 
# SampleN_aln <- SampleN %>% 
#   select(secondary_cluster, island_ID) %>% 
#   group_by(secondary_cluster) %>% 
#   summarise(IDs=list(island_ID)) %>% 
#   mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_13150_SampleN/pan/gifrop_flank/gifrop_out/my_islands/'), 
#          aln=map(.x=seqs, .f=msa, method='ClustalOmega'))
# 
# # this one is weird too...
# 
# 
# SampleO <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/shorts/6461_14229-5_SampleO/pan/gifrop_flank/', 
#          donor = '14229-5', 
#          recipient = '6461', 
#          result = 'SampleO')%>% 
#   write_tsv('SampleO_HGT_results.tsv')
# SampleO_aln <- SampleO %>% 
#   select(secondary_cluster, island_ID) %>% 
#   group_by(secondary_cluster) %>% 
#   summarise(IDs=list(island_ID)) %>% 
#   mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_14229-5_SampleO/pan/gifrop_flank/gifrop_out/my_islands/'), 
#          aln=map(.x=seqs, .f=msa, method='ClustalOmega'))
# 
# 
# 
# SampleP <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/shorts/11601MD_6631_SampleP/pan/gifrop_flank/', 
#          donor = '6631', 
#          recipient = '11601MD', 
#          result = 'SampleP')%>% 
#   write_tsv('SampleP_HGT_results.tsv')
# SampleP_aln <- SampleP %>% 
#   select(secondary_cluster, island_ID) %>% 
#   group_by(secondary_cluster) %>% 
#   summarise(IDs=list(island_ID)) %>% 
#   mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_14229-5_SampleP/pan/gifrop_flank/gifrop_out/my_islands/'), 
#          aln=map(.x=seqs, .f=msa, method='ClustalOmega'))
# # didnt work...
# 
# SampleQ <- 
#   HGT_ID(base_dir = '~/Vanina/June29_2020/shorts/6461_13150_SampleQ/pan/gifrop_flank/', 
#          donor = '6461', 
#          recipient = '13150', 
#          result = 'SampleQ')%>% 
#   write_tsv('SampleQ_HGT_results.tsv')
# 
# SampleQ_aln <- SampleQ %>% 
#   select(secondary_cluster, island_ID) %>% 
#   group_by(secondary_cluster) %>% 
#   summarise(IDs=list(island_ID)) %>% 
#   mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_13150_SampleQ/pan/gifrop_flank/gifrop_out/my_islands/'), 
#          aln=map(.x=seqs, .f=msa, method='ClustalOmega'))
# 
# 
# SampleQ_aln %>% mutate(aln=map(.x=aln, .f=as, 'BStringSet'))
# SampleQ_aln %>% mutate(aln=map(.x=aln, .f=as, 'DNAMultipleAlignment'))
# 
# 
# ggmsa(as(test[[1]], "DNAMultipleAlignment"))
# 
# # this is not bad....
# gff_test <- read_tsv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/14229-5_1_1.gff')
# gff_test %>% ggplot() + geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax=1), color='black') + 
#   geom_text(aes(x=(start+end)/2, y=.5, label=product), angle=90)
# 
# gff_test1 <- read_tsv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/14229-5_21_1.gff')
# gff_test2 <- read_tsv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/6461_173_1.gff')
# gff_test3 <- read_tsv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/SampleA_72_1.gff')
# 
# gff_test_all <- bind_rows(gff_test1, gff_test2, gff_test3) %>%
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
#   
# gpa_test <- read_csv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gene_presence_absence.csv')%>%
#   pivot_longer(cols = c(15:17), names_to='genome', values_to='locus_tag') %>% 
#   select(Gene, locus_tag, Annotation)
# 
# 
# gff_test_all <- gff_test_all %>% left_join(gpa_test, by = 'locus_tag')
# 
# 
# 
# gff_test_all %>% ggplot()+
#   geom_rect(aes(xmin = is_start, xmax=is_end, ymin=yval, ymax=yval+1, fill=genome), color='black') + 
#   geom_text(aes(x=(is_start+is_end)/2, y=yval+.5, label=Gene), angle=90)
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
# 
# 
# test <- SampleQ_aln %>% pull(aln)
# test2 <- unlist(test)
# Biostrings::DNAMultipleAlignment()




###### FUNCTIONIZE THIS #######


# 
# pii <- read_csv('/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan/gifrop_out/pan_with_island_info.csv')
# 
# LOOK <- read_tsv('/home/Julian.Trachsel/Vanina/Sept2020/6461_13150_6461x13150/pan/gene_presence_absence.Rtab') %>% 
#   column_to_rownames(var = 'Gene')
# 
# THESE <- LOOK[(LOOK[,1] == 0) & (LOOK[,3] > 0) & (LOOK[,2] > 0),]
# 
# THESE <- rownames(THESE)
# res <- pii %>% filter(Gene %in% THESE)
###################################

# 
# split_6461x13150 %>% 
#   select(genome_name, start, end, quat_cluster, HGT_role, seqid_len) %>% 
#   arrange(quat_cluster) %>% 
#   pivot_wider(names_from = HGT_role, values_from=c(start, end, genome_name, seqid_len)) %>% 
#   ggplot() + 
#     geom_segment(aes(x=1,
#                      xend=seqid_len_donor,
#                      y=genome_name_donor, 
#                      yend=genome_name_donor), color='grey') + 
#   geom_segment(aes(x=1,
#                    xend=seqid_len_result,
#                    y=genome_name_result, 
#                    yend=genome_name_result), color='grey') + 
#   geom_segment(aes(x=start_donor, xend=end_donor,
#                    y=genome_name_donor, yend=genome_name_donor, 
#                    color=quat_cluster), size=5) + 
#   geom_segment(aes(x=start_result, xend=end_result,
#                    y=genome_name_result, yend=genome_name_result, 
#                    color=quat_cluster), size=5) + 
#   geom_segment(aes(x=start_result, xend=end_donor,
#                    y=genome_name_result, yend=genome_name_donor, 
#                    color=quat_cluster)) 


###########