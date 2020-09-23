library(ggrepel)
library(ggmsa)
library(tidyverse)
library(msa)
library(Biostrings)

ani_tab <- read_tsv('/home/Julian.Trachsel/Vanina/June29_2020/ALL/pyani_res/ANIm_percentage_identity.tab') %>% 
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


### END ANI MDS ###

#
# cii coltypes
#island_ID,acc_frag,primary_cluster,secondary_cluster,seqid,genome_name,start,end,island_length,num_genes,
# locus_tags,seqid_len,percent_island,only_island,island_type,RESISTANCE,res_type,vir_type,plasmid_type,viro_type
'ccccccddddcddlcccccc'

#gpai coltypes
# Gene,num_islands,num_Sclusts,all_islands,Pcluster,all_Sclusters,Non-unique Gene name,Annotation,
#No. isolates,No. sequences,Avg sequences per isolate,Genome Fragment,Order within Fragment,
#Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,
#Avg group size nuc,11601MD,6631,SampleP

'cddcccccdddcdcdcdddccc'

HGT_ID <- 
  function(base_dir, donor, recipient, result){
    # this selects gifrop identified genomic islands (by secondary cluster)
    # that are present in donor and the result.
    # (they may also be present in the recipient)
    
  cii_path <- paste(base_dir, 'gifrop_out/clustered_island_info.csv', sep = '')
  # gpai_path <- paste(base_dir, 'gifrop_out/pan_with_island_info.csv', sep = '')
  
  cii <- read_csv(cii_path, col_types = 'ccccccddddcddlcccccc')
  # gpai <- read_csv(gpai_path, col_types = 'cddcccccdddcdcdcdddccc')
  
  # potentials <- gpai %>% filter(`No. isolates` == 2)
  
  donor_islands <- cii %>% filter(genome_name == donor) %>%
    pull(secondary_cluster)
  
  recipient_islands <- cii %>% filter(genome_name == recipient) %>%
    pull(secondary_cluster)
  
  result_islands <- cii %>% filter(genome_name == result) %>%
    pull(secondary_cluster)
  
  # donor_islands %in% result_islands
  potentials <- donor_islands %in% result_islands
  potentials2  <- !(donor_islands %in% recipient_islands)
  # these <- donor_islands[potentials & potentials2]
  these <- donor_islands[potentials]
  return(cii %>% filter(secondary_cluster %in% these))
  
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

###

SampleA <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleA/pan/', 
       donor = '14229-5', 
       recipient = '6461', 
       result = 'SampleA') %>% 
  write_tsv('SampleA_HGT_results.tsv')


# CHANGED TO WITH FLANKING
SampleA_aln <- 
  HGT_aln(HGT_ID_tab = SampleA,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_f5k/gifrop_out/my_islands/with_flanking/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleA', 
          out_folder = 'ALIGNMENTS/')



## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleA,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleA', 
          out_folder = 'DONOR_RESULT_ALIGNMENTS/',
          FILT = c('14229-5', 'SampleA'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleA,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleA', 
          out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
          FILT = c('6461', 'SampleA'))
###

SampleB <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleB/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleB') %>% 
  write_tsv('SampleB_HGT_results.tsv')

SampleB_aln <- 
  HGT_aln(HGT_ID_tab = SampleB,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleB/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleB', 
          out_folder = 'ALIGNMENTS/')



## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleB,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleB/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleB', 
        out_folder = 'DONOR_RESULT_ALIGNMENTS/',
        FILT = c('14229-5', 'SampleB'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleB,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleB/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleB', 
        out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
        FILT = c('6461', 'SampleB'))


########



SampleC <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleC/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleC')%>% 
  write_tsv('SampleC_HGT_results.tsv')

SampleC_aln <- 
  HGT_aln(HGT_ID_tab = SampleC,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleC/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleC', 
          out_folder = 'ALIGNMENTS/')



## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleC,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleC/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleC', 
        out_folder = 'DONOR_RESULT_ALIGNMENTS/',
        FILT = c('14229-5', 'SampleC'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleC,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleC/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleC', 
        out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
        FILT = c('6461', 'SampleC'))



######
######
# SAMPLE E #
SampleE <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleE/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleE')%>% 
  write_tsv('SampleE_HGT_results.tsv')

SampleE_aln <- 
  HGT_aln(HGT_ID_tab = SampleE,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleE/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleE', 
          out_folder = 'ALIGNMENTS/')



## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleE,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleE/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleE', 
        out_folder = 'DONOR_RESULT_ALIGNMENTS/',
        FILT = c('14229-5', 'SampleE'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleE,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleE/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleE', 
        out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
        FILT = c('6461', 'SampleE'))





#






#

SampleF <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleF/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleF') %>%
  write_tsv('SampleF_HGT_results.tsv')

SampleF_aln <- 
  HGT_aln(HGT_ID_tab = SampleF,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleF/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleF', 
          out_folder = 'ALIGNMENTS/')



## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleF,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleF/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleF', 
        out_folder = 'DONOR_RESULT_ALIGNMENTS/',
        FILT = c('14229-5', 'SampleF'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleF,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleF/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleF', 
        out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
        FILT = c('6461', 'SampleF'))


#

SampleG <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleG/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleG')%>% 
  write_tsv('SampleG_HGT_results.tsv')

SampleG_aln <- 
  HGT_aln(HGT_ID_tab = SampleG,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleG/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleG', 
          out_folder = 'ALIGNMENTS/')

## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleG,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleG/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleG', 
        out_folder = 'DONOR_RESULT_ALIGNMENTS/',
        FILT = c('14229-5', 'SampleG'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleG,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleG/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleG', 
        out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
        FILT = c('6461', 'SampleG'))

#

SampleH <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleH/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleH')%>% 
  write_tsv('SampleH_HGT_results.tsv')

SampleH_aln <- 
  HGT_aln(HGT_ID_tab = SampleH,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleH/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleH', 
          out_folder = 'ALIGNMENTS/')



## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleH,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleH/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleH', 
        out_folder = 'DONOR_RESULT_ALIGNMENTS/',
        FILT = c('14229-5', 'SampleH'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleH,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleH/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleH', 
        out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
        FILT = c('6461', 'SampleH'))


SampleI <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleI/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleI')%>% 
  write_tsv('SampleI_HGT_results.tsv')

SampleI_aln <- 
  HGT_aln(HGT_ID_tab = SampleI,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleI/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleI', 
          out_folder = 'ALIGNMENTS/')



## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleI,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleI/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleI', 
        out_folder = 'DONOR_RESULT_ALIGNMENTS/',
        FILT = c('14229-5', 'SampleI'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleI,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleI/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleI', 
        out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
        FILT = c('6461', 'SampleI'))


#


SampleJ <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/hybrids/14229-5_6461_SampleJ/pan/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleJ') %>% 
  write_tsv('SampleJ_HGT_results.tsv')



SampleJ_aln <- 
  HGT_aln(HGT_ID_tab = SampleJ,
          islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleJ/pan/gifrop_flank/gifrop_out/my_islands/',
          donor = '14229-5', 
          recipient = '6461', 
          result = 'SampleJ', 
          out_folder = 'ALIGNMENTS/')



## this one only aligns donor and result ##

HGT_aln(HGT_ID_tab = SampleJ,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleJ/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleJ', 
        out_folder = 'DONOR_RESULT_ALIGNMENTS/',
        FILT = c('14229-5', 'SampleJ'))


## this one only aligns recipient and result ##

HGT_aln(HGT_ID_tab = SampleJ,
        islands_path ='~/Vanina/June29_2020/hybrids/14229-5_6461_SampleJ/pan/gifrop_flank/gifrop_out/my_islands/',
        donor = '14229-5', 
        recipient = '6461', 
        result = 'SampleJ', 
        out_folder = 'RECIPIENT_RESULT_ALIGNMENTS/',
        FILT = c('6461', 'SampleJ'))




#############


### TODO ###
# Align other non-AMR regions to check for potential recombination




### NEED TO DETERMINE DONOR RECIPIENT HERE ###
# THESE NEED ATTENTION
SampleM2 <- 
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/June29_2020/shorts/6461_6067_SampleM/pan/gifrop_flank', 
         donor = '6067', 
         recipient = '6461', 
         result = 'SampleM')%>% 
  write_tsv('SampleM_HGT_results.tsv')
# this didn't work...
SampleM_aln <- SampleM %>% 
  select(secondary_cluster, island_ID) %>% 
  group_by(secondary_cluster) %>% 
  summarise(IDs=list(island_ID)) %>% 
  mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_6067_SampleM/pan/gifrop_flank/gifrop_out/my_islands/'), 
         aln=map(.x=seqs, .f=msa, method='ClustalOmega'))

### SampleM redo ###



SampleM <-
  HGT_ID(base_dir = '/home/Julian.Trachsel/Vanina/July16/SampleM_redo/pan/',
         donor = '6067',
         recipient = '6461',
         result = 'SampleM')%>%
  write_tsv('SampleM_HGT_results_REDO.tsv')

###

gpa <- read_csv('/home/Julian.Trachsel/Vanina/July16/SampleM_redo/pan/gene_presence_absence.csv')


###

SampleN <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/shorts/6461_13150_SampleN/pan/gifrop_flank/', 
         donor = '6461', 
         recipient = '13150', 
         result = 'SampleN')%>% 
  write_tsv('SampleN_HGT_results.tsv')

SampleN_aln <- SampleN %>% 
  select(secondary_cluster, island_ID) %>% 
  group_by(secondary_cluster) %>% 
  summarise(IDs=list(island_ID)) %>% 
  mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_13150_SampleN/pan/gifrop_flank/gifrop_out/my_islands/'), 
         aln=map(.x=seqs, .f=msa, method='ClustalOmega'))

# this one is weird too...


SampleO <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/shorts/6461_14229-5_SampleO/pan/gifrop_flank/', 
         donor = '14229-5', 
         recipient = '6461', 
         result = 'SampleO')%>% 
  write_tsv('SampleO_HGT_results.tsv')
SampleO_aln <- SampleO %>% 
  select(secondary_cluster, island_ID) %>% 
  group_by(secondary_cluster) %>% 
  summarise(IDs=list(island_ID)) %>% 
  mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_14229-5_SampleO/pan/gifrop_flank/gifrop_out/my_islands/'), 
         aln=map(.x=seqs, .f=msa, method='ClustalOmega'))



SampleP <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/shorts/11601MD_6631_SampleP/pan/gifrop_flank/', 
         donor = '6631', 
         recipient = '11601MD', 
         result = 'SampleP')%>% 
  write_tsv('SampleP_HGT_results.tsv')
SampleP_aln <- SampleP %>% 
  select(secondary_cluster, island_ID) %>% 
  group_by(secondary_cluster) %>% 
  summarise(IDs=list(island_ID)) %>% 
  mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_14229-5_SampleP/pan/gifrop_flank/gifrop_out/my_islands/'), 
         aln=map(.x=seqs, .f=msa, method='ClustalOmega'))
# didnt work...

SampleQ <- 
  HGT_ID(base_dir = '~/Vanina/June29_2020/shorts/6461_13150_SampleQ/pan/gifrop_flank/', 
         donor = '6461', 
         recipient = '13150', 
         result = 'SampleQ')%>% 
  write_tsv('SampleQ_HGT_results.tsv')

SampleQ_aln <- SampleQ %>% 
  select(secondary_cluster, island_ID) %>% 
  group_by(secondary_cluster) %>% 
  summarise(IDs=list(island_ID)) %>% 
  mutate(seqs=map(.x=IDs,.f=read_islands,base_path='~/Vanina/June29_2020/shorts/6461_13150_SampleQ/pan/gifrop_flank/gifrop_out/my_islands/'), 
         aln=map(.x=seqs, .f=msa, method='ClustalOmega'))


SampleQ_aln %>% mutate(aln=map(.x=aln, .f=as, 'BStringSet'))
SampleQ_aln %>% mutate(aln=map(.x=aln, .f=as, 'DNAMultipleAlignment'))


ggmsa(as(test[[1]], "DNAMultipleAlignment"))

# this is not bad....
gff_test <- read_tsv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/14229-5_1_1.gff')
gff_test %>% ggplot() + geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax=1), color='black') + 
  geom_text(aes(x=(start+end)/2, y=.5, label=product), angle=90)

gff_test1 <- read_tsv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/14229-5_21_1.gff')
gff_test2 <- read_tsv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/6461_173_1.gff')
gff_test3 <- read_tsv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gifrop_out/my_islands/SampleA_72_1.gff')

gff_test_all <- bind_rows(gff_test1, gff_test2, gff_test3) %>%
  mutate(genome=sub('(.*)_[0-9]+', '\\1', seqid), 
         gene_name=sub('.*','',attributes)) %>% 
  group_by(genome) %>% 
  mutate(yval=cur_group_id(), 
         is_start=start-min(start), 
         len=end-start, 
         is_end=is_start+len)

gff_test_all$attributes

  ### getting there...
  
gpa_test <- read_csv('June29_2020/hybrids/14229-5_6461_SampleA/pan/gifrop_flank/gene_presence_absence.csv')%>%
  pivot_longer(cols = c(15:17), names_to='genome', values_to='locus_tag') %>% 
  select(Gene, locus_tag, Annotation)


gff_test_all <- gff_test_all %>% left_join(gpa_test, by = 'locus_tag')



gff_test_all %>% ggplot()+
  geom_rect(aes(xmin = is_start, xmax=is_end, ymin=yval, ymax=yval+1, fill=genome), color='black') + 
  geom_text(aes(x=(is_start+is_end)/2, y=yval+.5, label=Gene), angle=90)

colnames(gff_test_all)
#
# gpa_test %>% pivot_longer(cols = c(15:17), names_to='genome', values_to='locus_tag')



'14229-5_21_1'
'6461_173_1'
'SampleA_72_1'




test <- SampleQ_aln %>% pull(aln)
test2 <- unlist(test)
Biostrings::DNAMultipleAlignment()