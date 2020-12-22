library(tidyverse)

gff_parse3 <- function(path){
  # only excludes 'gene' type annotations, might be better than to only allow CDS and trna...
  # added mutate statement to remove extra stuff that sometimes comes along with the locus tag
  gff <- read_delim(path,
                    delim = '\t',
                    col_names = c("seqid", "source", "type", "start", "end", "score", "strand","phase","attributes"),
                    comment = '#', progress = FALSE, col_types = c('cccddcccc')) %>%
    filter(type != 'gene') %>%
    tidyr::extract(attributes,
                   into = c('ID', 'locus_tag', 'product'),
                   regex ='ID=(.*);.*locus_tag=(.*_[0-9]+);.*product=(.*)',
                   remove = FALSE) %>%
    mutate(locus_tag = sub('([A-Za-z]_[0-9]+).*', '\\1', locus_tag))
  return(gff)
}



tst <- gff_parse3('./outputs/annotations/SampleA/SampleA.gff')

# islands should be split on the loc_tag_order column instead of num_loc_tag
tst1 <- tst %>%
  filter(seqid == 'SampleA_1') %>%
  filter(type == 'CDS') %>% 
  mutate(num_loc_tag=as.numeric(sub('(.*)_([0-9]+)','\\2',locus_tag))) %>% 
  arrange(num_loc_tag) %>% 
  mutate(loc_tag_order=seq_along(num_loc_tag))


seq_along(tst1$locus_tag)

tst1 %>% ggplot(aes(x=num_loc_tag, y=loc_tag_order)) + geom_point(size=.00001)
#





read_tsv('./outputs/abricate/SampleA.fna.abricate')


read_tsv('./outputs/pan_genomes/all_6461s/')



