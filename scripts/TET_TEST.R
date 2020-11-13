library(Biostrings)

# 
# 6461 + 13150 = 6461x13150
# 6461 + 13150 = 6461x13150exp137

gpa <- read_csv('outputs/pan_genomes/6461_13150_6461x13150/pan/gene_presence_absence.csv')

look <- gpa[grep('tetracycline ', gpa$Annotation),]

# 13150: PGHKFEJK_01855
# 6461:  EKEIMHEI_00148

# 6461x13150:
# AGLJFKNC_00132
# AGLJFKNC_01891








gpa2 <- read_csv('outputs/pan_genomes/6461_13150_6461x13150exp137/pan/gene_presence_absence.csv')

look2 <- gpa2[grep('tetracycline ', gpa2$Annotation),]

# PGHKFEJK_01855
# EKEIMHEI_00148
# MIPFMNEC_01862

# everything near tet seems to be mostly the same, but then downstream quite a ways there is a gene with very differnt
# length
#PGHKFEJK_01847
#PGHKFEJK_01897
NA
MIPFMNEC_01871



# THESE ARE THE LOCUS TAGS THAT ARE VERY DIFFERENT LENGTHS
PGHKFEJK_01897
NA
MIPFMNEC_01871


# THESE ARE TETO
# PGHKFEJK_01855
# EKEIMHEI_00148
# MIPFMNEC_01862


gff13150 <- gff_parse3('./outputs/pan_genomes/6461_13150_6461x13150exp137/13150/13150.gff')
gff6461 <- gff_parse3('./outputs/pan_genomes/6461_13150_6461x13150exp137/6461/6461.gff')
gff6461x13150exp137 <- gff_parse3('./outputs/pan_genomes/6461_13150_6461x13150exp137/6461x13150exp137/6461x13150exp137.gff')




# 6461x13150_1	tet(O)	chromosome
# 6461x13150_1	blaOXA-605	chromosome
# 6461x13150_2	tet(O)	plasmid
# 6461x13150_2	aph(3')-IIIa	plasmid
# 6461x13150exp137_1	blaOXA-605	chromosome
# 6461x13150exp137_2	tet(O)	plasmid
# 6461x13150exp137_2	aph(3')-IIIa	plasmid
#







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


gff_parse3('./outputs/annotations//SampleA.gff')

read_tsv('./outputs/abricate/SampleA.fna.abricate')


read_tsv('./outputs/pan_genomes/all_6461s/')






6461 + 13150 = 6461x13150exp137