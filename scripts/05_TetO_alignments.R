library(tidyverse)
library(Biostrings)

abric <- 
  read_tsv('./outputs/pan_genomes/tetO/pan_genome_abricate.ncbi') %>%
  mutate(Gene=SEQUENCE)

gpa <- read_csv('./outputs/pan_genomes/tetO/gene_presence_absence.csv') %>% 
  left_join(abric)


unique(gpa$RESISTANCE)
# look <- gpa %>% filter(`Genome Fragment` ==2)
look <- gpa[grep('TETRACYCLINE', gpa$RESISTANCE),]

THESE <- 
  look %>%
  gather(c(15:23), key = 'genome', value='locus_tag') %>% 
  select(Gene, Annotation, `Genome Fragment`, `Order within Fragment`, genome, locus_tag) %>% 
  na.omit() %>% 
  # mutate(genome_frag=ifelse(`Genome Fragment` == 77, 'chrom', 'plasmid')) %>% 
  mutate(ID=paste(genome, locus_tag, sep = '_'))

names(THESE$locus_tag) <- THESE$ID

name_swap <- THESE$ID
names(name_swap) <- THESE$locus_tag

gens <- THESE$genome %>% unique() %>% paste(collapse = '|')


faas <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*faa')
faas <- faas[!grepl('proteins', faas)]
faas <- faas[grepl(gens, faas)]


ffns <- list.files(path = './outputs/annotations/', recursive = T, full.names = T, pattern = '*ffn')
ffns <- ffns[!grepl('proteins', ffns)]
ffns <- ffns[grepl(gens, ffns)]


# faas <- list.files(path = './outputs/annotations', pattern='.*13150.*faa', recursive = TRUE, full.names = T)


# faas <- faas[-c(2:3)]
# faas <- c(faas, './outputs/annotations/6461/6461.faa')

nucseqs <- lapply(ffns, readDNAStringSet)
nucseqs <- do.call(c, nucseqs)

aaseqs <- lapply(faas, readAAStringSet)
aaseqs <- do.call(c, aaseqs)


names(aaseqs) <- sub('(.*_[0-9]+) .*','\\1',names(aaseqs))
names(nucseqs) <- sub('(.*_[0-9]+) .*','\\1',names(nucseqs))




TetO <- aaseqs[THESE$locus_tag]
TetOnuc <- nucseqs[THESE$locus_tag]




names(TetO) <- name_swap[names(TetO)]
names(TetOnuc) <- name_swap[names(TetOnuc)]


width(TetO)
width(TetOnuc)

691*3+3

# 2076 is TetO
# 1920 is TetO


# TetO[width(TetO) == 691]


# TetO[width(TetO) != 691]

library(msa)


# msa::msaClustalOmega(TetO[width(TetO) != 691])

# 
# TetO[width(TetO) == 691] %>% 
#   msa::msaClustalOmega() %>% 
#   msa::msaPrettyPrint(file = 'Elongation_factor_G.pdf')



TetOs <- TetO[width(TetO) %in% c(639, 561)] 



TetOnucs <- TetOnuc[width(TetOnuc) == 1920] 



names(TetOs)
tet_order <- c('13150_FAOCEDLB_01881', '6461x13150_DMIKNJMH_01884',
               'GCA_004947465_ELGLOCID_01772','GCA_004947585_PJNBKKGD_01276','GCA_005266575_DBOCFHAN_01484',
               '6461x13150exp137_DOLACMNP_01878','6461x13150_DMIKNJMH_00132','6461_EKEIMHEI_00148', 'JCC_BLFAMBHF_01924')


actives <- c('6461x13150exp137_DOLACMNP_01878','6461x13150_DMIKNJMH_00132','6461_EKEIMHEI_00148')

# GCA 004945185 IKGHOOBL 01821

TetOs <- TetOs[match(tet_order, names(TetOs))]
TetOnucs <- TetOnucs[match(tet_order, names(TetOnucs))]



# 13150_FAOCEDLB_01881 inactive
# 6461x13150_DMIKNJMH_01884 inactive
# GCA_004945185_IKGHOOBL_01821 inactive  # INCOMPLETE, REMOVED

# GCA_004947465_ELGLOCID_01772 inactive
# GCA_004947585_PJNBKKGD_01276 inactive
# GCA_005266575_DBOCFHAN_01484 inactive

# 6461x13150exp137_DOLACMNP_01878 active
# 6461x13150_DMIKNJMH_00132 active chrome
# 6461_EKEIMHEI_00148 active

TetOnucs %>% 
  msa::msaClustalOmega() %>% 
  msa::msaPrettyPrint(file = 'TetO_nuc.pdf')

TetOs %>% 
  msa::msaClustalOmega() %>% 
  msa::msaPrettyPrint(file = 'TetO_aa.pdf')

tst <- 
  stringDist(TetOs) %>%
  as.matrix() %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'from') %>% 
  pivot_longer(cols = -from, names_to='to') %>% 
  filter(from != to) %>% 
  filter(to %in% actives) %>% 
  filter(!(from %in% actives)) %>% 
  group_by(from) %>% 
  summarise(best_active=to)

rownames(tst)
colnames(tst)
# 
# 13150	 PGHKFEJK_00450
# 6461x13150	 ICGOJEKE_00453
# Exp137 	OHCLGAKK_00448


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


readAAStringSet('./outputs/annotations/6461x13150exp137/6461x13150exp137.faa')

BC07_gff <- gff_parse3('./outputs/annotations/6461x13150/6461x13150.gff') %>%
  mutate(genome='BC07')

exp137_gff <- gff_parse3('./outputs/annotations/6461x13150exp137/6461x13150exp137.gff') %>% 
  mutate(genome='exp137')

o13150_gff <- gff_parse3('./outputs/annotations/13150/13150.gff') %>% 
  mutate(genome='13150')



camp_plasmid <- rbind(BC07_gff, exp137_gff ,o13150_gff) %>% 
  filter(grepl('_2',seqid)) %>% 
  mutate(ymin=ifelse(genome =='13150', 0, 
                     ifelse(genome == 'BC07', 2, 4)), 
         ymax=ymin+1)

camp_plasmid %>%
  group_by(product,start, end) %>%
  tally() %>% ungroup() %>% 
  arrange((n))

camp_plasmid %>% group_by(genome, product) %>% 
  summarise(len=end-start)

camp_plasmid %>% ggplot() + 
  geom_rect(aes(ymin=start, ymax=end, xmin=ymin, xmax=ymax), color='black')


# genes are all identical lengths...
library(Biostrings)


exp137 <- readDNAStringSet('./outputs/annotations/6461x13150exp137/6461x13150exp137.fna')

exp137[2] %>% writeXStringSet(filepath = './outputs/exp137_50kb_plas.fna')



o13150 <- readDNAStringSet('./outputs/annotations/13150/13150.fna')

o13150[2] %>% writeXStringSet(filepath = './outputs/13150_50kb_plas.fna')


TETRES13150 <- readDNAStringSet('./outputs/annotations/6461x13150/6461x13150.fna')

TETRES13150[2] %>% writeXStringSet(filepath = './outputs/tetres_13150_50kb_plas.fna')

### where is tetO in the GCA genomes

teto_cii <- read_csv('./outputs/pan_genomes/tetO/gifrop_out/clustered_island_info.csv')


LOOK <- 
  teto_cii %>% filter(
  grepl('TETRACYCLINE', RESISTANCE) &
    grepl('GCA', genome_name)
)
LOOK %>% select(genome_name, RESISTANCE, only_island) %>% 
  transmute(genome=genome_name, 
            RESISTANCE=RESISTANCE, 
            location=ifelse(only_island, 'probably plasmid', 'chromosome'))

# 
# read_tsv('outputs/SNPs_13150/exp137/exp137.tab')
# 

### old below ####
# 
# 6461 + 13150 = 6461x13150
# 6461 + 13150 = 6461x13150exp137
# 
# gpa <- read_csv('outputs/pan_genomes/6461_13150_6461x13150/pan/gene_presence_absence.csv')
# 
# look <- gpa[grep('tetracycline ', gpa$Annotation),]

# 13150: PGHKFEJK_01855
# 6461:  EKEIMHEI_00148

# 6461x13150:
# AGLJFKNC_00132
# AGLJFKNC_01891




# 
# 
# 
# 
# gpa2 <- read_csv('./outputs/pan_genomes/6461_13150_6461x13150exp137/gene_presence_absence.csv')
# 
# look2 <- gpa2[grep('tetracycline ', gpa2$Annotation),]
# 
# 
# cii <- read_csv('./outputs/pan_genomes/6461_13150_6461x13150exp137/gifrop_out/clustered_island_info.csv')
# 
# 
# exp137 <- readDNAStringSet('./outputs/annotations/6461x13150exp137/6461x13150exp137.fna')
# 
# exp137[2] %>% writeXStringSet(filepath = './outputs/exp137_50kb_plas.fna')
# 
# 
# recipexp137 <- readDNAStringSet('./outputs/annotations/13150/13150.fna')
# recipexp137[2]%>% writeXStringSet(filepath = './outputs/RECIPexp137_50kb_plas.fna')
# 
# # PGHKFEJK_01855	EKEIMHEI_00148	LOEINDNJ_01862
# #218	5798	3008	PGHKFEJK_01897	NA	LOEINDNJ_01871 THESE LOC TAGS HAVE V DIF LENS
# 
# 
# # THESE ARE TETO
# # PGHKFEJK_01855
# # EKEIMHEI_00148
# # MIPFMNEC_01862
# 
# 
# gff13150 <- gff_parse3('./outputs/annotations/13150/13150.gff')
# gff6461 <- gff_parse3('./outputs/pan_genomes/6461_13150_6461x13150exp137/6461/6461.gff')
# gff6461x13150exp137 <- gff_parse3('./outputs/annotations/6461x13150exp137/6461x13150exp137.gff')
# 
# 
# # LOEINDNJ_01871 long one
# 
# # PGHKFEJK_01897 short one
# 
# # 6461x13150_1	tet(O)	chromosome
# # 6461x13150_1	blaOXA-605	chromosome
# # 6461x13150_2	tet(O)	plasmid
# # 6461x13150_2	aph(3')-IIIa	plasmid
# # 6461x13150exp137_1	blaOXA-605	chromosome
# # 6461x13150exp137_2	tet(O)	plasmid
# # 6461x13150exp137_2	aph(3')-IIIa	plasmid
# #
# 
# 
# 
# 
# 
# 
# 
# gff_parse3('./outputs/annotations//SampleA.gff')
# 
# read_tsv('./outputs/abricate/SampleA.fna.abricate')
# 
# 
# read_tsv('./outputs/pan_genomes/all_6461s/')
# 
# 
# 
# 
# 
# 
# 6461 + 13150 = 6461x13150exp137