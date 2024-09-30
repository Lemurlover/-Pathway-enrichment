#Udemy course - Gene Set Enrichment Analysis using ClusterProfiler, fgsea and gProfiler
#25-02-2024
#Author - George Sentis - DAMA the Seq hub

####Module 1 - Intro & Setting up your environment####
##Some basic commands about R
#Check your current directory
getwd()

#Navigate to the defined directory
#setwd('~/Desktop/udemy/pathway_enrichment/')

#Basic variable assignment
x1 <- 1
x2 = 2
x3 <- x2 <- 1 
x3 = x1 = 2

#vectors and data frames
genes <- c('IFN','GAPDH','TP53','FAM30A')
df = data.frame(newVariable = genes,
                newVariable2 = 1,
                newVariable3 = c(1,2),
                newVariable4 = c(1,2,3,4))

#clear any potential object in your environment
rm(x1)
rm(list=ls())

#Installing required libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('clusterProfiler') #done
BiocManager::install('gprofiler2') #DONE
BiocManager::install('fgsea') #DONE
BiocManager::install('msigdbr') #DONE
BiocManager::install('GSEABase') #done
BiocManager::install('org.Hs.eg.db') #done
BiocManager::install('DOSE') #not done
BiocManager::install('enrichplot') #not done
BiocManager::install("ComplexHeatmap") #ok
install.packages('circlize')
install.packages('openxlsx')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('tibble')
install.packages('ggVennDiagram')

#Load required libraries
library(openxlsx)
library(ggplot2)
library(ggVennDiagram)
library(dplyr)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(gprofiler2)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(ComplexHeatmap)
library(tibble)

#basic dplyr syntax
df = data.frame(newVariable = c('IFN','GAPDH','TP53','FAM30A'),
                newVariable2 = 1,
                newVariable3 = c(1,2),
                newVariable4 = c(1,2,3,4))

#equivalent actions
df$counts = 1

mutate(df, counts2 = 1)
df <- mutate(df, counts2 = 1)

df %>%
  mutate(counts3 = 1)

rm(df)

####Pathway enrichment####
#Pathway enrichment is performed mainly using 3 enrichment methods
#Overrepresentation analysis (ORA), Functional Class Scoring (FCS) and Topology-based methods

#/doi.org/10.1371/journal.pcbi.1002375

# ORA:
#   Input: a list of DE genes to determine which sets of DE genes 
#   are over-represented or under-represented.
#     Depends strongly on the criteria used to define significance of DE,
#     (i.e. LogFC & P-value thresholds, statistical test)

# FCS:
#   Input: complete list of genes (not just DE), ranked by
#   a measure (i.e. either expression values of logFC and 
#   statistical significance measure)
#     Reduces reliance on gene selection criteria

# TB:
#   Take into account pathway topology and interaction 
#   relationships between molecules
##reading data and results from GSE205748 DEA
#reading metadata
meta <- read.table('GSE205748_series_matrix_edit.txt',
                   header = TRUE)

#Description of the dataset
View(meta)

##reading expression data
exprcpm <- read.xlsx('GSE205748_cpm.xlsx',
                     rowNames = TRUE)

exprcpm <- read.xlsx('GSE205748_cpm.xlsx')
                     
 
View(exprcpm)

#removing meta and exprcpm for now
rm(meta,exprcpm)

#reading DE results
degs <- read.xlsx('GSE205748_DE_results_PsA_les_vs_healthy.xlsx')
                 
View(degs)


#### Module 2 - setting up ORA gene lists ####

#First thing we do is annotate up- and down-regulation

#we set thresholds
logFC_threshold = 0.58
padj_threshold=0.05

#and annotate
degs <- degs %>%
  mutate(regulation = case_when(
    logFC > logFC_threshold & fdr < padj_threshold ~ 'Upregulated',
    logFC < -logFC_threshold & fdr < padj_threshold ~ 'Downregulated',
    .default = 'Not significant'
  ))

##checking if annotation was performed correctly using a volcano plot
ggplot(degs, aes(x = logFC,
                 y=-log10(fdr),
                 color = regulation)) + 
  geom_point() +
  theme_minimal() + 
  geom_hline(yintercept = -log10(padj_threshold)) + 
  geom_vline(xintercept = logFC_threshold) + 
  geom_vline(xintercept = -logFC_threshold)

#check numbers of DE up- and down-regulated genes
table(degs$regulation)

#isolating
upgenes <- degs %>%
  filter(regulation == 'Upregulated', 
         !is.na(hgnc_symbol),
         !duplicated(hgnc_symbol)) %>%
  pull(hgnc_symbol)

downgenes <- degs %>%
  filter(regulation == 'Downregulated', 
         !is.na(hgnc_symbol),
         !duplicated(hgnc_symbol)) %>%
  pull(hgnc_symbol)

#getting background gene list
background <- degs %>%
  filter(!is.na(hgnc_symbol)) %>%
  filter(!duplicated(hgnc_symbol)) %>%
  pull(hgnc_symbol)


#### Module 3 - Gathering pathway gene sets ####
msigdbr::msigdbr_collections()
msigdbr::msigdbr_species()

genesets1 <- msigdbr(species = 'Homo sapiens',
                    category= 'H') %>%
  bind_rows(msigdbr(species = 'Homo sapiens',
                    category= 'C2',
                    subcategory = 'CP'))%>%
  bind_rows(msigdbr(species = 'Homo sapiens',
                    category= 'C5',
                    subcategory = 'GO:BP'))
  
#getting term sizes
termsizes <- genesets1 %>%
  distinct(gs_name,gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(count = n()) 

#filter out too low genesets
sets_to_keep <- termsizes %>%
  filter(count > 10) %>%
  pull(gs_name)

genesets <- genesets1 %>%
  filter(gs_name %in% sets_to_keep)

####Module 4 - ORA using gProfiler####
#up
ORA_gp_up <- gost(upgenes,
                  correction_method = 'fdr',
                  domain_scope = 'custom_annotated', #annotation refers to the data source (i.e. gene list).
                                                     #filtering out non-annotated genes in the data source is
                                                     #the default behaviour of gProfiler and results in less
                                                     #enrichment of very large "noise" terms
                  sources = c('GO:BP','KEGG','REAC','WP'),
                  custom_bg = background)
res_gp_up <- ORA_gp_up$result
#down
ORA_gp_down <- gost(downgenes,
                    correction_method = 'fdr',
                    domain_scope = 'custom_annotated', #annotation refers to the data source (i.e. gene list).
                                                       #filtering out non-annotated genes in the data source is
                                                       #the default behaviour of gProfiler and results in less
                                                       #enrichment of very large "noise" terms
                    sources = c('GO:BP','KEGG','REAC','WP'),
                    custom_bg = background)
res_gp_down <- ORA_gp_down$result

####Module 5 - ORA using clusterProfiler####
t2g_cp <- genesets %>%
  dplyr::distinct(gs_name,
                  gene_symbol) %>%
  as.data.frame()

#running the enrichment
#up
ORA_cp_up <- enricher(upgenes,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = 'BH',
                      universe = background,
                      qvalueCutoff = 0.2,
                      TERM2GENE = t2g_cp)

res_cp_ora_up <- ORA_cp_up@result %>%
  filter(p.adjust < 0.05)

#down
ORA_cp_down <- enricher(downgenes,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        universe = background,
                        qvalueCutoff = 0.2,
                        TERM2GENE = t2g_cp)

res_cp_ora_down <- ORA_cp_down@result %>%
  filter(p.adjust < 0.05)



####Module 7.5 checking overlap between ORA methods####
#checking ORA methods overlap
cpIDs <- stringr::str_to_lower(gsub('_',' ',gsub('^[[:alpha:]]+_','',res_cp_ora_up$ID)))
gpIDs <- stringr::str_to_lower(gsub('-',' ',gsub(',','',res_gp_up$term_name)))

ggVennDiagram(list(gProfiler = gpIDs,
                   clusterProfiler = cpIDs)) +
  coord_flip()

Enrich_compare <- data.frame(Term = cpIDs,
                             ClustPr = 1) %>%
  full_join(data.frame(Term = gpIDs,
                       gPr = 1), by = c('Term')) %>%
  mutate(Both = ifelse((ClustPr + gPr) == 2, TRUE, FALSE)) %>%
  arrange(Both, Term)


####Module 8 - preparing the ranked gene list####
#create a ranked named vector based on DE results
ranked <- degs %>%
  mutate(rank_metric = logFC*(-log10(fdr))) %>%
  arrange(desc(rank_metric)) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != '') %>%
  select(hgnc_symbol, rank_metric)

##deduplicating and keeping the gene with the highest absolute rank
duplicates <- ranked %>%
  filter(duplicated(hgnc_symbol)) %>%
  pull(hgnc_symbol) %>%
  unique()

for (i in 1:length(duplicates)) {
  i=1
  untangled <- ranked[ranked$hgnc_symbol %in% duplicates[i],] %>%
    arrange(desc(abs(rank_metric)))
  
  if(i == 1){
    drop <- c(rownames(untangled)[2:length(rownames(untangled))])
  } else {
    drop <- c(drop,rownames(untangled)[2:length(rownames(untangled))])
  }
  
}

#creating deduplicated named ranked vector
ranked_dedup <- ranked[!rownames(ranked) %in% drop,]
ranked_vector <- ranked_dedup$rank_metric
names(ranked_vector) <- ranked_dedup$hgnc_symbol

####Module 7 - FCS using GSEA of clusterProfiler####
gsea_cp <- GSEA(geneList = ranked_vector,
                pvalueCutoff = 0.05,
                pAdjustMethod = 'BH',
                eps = 0,
                seed = 149,
                minGSSize = 15,
                maxGSSize = 500,
                TERM2GENE = t2g_cp)


res_fcs_cp <- gsea_cp@result

####Module 8 - FCS using GSEA of fgsea####
msigdbr_list = split(x = genesets$gene_symbol, f = genesets$gs_name) # creating list of genesets

set.seed(149)

FCS_gsea_all <- fgsea(pathways = msigdbr_list,
                      stats = ranked_vector,
                      eps = 0,
                      minSize = 15,
                      maxSize = 500,
                      nproc = 6)


res_fcs_gsea <- FCS_gsea_all %>%
  filter(padj < 0.05)


##checking overlap between similar FCS methods
ggVennDiagram(list(gsea = res_fcs_gsea$pathway,
                   cp = res_fcs_cp$ID)) + 
  coord_flip()

####Module 9 - plotting gProfiler ORA ####
gPr_plot <- res_gp_up %>%
  mutate(regulation = 'Upregulation') %>% 
  bind_rows(res_gp_down %>%
              mutate(regulation = 'Downregulation')) %>%
  filter(term_size > 10) %>%
  mutate(Gene_ratio = intersection_size/query_size) %>% 
  slice_max(order_by = Gene_ratio,n = 30)
  
#plotting
ggplot(gPr_plot,
       aes(x = Gene_ratio,
           y = reorder(term_name,Gene_ratio),
           color = regulation,
           size = -log10(p_value))) + 
  geom_point() + 
  scale_y_discrete(breaks = gPr_plot$term_name,
                   labels = strtrim(gPr_plot$term_name,
                                    width = 30)) + 
  labs(x = 'Gene ratio',
       y = 'Term',
       color = 'Up/Downregulation',
       size = '-Log10(PValue)') + 
  scale_color_manual(breaks = c('Downregulation','Upregulation'),
                     values = c('steelblue1',
                                'firebrick3')) + 
  theme_minimal() + 
  facet_wrap(~ source,
             scales = 'free_y')

####Module 10 - plotting clusterProfiler ORA ####
cPr_plot <- res_cp_ora_up %>%
  mutate(regulation = 'Upregulation') %>% 
  bind_rows(res_cp_ora_down %>%
              mutate(regulation = 'Downregulation')) %>%
  mutate(source = gsub('^([[:alpha:]]+)_[[:print:]]+','\\1',ID)) %>%
  mutate(ID = gsub('_',' ',gsub('^[[:alpha:]]+_','',ID)))

#Converting GeneRatio to numeric
temp <- as.numeric(unlist(strsplit(cPr_plot$GeneRatio, split = '/')))
cPr_plot$GeneRatio2 <- temp[c(TRUE, FALSE)]/temp[c(FALSE, TRUE)]

top_cPr_plot <- cPr_plot %>%
  slice_max(order_by = GeneRatio2, n = 30)

#plotting
pdf('clustPr_ORA_bubble_plot.pdf', width = 7, height = 7)
ggplot(top_cPr_plot,
       aes(x = GeneRatio2,
           y = reorder(ID,GeneRatio2),
           color = regulation,
           size = -log10(p.adjust))) + 
  geom_point() + 
  scale_y_discrete(breaks = cPr_plot$ID,
                   labels = strtrim(cPr_plot$ID,
                                    width = 30)) + 
  labs(x = 'Gene ratio',
       y = 'Term',
       color = 'Up/Downregulation',
       size = '-Log10(Adj. PValue)') + 
  scale_color_manual(breaks = c('Downregulation','Upregulation'),
                     values = c('steelblue1',
                                'firebrick3')) + 
  theme_minimal() + 
  facet_wrap(~ source,
             scales = 'free_y')
dev.off()

####Module 11 - plotting clusterProfiler FCS ####
cP_fcs_plot <- res_fcs_cp %>%
  mutate(source = gsub('^([[:alpha:]]+)_[[:print:]]+','\\1',ID),
         regulation = ifelse(NES > 0, 'Upregulation','Downregulation')) %>%
  mutate(ID = gsub('_',' ',gsub('^[[:alpha:]]+_','',ID))) %>%
  slice_max(order_by = abs(NES), n = 30)
  
#plotting
pdf('clustPr_FCS_bubble_plot.pdf', width = 13, height = 7)
ggplot(cP_fcs_plot,
       aes(x = NES,
           y = reorder(ID,NES),
           color = regulation,
           size = -log10(p.adjust))) + 
  geom_point() + 
  scale_y_discrete(breaks = cPr_plot$ID,
                   labels = strtrim(cPr_plot$ID,
                                    width = 50)) + 
  labs(x = 'NES',
       y = 'Term',
       color = 'Up/Downregulation',
       size = '-Log10(Adj. PValue)') + 
  scale_color_manual(breaks = c('Downregulation','Upregulation'),
                     values = c('steelblue1',
                                'firebrick3')) + 
  theme_minimal() + 
  facet_wrap(~ source,
             scales = 'free_y')
dev.off()

#plotting barplot
pdf('clustPr_FCS_barplot_plot.pdf', width = 12, height = 7)
ggplot(cP_fcs_plot,
       aes(x = NES,
           y = reorder(ID,NES),
           fill = regulation,
           label = p.adjust)) + 
  geom_bar(stat = 'identity') + 
  scale_y_discrete(breaks = cPr_plot$ID,
                   labels = strtrim(cPr_plot$ID,
                                    width = 30)) + 
  labs(x = 'NES',
       y = 'Term',
       fill = 'Up/Downregulation') + 
  scale_fill_manual(breaks = c('Downregulation','Upregulation'),
                     values = c('steelblue1',
                                'firebrick3')) + 
  theme_minimal() + 
  facet_wrap(~ source,
             scales = 'free_y')
dev.off()

####Module 12 - plotting fgsea FCS ####
fgsea_fcs_plot <- res_fcs_gsea %>%
  mutate(source = gsub('^([[:alpha:]]+)_[[:print:]]+','\\1',pathway),
         regulation = ifelse(NES > 0, 'Upregulation','Downregulation')) %>%
  mutate(pathway = gsub('_',' ',gsub('^[[:alpha:]]+_','',pathway))) %>%
  slice_max(order_by = abs(NES), n = 30)

#plotting
ggplot(fgsea_fcs_plot,
       aes(x = NES,
           y = reorder(pathway,NES),
           color = regulation,
           size = -log10(padj))) + 
  geom_point() + 
  scale_y_discrete(breaks = fgsea_fcs_plot$pathway,
                   labels = strtrim(fgsea_fcs_plot$pathway,
                                    width = 30)) + 
  labs(x = 'NES',
       y = 'Term',
       color = 'Up/Downregulation',
       size = '-Log10(Adj. PValue)') + 
  scale_color_manual(breaks = c('Downregulation','Upregulation'),
                     values = c('steelblue1',
                                'firebrick3')) + 
  theme_minimal() + 
  facet_wrap(~ source,
             scales = 'free_y')


#Enrichment running score plot
plotEnrichment(msigdbr_list[['GOBP_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA']], 
               stats = ranked_vector)

#GSEA table
plotGseaTable(msigdbr_list[grepl('INTERFERON_GAMMA',names(msigdbr_list))],
              stats = ranked_vector,
              FCS_gsea_all,
              pathwayLabelStyle = list(size = 8))

####Module 12a - plotting heatmap fgsea FCS####
genes_to_plot <- unlist(FCS_gsea_all$leadingEdge[FCS_gsea_all$pathway == 'GOBP_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA'])

##reading data and results from GSE205748 DEA
#reading metadata
meta <- read.table('GSE205748_series_matrix_edit.txt',
                   header = TRUE)

##reading expression data
exprcpm <- read.xlsx('GSE205748_cpm.xlsx',
                     rowNames = TRUE)

exprcpm$ensembl_gene_id <- rownames(exprcpm)

#getting table for renaming expression to gene symbols
transition_table <- degs %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol)

#isolating leading edge gene expression
exprisol <- data.frame(hgnc_symbol = genes_to_plot) %>%
  left_join(transition_table, by = 'hgnc_symbol') %>%
  left_join(exprcpm, by = 'ensembl_gene_id') %>%
  dplyr:: select(-ensembl_gene_id) %>%
  column_to_rownames('hgnc_symbol')

#z-scoring
exprisol <- t(scale(t(exprisol))) # The scale function scales across columns. We want to scale the genes, therefore we transpose, scale and retranspose to the original layout.

#checking if order of samples is the same in metadata and counts
identical(meta$Sample_code, colnames(exprisol))

#Heatmap
pdf('Heatmap_IFN_g_leading_edge_genes.pdf', width = 8, height = 10)
Heatmap(exprisol, 
        heatmap_legend_param = list(title = 'Expression\nZ-score'),
        column_title = 'IFN-g Leading Edge genes',
        row_names_gp = gpar(fontsize = 9),
        top_annotation = HeatmapAnnotation(Tissue_type = meta$Tissue_type,
                                           which = 'column',
                                           col = list(Tissue_type = c('Healthy' = 'chartreuse3',
                                                                      'PsA_les' = 'sienna3',
                                                                      'PsA_uninv' = 'honeydew2'))
                                           ))
dev.off()

####Final module - Reproducibility####
save.image('Enrichment_results_Udemy.RData')

sink('Enrichment_results_Udemy.txt')
sessionInfo()
sink()


# ORA:
#   Input: a list of DE genes to determine which sets of DE genes 
#   are over-represented or under-represented.
#     Depends strongly on the criteria used to define significance of DE,
#     (i.e. LogFC & P-value thresholds, statistical test)

# FCS:
#   Input: complete list of genes (not just DE), ranked by
#   a measure (i.e. either expression values of logFC and 
#   statistical significance measure)
#     Reduces reliance on gene selection criteria

# TB:
#   Take into account pathway topology and interaction 
#   relationships between molecules