#!/usr/bin/env Rscript

library(tidyverse, quietly = TRUE)
library(formattable, quietly = TRUE)
library(AnnotationForge, quietly = TRUE)
library(seqinr, quietly = TRUE)
library(clusterProfiler, quietly = TRUE)
library(ggsci)


# read emapper result
library(readr)
emapper <- query_seqs_fa_emapper <- read_delim("~/2.compara_analysis_2021.5.20/5.clusterprofiler/0.empper_online/1.increase_gene.eggmapper/query_seqs.fa.emapper.annotations", 
                                               "\t", escape_double = FALSE, col_names = FALSE, 
                                               comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                Gene_Symbol = X6, 
                GO = X7, 
                KO = X9, 
                Pathway = X10, 
                OG = X21, 
                Gene_Name = X22)

  
  
  # Pathway statistics and plot ---------------------------------------------
  gene2pathway <- dplyr::select(emapper, GID, Pathway) %>%
    separate_rows(Pathway, sep = ',', convert = F) %>%
    filter(str_detect(Pathway, 'ko'))
  
  load("~/2.compara_analysis_2021.5.20/5.clusterprofiler/3.genefamily_anno/kegg_info.RData")
  
  gene2pathway <- gene2pathway %>%
    left_join(pathway2name, by = "Pathway") %>%
    dplyr::select(GID, Pathway, Pathway_Name, Pathway_Class, Pathway_Subclass) %>%
    distinct() %>%
    na.omit()
  
  # write.table(gene2pathway, file = "pathway.txt", sep = "\t", quote = F)
  
  # library(readr)
  # gene2pathway <- read_delim("pathway1.txt", "\t", 
  #                       escape_double = FALSE, trim_ws = TRUE)

  
  pathway_anno = length(unique(gene2pathway$GID))
  pathway_stat <- dplyr::select(gene2pathway, GID, Pathway_Class, Pathway_Subclass) %>% 
    distinct() %>% 
    group_by(Pathway_Class, Pathway_Subclass) %>%
    summarise(Count = n(), Percentage = percent(n()/pathway_anno))
  
  pathway_stat$Pathway_Subclass <- ordered(pathway_stat$Pathway_Subclass, levels = pathway_stat$Pathway_Subclass) 
  
 # 按照Count从高到低排序
 library(forcats)
 h <- filter(pathway_stat, Pathway_Class != "Human Diseases") %>%
  arrange(desc(Count)) 
 
  h_arrange <- mutate(h, Pathway_Subclass = factor(Pathway_Subclass, levels = h$Pathway_Subclass)) %>%
  filter(Pathway_Subclass != 'Cellular community - prokaryotes')
  
  ggplot(h_arrange, aes(x=Pathway_Subclass, y = Count)) +
    geom_col(aes(fill=Pathway_Class),
             width=0.6) +
    geom_text(aes(label=Count), nudge_y = 3.3) +
    # scale_y_continuous(labels=percent) + 
    scale_fill_npg() +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45,
                                      hjust=1,
                                      vjust = 1),
          legend.position = c(0.9, 0.6))
          #legend.direction = "horizontal")

  # 张老师的图
 ggplot(h_arrange, aes(x = Pathway_Subclass, y = Percentage)) +
    geom_bar(aes(fill = Pathway_Class), stat = 'identity') +
    geom_text(aes(label = Count), nudge_y = 0.005) +
    scale_y_continuous(labels=percent) + 
    labs(y = "Percent of genes(%)", x ="", fill = "Class") +
    coord_flip() +
    theme_classic()
  
  ggsave("pathway.pdf", p, width = 20, height = 7)
  
  write.table(gene2pathway, file = "pathway.txt", sep = "\t", quote = F)
  write.table(pathway2name, file = 'pathway_name.txt', sep = '\t', quote = F)
  write.table(pathway_stat, file = "pathway_stat.txt", sep = "\t", quote = F, row.names = F)
  
  
  # COG statistics and plot -------------------------------------------------
  library(readr)
  cog_funclass <- read_delim(paste(script_dir, "cog_funclass.tab", sep = "/"), 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
  
  insert_comma <- function(x){
    str_c(x, sep = '', collapse = ',')
  }
  
  gene2cog <- dplyr::select(emapper, GID, COG) %>%
    filter(!is.na(COG)) %>%
    mutate(COG = sapply(str_split(COG, ''), insert_comma)) %>%
    separate_rows(COG, sep = ',', convert = F) %>%
    left_join(cog_funclass, by = c('COG' = 'COG'))
  
  cog_anno = length(unique(gene2cog$GID))
  
  write.table(gene2cog, file = "cog.txt", sep = "\t", quote = F, row.names = F)
  
  p <- ggplot(data = gene2cog) + 
    geom_bar(aes(x = COG, 
                 fill = COG_Name)) +
    labs(title = "COG/KOG Function Classification ", 
         x = "",
         y = "Number of genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.key.size=unit(1,"line"),
          legend.text = element_text(size = 7.5)) +
    guides(fill=guide_legend(ncol=1))
  ggsave("cog.pdf", p, width = 16, height = 7)
  
  # number and percentage ---------------------------------------------------
  anno_stat <- tibble(
    database = c("EggNOG", "GO", "COG/KOG", "KEGG Pathway"),
    number = comma(c(eggnog_anno, go_anno, cog_anno, pathway_anno), digits = 0),
    percentage = percent(c(eggnog_anno, go_anno, cog_anno, pathway_anno)/total_gene)
  )
  
  write.table(anno_stat, "anno_stat.txt", quote = F, row.names = F, sep = "\t")


# set script dir
script_dir <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

# get total gene number
all_gene <- getName.list(read.fasta(file = argv$proteins, 
                                    seqtype = 'AA'))
total_gene = length(all_gene)
emapper <- read_emapper(argv$emapper_anno)

# gene name
gene_info <- dplyr::select(emapper,  GID, Gene_Name) %>%
  dplyr::filter(!is.na(Gene_Name))
eggnog_anno = length(gene_info$GID)

# gene to gene ontology
gene2go <- dplyr::select(emapper, GID, GO) %>%
  separate_rows(GO, sep = ',', convert = F) %>%
  filter(!is.na(GO)) %>%
  mutate(EVIDENCE = 'IEA') 


