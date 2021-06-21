library(tidyverse)
library(clusterProfiler)
library(enrichplot)

library('clusterProfiler')

sessionInfo()

pkgbuild::build('./org.Achi.eg.db/', dest_path = "./Enrichment-clusterprofiler")

# 创建一个文件夹
dir.create('R_Library', recursive = T)

# 将包安装在该文件夹下
install.packages('./org.My.eg.db_1.0.tar.gz', 
                 repos = NULL, #从本地安装
                 lib = 'R_Library') # 安装文件夹

# 加载 OrgDB
library(org.My.eg.db, lib = 'R_Library')

# 基因列表
library(readr)
increase_gene <- read_table2("increase.geneID", 
                             col_names = FALSE)

##提取一列
gene <- pull(increase_gene, X1)

## 富集分析
de_ego <- enrichGO(gene = gene,
                   OrgDb = org.My.eg.db,
                   keyType = 'GID', #重要，geneID输入和orgdb一致
                   ont = 'BP', #GO分子功能富集MF，CC
                   qvalueCutoff = 0.05,
                   pvalueCutoff = 0.05)

de_ego_df <- as.data.frame(de_ego)
head(de_ego_df)

## 可视化
# barplot(de_ego, showCategory = 20)
enrichplot::dotplot(v, showCategory = 20)

emapplot(de_ego, showCategory = 10)

x1 <- pairwise_termsim(de_ego)
emapplot(x1, showCategory = 5)

## 对GO富集结果进行优化
library(DOSE)

de_ego_df$GeneRatio  <- parse_ratio(de_ego_df$GeneRatio) 
ego_arrange <-  mutate(de_ego_df, Description=fct_reorder(Description, GeneRatio, .desc = T)) %>%
  arrange(desc(de_ego_df$p.adjust)) %>%
  head(32)
 

#### 修改合并富集内容
t <- mutate(ego_arrange, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
m <- dplyr::select(t, ID, Description, GeneRatio, p.adjust, Count, richFactor) 
####删除某一行
library(stringr)

v <- filter(t, !(Description %in% c('detection of mechanical stimulus involved in sensory perception',
                                    'sodium ion transmembrane transport',
                                    'detection of stimulus involved in sensory perception of pain',
                                    'larval locomotory behavior',
                                    'sensory perception of mechanical stimulus',
                                    'nucleotide metabolic process',
                                    'carbohydrate derivative metabolic process',
                                    'ibose phosphate metabolic process',
                                    'nucleobase-containing small molecule metabolic process',
                                    'positive regulation of cellular metabolic process'
                                    )))
                                    

### 对优化结果做图
library(ggsci)
library(ggplot2)
library(cowplot)
v <- head(v ,20) 
 ggplot(v, aes(x=Description, y=GeneRatio)) +
        geom_point(shape=21,
                   alpha = 0.5,
                   aes(size=Count, fill = p.adjust)) +
        scale_size(range = c(3,8)) +
        ggtitle("GO Term Cluster") +
        scale_fill_gsea() +
        scale_y_continuous(breaks = seq(0,0.5,0.1),
                           labels = c('0','0.1','0.2','0.3','0.4','0.5')) +
        #scale_fill_gradient2(low='red', high='blue',
        #                   midpoint = 4e-05) +
        theme_classic() +
        theme(axis.text.x = element_text(colour = "black",
                                         angle = 45,
                                       size = 12, hjust =1, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = 12, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = 12),
            legend.position = c(0.8, 0.8))  
    

# 用clusterprofiler帮助文档里的代码
 
library(forcats)
library(enrichplot)
 
 ggplot(v, showCategory = 1, 
        aes(richFactor, fct_reorder(Description, richFactor))) + 
   geom_segment(aes(xend=0, yend = Description)) +
   geom_point(aes(color=p.adjust, size = Count)) +
   scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
   scale_size_continuous(range=c(2, 8)) +
   theme_minimal() + 
   xlim(0,0.225) +
   xlab("Rich Factor") +
   ylab(NULL) + 
   ggtitle("Function Enrichment of GO") 

# pathway
## 如果研究的物种再KEGG数据库中包含，用enrichKEGG做富集，否则，用enricher()

### 准备TERM2GENE
library(readr)
emapper <- query_seqs_fa_emapper <- read_delim("~/2.compara_analysis_2021.5.20/5.clusterprofiler/0.empper_online/query_seqs.fa.emapper.annotations", 
                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                    comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                Gene_Symbol = X6, 
                GO = X7, 
                KO = X9, 
                Pathway = X10, 
                OG = X21, 
                Gene_Name = X22)

pathway2gene <- dplyr::select(emapper, Pathway, GID) %>%
  separate_rows(Pathway, sep=',', convert = F) %>%
  filter(str_detect(Pathway,'ko')) %>%
  mutate(Pathway = str_remove(Pathway, 'ko'))
 
### 准备TERM2NAME,得到pathway与KEGG ID对应关系
library(magrittr)
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
  colnames(keggpathid2name.df) <- c("path_id","path_name")
  return(keggpathid2name.df)
}
pathway2name <- get_path2name()

### 富集分析
library(clusterProfiler)
de_ekp <- enricher(gene,
                   TERM2GENE = pathway2gene,
                   TERM2NAME = pathway2name,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
de_ekp_df <- as.data.frame(de_ekp)
head(de_ekp_df)

### 可视化
barplot(de_ekp, showCategory = 10)
enrichplot::dotplot(de_ekp, showCategory = 20)

### 用clusterprofiler给的代码作图
library(forcats)
library(enrichplot)
library(DOSE)

ekp_arrange <-  mutate(de_ekp_df, Description=fct_reorder(Description, GeneRatio, .desc = T)) %>%
  arrange(desc(de_ekp_df$p.adjust)) %>%
  head(20) %>%
  filter(!is.na(Description))
ekp_arrange$GeneRatio  <- parse_ratio(ekp_arrange$GeneRatio) 
l <- mutate(ekp_arrange, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
view(l)

n <- filter(l, !(Description %in% c('Chemical carcinogenesis - DNA adducts',
                                    'Porphyrin and chlorophyll metabolism',
                                    'Alzheimer disease',
                                    'Biofilm formation - Vibrio cholerae',
                                    'Microbial metabolism in diverse environments',
                                    'Carbon fixation in photosynthetic organisms')))

ggplot(n, showCategory = 1, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 8)) +
  theme_minimal() + 
  xlim(0,0.3) +
  xlab("Rich Factor") +
  ylab(NULL) + 
  ggtitle("Function Enrichment of KEGG")  
  



