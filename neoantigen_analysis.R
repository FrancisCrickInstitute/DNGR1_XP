options(stringsAsFactor = F)

library(tidyverse)
library(ggpubr)
library(cowplot)
library(ComplexHeatmap)

#### Load filtered MHC I files ####

paths<-list.files("neoantigens/results",full.names = T)

files_list<-list()

for (i in 1:length(paths)){
  files_list[[i]]<-read.delim(paths[i])
  name<-paths[i] %>% strsplit(.,"/") %>% sapply(.,"[",3) %>% gsub(".tsv","",.)
  files_list[[i]]$Sample.Name<-name
}

raw_df<-do.call(rbind.data.frame,files_list)

phenotype<-read.delim("annotations/phenotype.txt")

FABP_ensembl<-read.delim("annotations/FABP_ensembl.txt",header = F)$V1 %>% as.character()
MBP_genes<-read.delim("annotations/MBP_mouse_genes.txt",header = T)$Mouse_gene_symbol %>% as.character()

# Collect all data 

all_neoags<-raw_df %>% 
  mutate(Cell.line=strsplit(Sample.Name,"_") %>% sapply("[",1)) %>% 
  left_join(.,phenotype) %>% 
  mutate(Genotype=substr(Sample.Name,3,3) %>% str_replace_all(c("^W$" = "WT", "^D$" = "DNGR-1 KO", "^R$" = "Rag1 KO"))) %>% 
  mutate(BA_PASS=ifelse(Best.MT.IC50.Score<=150,T,F),
         FABP=ifelse(Ensembl.Gene.ID%in%FABP_ensembl,"FABP","not FABP"),
         MBP=ifelse(Gene.Name%in%MBP_genes,"MBP","not MBP"))



#### Progressor vs Regressor - n of neoags ####  

##### All neoantigens ##### 

counts_all<-all_neoags %>% 
  filter(BA_PASS==T) %>% 
  count(Sample.Name,Phenotype,Genotype)

ggplot(counts_all,aes(Phenotype,n,color=Phenotype))+
  geom_boxplot(aes(fill=Phenotype),alpha=0.2,width=0.6, outlier.shape = NA)+
  geom_jitter(height = 0,width = 0.2,size=2.5)+
  scale_color_manual(values=c("Progressor"="#A32A31","Regressor"="#3658A7"))+
  scale_fill_manual(values=c("Progressor"="#A32A31","Regressor"="#3658A7"))+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")+
  ylab("Predicted neoantigens (BA<=150)")+
  xlab("")+
  stat_compare_means(method.args = list(alternative = "two.sided"), aes(label = sprintf("p = %1.g", as.numeric(..p.format..))))


# Empty counts to fill in samples removed by other filters and avoid dropping zeros
empty<-counts_all %>% mutate(n=0)


#### WT vs DNGR1 KO vs Rag1 KO - n of neoags ####  

my_comparisons<-list(c("DNGR-1 KO","WT"),c("Rag1 KO","WT"),c("DNGR-1 KO","Rag1 KO"))


counts_all$Genotype<-factor(counts_all$Genotype,
                                levels=c("WT","DNGR-1 KO","Rag1 KO"),
                                ordered=T)

ggplot(counts_all,aes(Genotype,n,color=Genotype))+
  geom_boxplot(aes(fill=Genotype),alpha=0.2,width=0.6, outlier.shape = NA)+
  geom_jitter(height = 0,width = 0.2,size=2.5)+
  scale_color_manual(values=c("DNGR-1 KO"="#4575b4","Rag1 KO"="#d73027","WT"="black"))+
  scale_fill_manual(values=c("DNGR-1 KO"="#4575b4","Rag1 KO"="#d73027","WT"="black"))+
  theme_bw(base_size = 14)+
  ylab("Neoantigens on all genes (BA<=150)")+
  xlab("")+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons)


##### FABP #####

counts_FABP<-all_neoags %>% 
  filter(BA_PASS==T, FABP=="FABP") %>% 
  count(Sample.Name,Phenotype,Genotype)

if(nrow(counts_FABP)<nrow(phenotype)){
  missing_rows<-filter(empty,!Sample.Name%in%counts_FABP$Sample.Name)
  counts_FABP<-rbind.data.frame(counts_FABP,missing_rows)
}

counts_FABP$Genotype<-factor(counts_FABP$Genotype,
                                 levels=c("WT","DNGR-1 KO","Rag1 KO"),
                                 ordered=T)

ggplot(counts_FABP,aes(Genotype,n,color=Genotype))+
  geom_boxplot(aes(fill=Genotype),alpha=0.2,width=0.6, outlier.shape = NA)+
  geom_jitter(height = 0,width = 0.2,size=2.5)+
  scale_color_manual(values=c("DNGR-1 KO"="#4575b4","Rag1 KO"="#d73027","WT"="black"))+
  scale_fill_manual(values=c("DNGR-1 KO"="#4575b4","Rag1 KO"="#d73027","WT"="black"))+
  theme_bw(base_size = 14)+
  ylab("Neoantigens on FABP genes (BA<=150)")+
  xlab("")+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons)


##### MBP #####

counts_MBP<-all_neoags %>% 
  filter(BA_PASS==T, MBP=="MBP") %>% 
  count(Sample.Name,Phenotype,Genotype)

if(nrow(counts_MBP)<nrow(phenotype)){
  missing_rows<-filter(empty,!Sample.Name%in%counts_MBP$Sample.Name)
  counts_MBP<-rbind.data.frame(counts_MBP,missing_rows)
}

counts_MBP$Genotype<-factor(counts_MBP$Genotype,
                                levels=c("WT","DNGR-1 KO","Rag1 KO"),
                                ordered=T)


ggplot(counts_MBP,aes(Genotype,n,color=Genotype))+
  geom_boxplot(aes(fill=Genotype),alpha=0.2,width=0.6, outlier.shape = NA)+
  geom_jitter(height = 0,width = 0.2,size=2.5)+
  scale_color_manual(values=c("DNGR-1 KO"="#4575b4","Rag1 KO"="#d73027","WT"="black"))+
  scale_fill_manual(values=c("DNGR-1 KO"="#4575b4","Rag1 KO"="#d73027","WT"="black"))+
  theme_bw(base_size = 14)+
  ylab("Neoantigens on MBP genes (BA<=150)")+
  xlab("")+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = my_comparisons)



#### Heatmaps of FABP genes ####


FABP_n_table<-all_neoags %>% 
  filter(FABP=="FABP") %>% 
  select(Cell.line,Phenotype,Genotype,Gene.Name,BA_PASS) %>% 
  mutate(isCalled=1,BA_PASS_TF=ifelse(BA_PASS,1,0)) %>% 
  group_by(Cell.line,Phenotype,Genotype,Gene.Name) %>% 
  mutate(nCalled=sum(isCalled),n_BA_PASS=sum(BA_PASS_TF)) %>% 
  select(Cell.line,Phenotype,Genotype,Gene.Name,nCalled,n_BA_PASS) %>% 
  ungroup()

FABP_n_table_PASS<-FABP_n_table %>% 
  filter(n_BA_PASS>0)


samples<-unique(all_neoags$Cell.line) 
missing_samples_b<-samples[!samples%in%FABP_n_table_PASS$Cell.line]

missing_df_b<-empty %>% 
  mutate(Cell.line=Sample.Name %>% strsplit(.,"_") %>% sapply(., "[",1) ) %>% 
  filter(Cell.line %in% missing_samples_b) %>% 
  select(Cell.line,Phenotype,Genotype) %>% 
  mutate(Gene.Name=0,nCalled=0,n_BA_PASS=0)

# Order
FABP_n_table_ordered<-rbind.data.frame(FABP_n_table_PASS,missing_df_b) %>% 
  # order here
  mutate(ord1=ifelse(Genotype=="WT",1,ifelse(Genotype=="DNGR-1 KO",10,100)),
         ord2=ifelse(Phenotype=="Progressor",1,2)) %>%
  arrange(ord1+ord2,Cell.line) %>%
  select(-ord1,-ord2) %>%
  unique %>%
  ungroup()


## Make wide

wide_FABP<-FABP_n_table_ordered %>% 
  select(Cell.line,n_BA_PASS,Gene.Name) %>% 
  pivot_wider( names_from = Cell.line, values_from=n_BA_PASS,values_fill=0) %>% 
  filter(Gene.Name!=0) %>% 
  # col to rownames
  column_to_rownames("Gene.Name")


# For annotations

anno_df<-FABP_n_table_ordered %>% 
  select(Cell.line,Phenotype,Genotype) %>% 
  unique

vector_anno_top<-HeatmapAnnotation(
  #Vectors of variables to annotate
  Phenotype=anno_df$Phenotype,
  #List of named vectors for colours
  col = list(Phenotype = c("Progressor"="#A32A31","Regressor"="#3658A7")),
  # #Add white lines in between cells
  gp = gpar(col = "white", lwd = 0.01),
  #Legend font sizes
  annotation_legend_param = list(labels_gp = gpar(fontsize = 12), 
                                 title_gp = gpar(fontsize = 14)))

library(circlize)
col_fun = colorRamp2(c(0, 1, 2), c("white", "grey", "black"))



#### Add barplots with totals per gene on the right, row anno

#Create a HeatmapAnnotation object with the anno_barplot function 
barplot_anno_right<-HeatmapAnnotation(n = anno_barplot(cbind.data.frame(Total=wide_FABP %>% rowSums),
                                                       #Assign colors with a named vector
                                                       gp = gpar(fill = "darkgrey"),
                                                       bar_width = 0.8),
                                      annotation_legend_param = list(labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 8)),
                                      #Row annotation
                                      which="row")


FABP_wide_heatmap<-Heatmap(wide_FABP,
                       col=col_fun,
                       name="neoags",
                       #Adds white lines in between matrix cells
                       rect_gp = gpar(col = "grey", lwd = 0.3),
                       #Add annotations
                       top_annotation = vector_anno_top,
                       # right_annotation = barplot_anno_right,
                       cluster_rows = F,
                       cluster_columns = F,
                       #Split columns by variable
                       column_split = factor(anno_df$Genotype,levels=unique(anno_df$Genotype),ordered=T),
                       column_gap=unit(5,"mm"),
                       #Size of column names
                       column_names_gp = gpar(fontsize = 10),
                       #Size of row names
                       row_names_gp = gpar(fontsize = 12),
                       row_names_side = "left",
                       #Column titles (in this case, Patient1 etc)
                       column_title_gp = gpar(size=1))



#### New heatmap with only regressors and order alphabetically ####


FABP_n_table_ordered_regressors<-FABP_n_table_ordered %>% 
  filter(Phenotype=="Regressor") %>% 
  mutate(ncellline=Cell.line %>% substr(.,4,5) %>%  as.numeric) %>% 
  arrange(ncellline) %>% 
  mutate(ord1=ifelse(Genotype=="WT",1,
                     ifelse(Genotype=="DNGR-1 KO",10,100)),
         ord2=ifelse(Phenotype=="Progressor",1,2)) %>%
  arrange(ord1+ord2) %>% 
  select(-c(ord1,ord2,ncellline))

## Make wide

wide_FABP_R<-FABP_n_table_ordered_regressors %>% 
  select(Cell.line,n_BA_PASS,Gene.Name) %>% 
  pivot_wider( names_from = Cell.line, values_from=n_BA_PASS,values_fill=0) %>% 
  filter(Gene.Name!=0) %>% 
  # col to rownames
  column_to_rownames("Gene.Name") %>% 
  arrange(rownames(.))


# For annotations

anno_df_R<-FABP_n_table_ordered_regressors %>% 
  select(Cell.line,Genotype) %>% 
  unique

vector_anno_top_R<-HeatmapAnnotation(
  #Vectors of variables to annotate
  Genotype=anno_df_R$Genotype,
  #List of named vectors for colours
  col = list(Genotype = c("DNGR-1 KO"="#4575b4","Rag1 KO"="#d73027","WT"="black")),
  # #Add white lines in between cells
  gp = gpar(col = "white", lwd = 0.01),
  #Legend font sizes
  annotation_legend_param = list(labels_gp = gpar(fontsize = 12), 
                                 title_gp = gpar(fontsize = 14)))



library(circlize)
col_fun = colorRamp2(c(0, 1, 2), c("white", "grey", "black"))


##### With space #####

final_heatmap<-Heatmap(wide_FABP_R,
                       col=col_fun,
                       name="neoags",
                       #Adds white lines in between matrix cells
                       rect_gp = gpar(col = "grey", lwd = 0.3),
                       #Add annotations
                       top_annotation = vector_anno_top_R,
                       # right_annotation = barplot_anno_right,
                       cluster_rows = F,
                       cluster_columns = F,
                       #Split columns by variable
                       column_split = factor(anno_df_R$Genotype,levels=unique(anno_df_R$Genotype),ordered=T),
                       column_gap=unit(2,"mm"),
                       #Size of column names
                       column_names_gp = gpar(fontsize = 10),
                       #Size of row names
                       row_names_gp = gpar(fontsize = 12),
                       row_names_side = "left",
                       #Column titles (in this case, Patient1 etc)
                       column_title_gp = gpar(size=1), 
                       border = T)


