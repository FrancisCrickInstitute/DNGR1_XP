options(stringsAsFactor = F)

library(tidyverse)
library(ggpubr)

#### 1. Import and merge data ####

# Importing cns files
seg_paths<-list.files("SCNAs/segments_cnvkit",full.names = T)
samples<-seg_paths %>% strsplit(.,"/") %>% sapply(.,"[",3) %>% gsub(".cns","",.)

# Group is D/R/W - ie the host genotypes
group<-samples %>% substr(.,3,3)

# Table with the regressor/progressor phenotype per sample 
phenotype_df <- read.delim("annotations/phenotype.txt")

# Vector with order matching the file order
phenotype<-phenotype_df$Phenotype[match(samples,phenotype_df$Cell.line)]

# Read segments into a list
segs_list<-lapply(seg_paths, read.delim)

# Add sample, host and phenotype to each df
for (i in 1:length(samples)){
  segs_list[[i]]$Sample<-samples[i]
  segs_list[[i]]$Group<-group[i]
  segs_list[[i]]$Phenotype<-phenotype[i]
}

# Merge into a single table, filter out X and Y chromosomes and sort
segs_df<-do.call(rbind.data.frame,segs_list) %>% 
  filter(!chromosome%in%c("X","Y")) %>% 
  arrange(Group,Phenotype,chromosome)





#### 2. Plot logR along the genome per sample ####


# chromosome as factor
segs_df$chromosome<-factor(segs_df$chromosome,levels=1:19,ordered = T)

# Vector to transform samples as numbers for plotting
samples_n<-cbind.data.frame(Sample=unique(segs_df$Sample),
                            n=1:length(unique(segs_df$Sample)))

# Add pehnotype
samples_n_phe<-phenotype_df %>% 
  rename(Sample=Cell.line) %>% 
  left_join(samples_n,.) %>% 
  mutate(phenoColour=ifelse(Phenotype=="Progressor","#a50026","#006837"))


# Blank out anything between 0.6 and -0.6, the std threshold for calling CN
segs_capped<-segs_capped %>% 
  mutate(log2R=ifelse(log2R>(-0.6) & log2R<0.6,0,log2R))

segs_capped %>% 
  ggplot(.)+geom_rect(aes(xmin=start,xmax=end,ymin=n,ymax=n1,fill=log2R))+
  scale_fill_gradientn(limits = c(-1.1,1.1),
                       colours=c("navyblue", "white", "darkred"))+
  facet_grid(.~chromosome,space="free_x",scales="free_x")+
  scale_y_continuous(breaks=samples_n$n+0.5,labels=samples_n$Sample,expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        axis.text.y = element_text(colour = samples_n_phe$phenoColour))+
  geom_hline(yintercept = c(15,26))


### MYC seems commonly gained - add a dashed line 

myc_coords<-cbind.data.frame(chromosome="15",
                             mid=61987915) %>% 
  mutate(chromosome=factor(chromosome,levels=1:19,ordered = T))
  
segs_capped %>% 
  ggplot(.)+geom_rect(aes(xmin=start,xmax=end,ymin=n,ymax=n1,fill=log2R))+
  scale_fill_gradientn(limits = c(-1.1,1.1),
                       colours=c("navyblue", "white", "darkred"))+
  facet_grid(.~chromosome,space="free_x",scales="free_x")+
  scale_y_continuous(breaks=samples_n$n+0.5,labels=samples_n$Sample,expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x=unit(0.1, "lines"),
        axis.text.y = element_text(colour = samples_n_phe$phenoColour))+
  geom_hline(yintercept = c(15,26))+
  geom_vline(data=myc_coords,aes(xintercept =mid),linetype="dashed",alpha=0.3)



#### 3. Proportion of the genome altered by group ####

# Call SCNAs
segs_df_scnas<-segs_capped %>% 
  mutate(Call=ifelse(log2>0.6,"Gain",
                     ifelse(log2<(-0.6),"Loss","Neutral")))

## Gains

# Calculate proportion of the genome gained
props_gain<-segs_df_scnas %>% 
  mutate(Length=end-start,
         Altered=ifelse(Call=="Gain","SCNA","Neutral")) %>%
  group_by(Sample) %>% 
  mutate(TotLen=sum(Length), LenProp=Length/TotLen) %>% 
  ungroup() %>%  
  group_by(Sample,Group,Altered) %>% 
  summarise(Proportion= sum(LenProp))

# Add samples with no gains to avoid dropping them
complete_props<-props %>% 
  ungroup %>% 
  expand(Sample,Altered) %>% 
  left_join(.,props)

# Transform NAs into 0 and add host GT column
complete_props$Proportion[is.na(complete_props$Proportion)]<-0  
complete_props$Group<-complete_props$Sample %>% substr(.,3,3)

# Comparisons for ggpubr
my_comp=list(c("D","R"),c("D","W"),c("R","W"))

# Plot proportion of the genome gained
complete_props %>% 
  filter(Altered=="SCNA") %>% 
  ggplot(.,aes(Group,Proportion))+
  geom_boxplot(aes(fill=Group),alpha=0.2,width=0.6, outlier.shape = NA)+
  geom_jitter(height = 0,width = 0.2,size=1.5)+
  theme_bw(base_size = 16)+
  ylab("Proportion of the genome gained")+
  stat_compare_means(comparisons = my_comp)+
  scale_x_discrete(labels=c("DNGR-1 KO","RAG1 KO","WT"))+
  scale_fill_manual(values=c("#4575b4","#d73027","black"))+
  xlab("")+
  theme(legend.position = "none")


## Losses

# Calculate proportion of the genome lost
props_gain<-segs_df_scnas %>% 
  mutate(Length=end-start,
         Altered=ifelse(Call=="Loss","SCNA","Neutral")) %>%
  group_by(Sample) %>% 
  mutate(TotLen=sum(Length), LenProp=Length/TotLen) %>% 
  ungroup() %>%  
  group_by(Sample,Group,Altered) %>% 
  summarise(Proportion= sum(LenProp))

# Add samples with no losses to avoid dropping them
complete_props<-props %>% 
  ungroup %>% 
  expand(Sample,Altered) %>% 
  left_join(.,props)

# Transform NAs into 0 and add host GT column
complete_props$Proportion[is.na(complete_props$Proportion)]<-0  
complete_props$Group<-complete_props$Sample %>% substr(.,3,3)

# Plot proportion of the genome lost
complete_props %>% 
  filter(Altered=="SCNA") %>% 
  ggplot(.,aes(Group,Proportion))+
  geom_boxplot(aes(fill=Group),alpha=0.2,width=0.6, outlier.shape = NA)+
  geom_jitter(height = 0,width = 0.2,size=1.5)+
  theme_bw(base_size = 16)+
  ylab("Proportion of the genome lost")+
  stat_compare_means(comparisons = my_comp)+
  scale_x_discrete(labels=c("DNGR-1 KO","RAG1 KO","WT"))+
  scale_fill_manual(values=c("#4575b4","#d73027","black"))+
  xlab("")+
  theme(legend.position = "none")



#### 4. Plot frequency of SCNAs along the genome by group ####

source("TuCNtools_functions.R") 

# Transform calls to TuCNtools format
segs_calls<-segs_df_scnas %>% 
  mutate(Call=round((2^log2)*2)-2,
         LOH=ifelse(Call%in%c(-1),1,0),
         AI=ifelse(Call/2==round(Call/2),0,1)) %>% 
  select(chromosome,start,end,Call,LOH,AI,Sample) %>% 
  rename(Chromosome=chromosome,Start=start,End=end)

# Perform MCS
MCS_list<-MCS_calls_LOH_AI(segs_calls)


## Per genotype

# Vector with genotype per sample
genotype_vec<-samples_n_phe %>% 
  filter(Sample%in%colnames(MCS_list$Call)) %>% 
  mutate(Genotype=substr(Sample,3,3)) %>% 
  pull(Genotype)

# Calculate CN rates per group
genotype_rates<-getGroupRates(MCS_list,genotype_vec)

# Chromosomes as factors for plotting
genotype_rates$W$Chromosome<-factor(genotype_rates$W$Chromosome,levels=1:19,ordered=T)
genotype_rates$D$Chromosome<-factor(genotype_rates$D$Chromosome,levels=1:19,ordered=T)
genotype_rates$R$Chromosome<-factor(genotype_rates$R$Chromosome,levels=1:19,ordered=T)

# Plot CN freq plots
plot_singleSample_SCNA_LOH_AI(genotype_rates$W,"WT host")
plot_singleSample_SCNA_LOH_AI(genotype_rates$D,"DNGR-1 KO host")
plot_singleSample_SCNA_LOH_AI(genotype_rates$R,"RAG1 KO host")


## Per phenotype

# Vector with phenotype per sample 
phenotype_vec<-samples_n_phe %>% 
  filter(Sample%in%colnames(MCS_list$Call)) %>% 
  pull(Phenotype)

# Calculate CN rates per group
phenotype_rates<-getGroupRates(MCS_list,phenotype_vec)

# Chromosomes as factors for plotting
phenotype_rates$Progressor$Chromosome<-factor(phenotype_rates$Progressor$Chromosome,levels=1:19,ordered=T)
phenotype_rates$Regressor$Chromosome<-factor(phenotype_rates$Regressor$Chromosome,levels=1:19,ordered=T)

# Plot CN freq plots
plot_singleSample_SCNA_LOH_AI(phenotype_rates$Progressor,"Progressors")
plot_singleSample_SCNA_LOH_AI(phenotype_rates$Regressor,"Regressors")

