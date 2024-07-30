## R scripts
# Scripts for manuscript Smit et al., 2024. Plasma Proteomes Of Acute Myeloid Leukemia Patients treated with Transfusions Reveal Signatures Of Inflammation And Hemostatic Dysregulation
### load necessary packages
library(tidyverse)
library(eulerr)
library(Hmisc)
library(ggrepel)
library(limma)
library(ComplexHeatmap)
library(corrplot)
library(data.table)
library(RColorBrewer)
library(TissueEnrich)
library(ggpubr)


##########################################################
###### PREPROCESS AND IMPUTING OF PROTEOMIC DATA #########
##########################################################

### Load data file from DIA-NN filtered for minimum of 2 precursors, anonymyzed and raw LFQ data was log transformed to obtain normally distributed data
df_prot_raw <- read.csv("../Manuscript_scripts_treatment_Window/Scripts_updateFeb24/Supplement_tables_Feb24/df_raw_manuscript.csv") # load file from Supplemental table 1B
  
# Prepare a cleaner object with protein and sample data
df_prot_clean <- df_prot_raw %>%
  dplyr::select(-Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description) %>%
  dplyr::select(UID, everything())

# create protein dictionary based on protein data file
df_uid_information <- df_prot_raw %>%
  dplyr::select(UID, Protein.Group, Protein.Ids, Protein.Names, Genes, First.Protein.Description)

### Load sample data
#load sample data
pheno <- read.csv("../Manuscript_scripts_treatment_Window/Scripts_updateFeb24/Supplement_tables_Feb24/pheno_manuscript.csv") # load file from Supplemental table 1A

### set a file without QC samples for protein calculations in study groups
pheno_without_qc <- pheno %>% 
  filter(!grepl(Samples, pattern = "pool")) %>%
  mutate(Group = as.character(Group)) %>%
  mutate(Group_fac = as.factor(Group)) %>%
  mutate(DIAGNOSIS = factor(Group, levels = unique(Group))) %>%
  mutate(Date_Trial = as.numeric(str_extract(Samples, "(?<=D)\\d+")))


#Calculate valid values across proteins detected, proteins considered quantified only in >40% of samples per condition (controls or AML cohort)
vvs = df_prot_clean %>% 
  column_to_rownames('UID') %>% 
  dplyr::select(pheno_without_qc$Samples) %>%
  mutate(across(.cols = everything(), ~ !is.na(.))) %>% 
  t() %>% 
  as.data.frame() %>% 
  add_column(ct = pheno_without_qc$DIAGNOSIS) %>% 
  group_by(ct) %>% 
  summarise_all(sum) %>% 
  column_to_rownames('ct') %>%
  mutate(across(.cols = everything(), ~ . / as.numeric(table(list(pheno_without_qc$DIAGNOSIS))))*100) %>% 
  t() %>%
  as.data.frame() %>% 
  mutate(No_pass = rowSums(across(.cols = everything(), ~ . >= 40 ))) %>% # valid value threshold 
  rownames_to_column('UID')

#filter data based on vvs calculated
#IALO # percent
df_prot_valid = df_prot_clean %>% 
  dplyr::select(UID, pheno_without_qc$Samples) %>% 
  subset(UID %in% vvs$UID[vvs$No_pass >=1])


#### Imputation of missing values in protein data
# define parameters
mann_downshift = 1.8
mann_width = 0.3

# create imputed values
imputed_values = df_prot_valid %>% 
  dplyr::select(pheno_without_qc$Samples[order(pheno_without_qc$ID)]) %>% 
  mutate(across(.cols = everything(), ~ rnorm(length(.), (mean(na.omit(.) -  mann_downshift*sd(na.omit(.)))), (sd(na.omit(.))*mann_width))))

# replace NA's in DF with imputed values using coalesce (no need for a mask in tidyverse)  
df_prot_imp = df_prot_valid %>% 
  dplyr::select(colnames(imputed_values)) %>%
  map2_df(., imputed_values, coalesce) %>% 
  bind_cols(UID = df_prot_valid$UID) 


######################################################################
######## OVERVIEW OF PLASMA PROTEOMIC PROFILES BETWEEN GROUPS ########
######################################################################

## Figure 2A Coefficient of variations from QC samples
df_prot_clean %>%
  dplyr::select(UID, pheno$Samples) %>%
  tidyr::gather("Samples", "Intensity", -UID) %>%
  mutate(Intensity = 2^(Intensity)) %>%
  left_join(pheno, by = "Samples") %>%
  dplyr::group_by(UID, Group) %>%
  dplyr::summarise(mean = mean(Intensity, na.rm = T),
                   sd =  sd(Intensity, na.rm =T),
                   coefv = 100* sd/mean) %>%
  ungroup() %>% 
  left_join(df_uid_information, "UID") %>%  
  filter(., Group == "QC") %>% 
  ggplot(.,  aes(x = coefv)) + 
  geom_density(aes(y=..density.. *488*3), lwd = 0.3, colour = "grey45",
               fill =24, alpha = 0.25) +  # 483 proteins quantified in the QC
  geom_histogram(aes(y = ..count..), colour="grey45", fill= "grey85" , alpha=0.5, binwidth=3) +
  geom_vline(xintercept = 32.5, linetype="dashed", colour="orangered") +
  xlab("Coefficient of variation (%)") +
  ylab("Frequency") +
  scale_x_continuous(breaks = seq(0,85,20), limits = c(0,85), expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x= element_text(size=8),
        axis.text.y= element_text(size=8),
        axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),
        plot.title = element_text(size=10)) 

## Fig 2B Ranking of the abundance of all quantified proteins per group, AML or controls
dynamic_range <- df_prot_valid %>%
  dplyr::select(pheno_without_qc$Samples, UID) %>% 
  gather("Samples", "Intensity", -UID) %>%
  filter(Intensity != 0) %>%
  filter(is.finite(Intensity)) %>%
  left_join(df_uid_information, "UID") %>%
  left_join(pheno_without_qc, "Samples") %>% 
  dplyr::select(Genes, Intensity, Protein.Group, UID, Group) %>% 
  group_by(UID, Genes, Protein.Group, Group) %>%
  summarise(Med_Intensity = median(Intensity)) %>% ### median intensity
  ungroup() %>%
  left_join(vvs %>% 
              dplyr::select(-No_pass) %>% 
              gather("Group", "VVSValue", -UID), NULL) %>% 
  filter(VVSValue >= 40) %>% # filter out not valid with a value
  group_by(Group) %>%
  dplyr::mutate(sum = sum(Med_Intensity),
                relative = Med_Intensity / sum) %>%
  arrange(relative) %>%
  mutate(ranks = order(order(relative, decreasing=TRUE))) %>% 
  ggplot(aes(x = ranks, y = Med_Intensity, 
             colour=Group, group=Group)) +
  geom_point(aes(colour=Group, group=Group), size=0.6
  ) +
  scale_colour_manual(values=c("hotpink", "darkblue", "grey45")) +
  scale_y_continuous(breaks=seq(15, 30, 5), limits=c(12,32)) +
  scale_x_continuous(breaks=c(0,414,425), limits=c(0,430)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x= element_text(size=8),
        axis.text.y= element_text(size=8),
        axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),
        plot.title = element_text(size=10)) +
  labs(
    x = "Relative rank",
    y = expression('log'[2]*' LFQ intensity'),
    col = "",
    title = ""
  )

#Figure 2C Number of shared and unique proteins between groups, AML and controls
list_input_genes <- df_prot_valid %>%
  dplyr::select(pheno_without_qc$Samples, UID) %>%
  gather("Samples", "Intensity", -UID) %>%
  filter(Intensity != 0) %>%
  filter(is.finite(Intensity)) %>%
  left_join(df_uid_information, "UID") %>%
  left_join(pheno_without_qc, "Samples") %>%
 dplyr::select(Genes, Intensity, Protein.Group, UID, Group) %>%
  group_by(Group, UID, Genes, Protein.Group) %>%
  summarise(Intensity = mean(Intensity))  %>%
  left_join(vvs %>%
              dplyr::select(-No_pass) %>%
              gather("DIAGNOSIS", "VVSValue", -UID), NULL) %>%
  filter(VVSValue >= 40) %>% # filter out not valid with a value
    dplyr::select(-Intensity) %>% 
  pivot_wider(names_from = DIAGNOSIS, values_from = UID, values_fill = NA) %>% 
  ungroup() %>% 
  dplyr::select(CTRL, AML)  %>%
  as.list()

# Remove NA
list_input_genes <- lapply(list_input_genes, function(x) x[!is.na(x)])
fit <- euler(list("CTRL" = unique(list_input_genes$CTRL), "AML" = unique(list_input_genes$AML)))
# Plot venn diagram
venn_diagram <- plot(fit, 
                     counts = TRUE,  
                     fills = c("hotpink", "darkblue"), 
                     alpha = 0.75,
                     quantities = list(type = c('counts')),
                     edges = FALSE)
plot(venn_diagram)

### Figure 2 D-F Principal Component Analysis results

##### PCA calculation
pci <- df_prot_valid %>% 
  remove_rownames %>%
  column_to_rownames("UID") %>% 
  dplyr::select(pheno_without_qc$Samples) %>% 
  na.omit() %>%
  t() %>%
  prcomp(rank = 5)

## Figure 2D plot
PCA <- pci$x %>%
  as.data.frame() %>%
  rownames_to_column(var = "Samples") %>%
  left_join(pheno_without_qc, "Samples") %>%
  mutate(., New_lab_fac = if_else(Group == "AML", as.factor(ID), as.factor(Group))) %>% # control as one group
  ggplot(aes(x = PC1, y = PC2, fill = New_lab_fac)) +
  geom_point(aes(fill= New_lab_fac), shape=21, size = 2) + 
  theme_classic() +
  theme(plot.title = element_blank(), 
        plot.subtitle = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=7),
        legend.title = element_text(size=7),
        axis.text.x= element_text(size=8),
        axis.text.y= element_text(size=8),
        axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),
        ) +
 scale_y_continuous(limits = c(-15,15))+
 scale_x_continuous(limits = c(-15,15))+
 scale_fill_brewer( palette="Paired") +
 guides(fill = guide_legend(override.aes = list(shape = 21),
                             keywidth=0.1,
                             keyheight=0.1,
                             default.unit="inch")) +
 labs(
    fill = "Samples",
    shape = "",
    x=paste0("PC1: ",round(as.data.frame(as.matrix(summary(pci)$importance))[2,1]*100,2),"%"),
    y=paste0("PC2: ",round(as.data.frame(as.matrix(summary(pci)$importance))[2,2]*100,2),"%")
    ) 

plot(PCA)

# Figure 2E Extract the top 10 proteins contributing to PC1
PCA_loads_bars <-  pci$rotation %>%
  as.data.frame() %>%
  rownames_to_column(var = "UID") %>%
  left_join(df_uid_information, by = "UID") %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n=10) %>% # top 10 proteins with highest absolute loading plot value along PC1
  mutate(., UID = reorder(UID, abs(PC1)))%>%
  ggplot(., aes(x=abs(PC1), y=UID)) + 
  geom_bar(stat="identity", alpha=0.7, width = 0.6, fill="grey")+
  geom_text(aes(label=Genes, hjust = "right"),size=2) +
  labs(
    y='Top 10 proteins', 
    x="PC1"
  ) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.y.right =  element_blank(),
        axis.title.x = element_text( size= 9),
        axis.text.x = element_text(size=8), 
        axis.text.y = element_blank(),
        axis.title.y =  element_blank(), 
        axis.ticks.y =element_blank(),
        legend.text = element_text(size=8),
        legend.position = "none", 
        legend.title = element_text(size=8)
  )
plot(PCA_loads_bars)


# Figure 2F  SAA1 dynamics 
# Calculate mean of control group for example protein SAA1
### 18.75779 mean calculated for control group

# Plot dynamics of example protein
SAA1_dynamics <- 
  df_prot_valid %>%
  dplyr::select(UID, pheno_without_qc$Samples) %>%
  tidyr::gather("Samples", "Intensity", -UID) %>%
  filter(., !is.na(Intensity)) %>%
  merge(.,pheno_without_qc , by ="Samples") %>%
  ungroup() %>%
  left_join(df_uid_information, "UID") %>% 
  filter(., UID == "UID242") %>% # SAA1;SAA2
  group_by(Group, ID) %>%
  arrange(ID) %>%
  filter(Group == "AML") %>% 
  mutate(., New_lab_fac = if_else(Group == "AML", as.factor(ID), as.factor(Group))) %>% 
  ggplot() +
  geom_line(aes(x= Date_Trial, y = Intensity, colour = New_lab_fac), size=1) +
  scale_colour_brewer (palette="Paired") + 
  geom_hline(yintercept=18.75779,linetype='dashed', colour= "gold") + 
  theme_classic() +
  theme(legend.position = "none",
         axis.text.x = element_text(size = 8),
         axis.text.y= element_text(size=8),
         axis.title.x = element_text(size=9),
         axis.title = element_text(size=9)) +
  labs(
    y = expression('log'[2]*' LFQ intensity'), 
    x = "Trial day",
    color = ""
  ) +
  guides(fill = guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch"))
plot(SAA1_dynamics)

######################################################################
#### STATISTICAL COMPARISON OF QUANTIFIED PROTEINS BETWEEN GROUPS #### 
######################################################################
# Set.seed
set.seed(18543)

# Prepare a pheno file and compare  AML and control samples to define statistically significant different proteins between groups 
pheno_statistics <- pheno_without_qc %>%  
  filter(., !is.na(Samples)) %>%
  mutate(DIAGNOSIS = factor(Group, levels = unique(Group))) %>%
  mutate(ct = if_else(Group =="CTRL", "CTRL", as.character(ID))) %>%
  mutate(ID_AML = if_else(Group =="AML", ID, NA))

# make a design
design = model.matrix(~0+as.factor(as.numeric(pheno_statistics$DIAGNOSIS)))
colnames(design) = levels(pheno_statistics$DIAGNOSIS)

# Duplicate correlation
df_corfit <- df_prot_imp %>%
  dplyr::select(pheno_statistics$Samples, UID) %>%
  column_to_rownames("UID")
corfit <- duplicateCorrelation(df_corfit, block = pheno_statistics$ct, design = design)

# fit in limma
fit = df_prot_imp %>%
  dplyr::select(pheno_statistics$Samples, UID) %>%
  column_to_rownames("UID") %>%
  lmFit(.,design, block=pheno_statistics$ct, cor=corfit$consensus)

# define all possible comparisons with a contrast matrix
contrast_list= combn(colnames(design), 2)
contrasts = NULL
for(i in seq(1, length(contrast_list), 2)){
  contrasts = c(contrasts, sprintf("%s - %s", contrast_list[i], contrast_list[i + 1]))
}

# Make contrast matrix
contrast_matrix = makeContrasts(contrasts = contrasts,
                                levels = design)

# fit contrasts
fit2 = contrasts.fit(fit, contrast_matrix)
ebayes = eBayes(fit2)

# get the results
m_results = decideTests(ebayes, p.value = 0.05, adjust.method = "BH", lfc = 1)

# make a dataframe with significant results
df_sigresults <- m_results %>%
  as.data.frame() %>% 
  filter(if_any(where(is.numeric), ~ (.x >= 1) | (.x <= -1))) %>%
  rownames_to_column("UID") %>%
  left_join(df_uid_information, "UID")

# evaluate tissue enrichment of significant proteins
inputGenes_DEGs_all <- df_sigresults$Genes %>% na.omit() ## proteings found significantly different between AML and Control group
inputGenes_background <- df_uid_information$Genes %>% na.omit() ## all proteins

gs_DEGs_all <-GeneSet(geneIds=inputGenes_DEGs_all,organism="Homo Sapiens",geneIdType=SymbolIdentifier()) ## Define input gene set of interest with Gene Symbols
gs_background <-GeneSet(geneIds=unique(inputGenes_background),organism="Homo Sapiens",geneIdType=SymbolIdentifier()) ## Define background input gene set with Gene Symbols

output_DEGs_all <-teEnrichment(inputGenes = gs_DEGs_all, backgroundGenes = gs_background) # calculate enrichment

seEnrichmentOutput_DEGs_all<-output_DEGs_all[[1]] # Extract tissue enrichment summary results
enrichmentOutput_DEGs_all <-setNames(data.frame(assay(seEnrichmentOutput_DEGs_all),row.names = rowData(seEnrichmentOutput_DEGs_all)[,1]), colData(seEnrichmentOutput_DEGs_all)[,1])
enrichmentOutput_DEGs_all$Tissue<-row.names(enrichmentOutput_DEGs_all)

### Figure 3A Overview of differentially quantified proteins - results from limma analysis
df_DEGs_corplot <- df_prot_imp %>%
  dplyr::select(UID, pheno_without_qc$Samples) %>% 
  left_join(df_uid_information, "UID") %>%
  filter(., UID %in% df_sigresults$UID) %>%
  mutate(ID = paste(UID, Genes, sep = "_")) %>%
  distinct(ID, .keep_all = T) %>%
  filter(!is.na(ID)) %>%
  `rownames<-`(NULL) %>%
  column_to_rownames("ID") 

DEGs_heatmap= ComplexHeatmap::Heatmap(t(scale(t(df_DEGs_corplot[2:46]))),
                                      cluster_columns=F, column_order = pheno_statistics$Samples[order(pheno_statistics$ct)], 
                                      column_split = pheno_statistics$ct,
                                      column_names_side='top',
                                      show_row_names =T,
                                      show_column_names= F, row_labels =  df_DEGs_corplot$Genes,
                                      heatmap_legend_param = list(title= "z-score", 
                                                                  legend_direction= "vertical",
                                                                  heatmap_legend_side = "bottom", 
                                                                  annotation_legend_side = "bottom", 
                                                                  legend_width= unit(2, "cm")))

ComplexHeatmap::draw(DEGs_heatmap, heatmap_legend_side='left')

### Figure 3b Dynamics of top 2 proteins with highest effect size and 2 top with lowest pvalue from limma results
myProt_colors <- brewer.pal(n=4,"PuOr")
names(myProt_colors) <- c("CAMP", "C9", "HAMP", "SAA2")

Plot_AML_sigdf_prot <-
  df_prot_valid %>%
  dplyr::select(UID, pheno_without_qc$Samples) %>%
  tidyr::gather("Samples", "Intensity", -UID) %>%
  filter(., !is.na(Intensity)) %>%
  merge(.,pheno_without_qc , by ="Samples") %>%
  ungroup() %>%
  left_join(df_uid_information, "UID") %>% 
  filter(., Genes %in% c("CAMP", "HAMP", "SAA2", "C9")) %>% 
  filter(., Group =="AML") %>% 
  arrange(ID, Date_Trial) %>% 
  ggplot(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes)) +
  geom_point(size=2) +
  geom_line(aes(group=Genes),size=0.3)+
  scale_color_manual(values=myProt_colors ) +
  scale_y_continuous(limits=c(15, 27.5, 2.5)) +
  facet_grid(~ID, scales ="free", space="free_x" ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.text.y= element_blank(), 
    axis.ticks.y = element_blank(),
    axis.title = element_text(size=9),
    strip.text.x = element_text(size=9, angle = 0, 
                                margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))) +
  labs(
    y = expression('log'[2]*' LFQ intensity'), 
    x = "Trial day",
    color = ""
  )

Plot_CTRL_sigdf_prot <-
  df_prot_valid %>%
  dplyr::select(UID, pheno_without_qc$Samples) %>%
  tidyr::gather("Samples", "Intensity", -UID) %>%
  filter(., !is.na(Intensity)) %>%
  merge(.,pheno_without_qc , by ="Samples") %>%
  ungroup() %>%
  left_join(df_uid_information, "UID") %>% 
  filter(., Genes %in% c("CAMP", "HAMP", "SAA2", "C9")) %>% 
  filter(., Group =="CTRL") %>% 
  ggplot(aes(x= ID, y = Intensity, col = Genes)) +
  geom_point(size=2) +
  scale_y_continuous(limits=c(15, 27.5, 2.5)) +
  scale_color_manual(values=myProt_colors ) +
  facet_grid(~Group, scales ="free", space="free_x" ) +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.x  = element_blank(),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=9),
        axis.title = element_text(size=9),
        strip.text.x = element_text(size=9, angle = 0, 
                                    margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))) +
  labs(
    y = expression('log'[2]*' LFQ intensity'), 
    x = "Sample",
    color = ""
  )

# combine figures
Fig_top_DEGs <- ggarrange (Plot_CTRL_sigdf_prot, Plot_AML_sigdf_prot, nrow=1, ncol=2,
                           heights = c(1, 1),
                           widths=c(0.2, 1), align = "h")


plot(Fig_top_DEGs)

######################################################################
######### CORRELATION OF EXPRESSION PATTERNS IN PROTEOMIC DATA ########## 
######################################################################

# Calculate Pearson correlation coeffcients between all proteins
df_prot_correlation_all <- df_prot_imp %>%
  left_join(df_uid_information, "UID") %>%
  mutate(ID = paste(UID, Genes, sep = "_")) %>%
  distinct(ID, .keep_all = T) %>%
  filter(!is.na(ID)) %>%
  `rownames<-`(NULL) %>%
  column_to_rownames("ID") %>%
  dplyr::select(pheno_without_qc$Samples) %>%
  t() %>% 
  Hmisc::rcorr(type = "pearson")

# NA correlation to 0
df_prot_correlation_all$r[is.na(df_prot_correlation_all$r)] <- 0

# Filter for Correlations > 0.8 to define clusters of correlation in Cytoscape
cytoscape_corr <- df_prot_correlation_all$r %>%
  as.data.frame() %>%
  rownames_to_column("Node") %>% 
  gather("Target", "Correlation",-Node) %>%
  mutate(Correlation = Correlation) %>% 
  mutate(Target = gsub(Target, pattern = ".*_", replacement = "TEST_")) %>%
  mutate(Target = gsub(Target, pattern = "TEST_", replacement = "")) %>%
  mutate(Node = gsub(Node, pattern = ".*_", replacement = "TEST_")) %>%
  mutate(Node = gsub(Node, pattern = "TEST_", replacement = "")) %>%
  filter(Correlation >= 0.80) %>%
  filter(Node != "" & Target != "") %>%
  filter(Node != Target)

cytoscape_corr = cytoscape_corr[!duplicated(t(apply(cytoscape_corr[c("Node", "Target", "Correlation")], 1, sort))), ] 

### With cytoscape_corr define correlated clusters in Cytoscape and import cluster information, list of genes of each cluster are detailed in Supplemental table 1H
clusterPlt <- c("PPBP", "THBS1", "PF4;PF4V1") # Platelet associated
clusterHb <- c( "HBB", "HBA1", "HBD", "HBB;HBD;HBG2") # Hemoglobin associated 

cluster1 <- c("IGKC","IGKV1-33;IGKV1D-33","C1R","IGLC3","IGKV3-20","IGKV2-28;IGKV2-40;IGKV2D-28;IGKV2D-40",
  "C9","C6","HP","LYVE1","IGHV3-15;IGHV3-72","ORM1","F10","IGHG1;IGHG2","KNG1","IGHG1;IGHG3","C8A",
  "A1BG","CFI","MASP2","CD14","MRC1","SAA1;SAA2","CFH;CFHR1","CRP","F11","SAA2","C1S","LRG1","RNASE4",
  "CD44","APOF","IGFBP6","SAA1","IGFBP2","CFB","AZGP1","CTSD","HAMP","ORM1;ORM2","CFHR3;CFHR4","LBP",
  "EFEMP1","F5","CFH","RNASE1","VASN","FGL1","CFHR5","CST3","REG1B","VCAM1","FCGR3A","CFHR3","PROS1","JCHAIN") # Inflammation assocaited

cluster2 <- c("GPLD1","IGHV3-64D","IGHV5-51","CPN1","IGHG1;IGHG3;IGHG4","IGHV3-7","ALB","ITIH1","TF","C1QA",
"BCHE","SERPIND1","AGT","IGHV4-34;IGHV4-38-2","ITIH2","APOL1","PGLYRP2","C5","SERPINA6","SERPINA3",
"A2M","GPX3","SERPINA7","A2M;PZP","IGHG1;IGHG2;IGHG3","IGHV4-28","IGFALS","C3","APOA1","CPN2",
"IGLC3;IGLC6;IGLC7","IGHG2;IGHG3","AFM","TTR","IGKV4-1","GSN","LCAT","IGHV3-38","PON1","SERPING1", # # endopeptidase inhibition associated
"IGHG1","SERPINA4","IGHV3-38-3","CNDP1")

# Figure 4A Pearson correlation heatmap with highlighted clusters of high correlation
annotated_correlation <- df_prot_correlation_all$r %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  separate(ID, into = c("UID", "Genes"), sep ="_") %>%
  mutate(Cluster = ifelse(Genes %in% clusterPlt, "1",  
                          ifelse(Genes %in% clusterHb, "2",
                                 ifelse(Genes %in% cluster1, "3",
                                        ifelse(Genes %in% cluster2, "4",NA)))))

annotated_correlation[is.na(annotated_correlation)] <- 0

ht_clusters= ComplexHeatmap::Heatmap(as.matrix(annotated_correlation[3:452]), 
                                     cluster_columns=T,
                                     column_names_side='top', 
                                     column_split = annotated_correlation$Cluster, 
                                     row_split = annotated_correlation$Cluster,
                                     show_row_names = F,
                                     show_column_names= F, 
                                     show_heatmap_legend = F,
                                     column_gap=unit(.000001, "npc"),
                                     row_gap=unit(.000001, "npc"),
                                     raster_quality =1,use_raster = T,
                                     heatmap_legend_param = list(title= "Pearson\ncorrelation", 
                                                                 legend_direction= "vertical",
                                                                 heatmap_legend_side = "bottom",
                                                                 annotation_legend_side = "bottom",
                                                                 legend_width= unit(2, "cm")))

ComplexHeatmap::draw(ht_clusters, heatmap_legend_side='right') 


# Figure 4B - Expand correlations from Hemoglobin associated cluster
# plot correlations of Hemoglobin associated proteins 
corr_matrix_Hb <- df_prot_correlation_all$r %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  separate(ID, into = c("UID", "Genes"), sep ="_") %>%
  mutate(Cluster = ifelse(Genes %in% clusterHb, "keep", "remove")) %>%
  filter(Cluster == "keep") %>%
  mutate(ID = paste(UID, Genes, sep = "_"))  

# remove UID from column names
names(corr_matrix_Hb) = gsub(pattern = ".*._", replacement="", x=names(corr_matrix_Hb))

# select only  data  from Hemoglobin associated cluster
corr_matrix_Hb_2 <- corr_matrix_Hb %>%
  select(c(Genes,  HBB, HBA1, HBD, `HBB;HBD;HBG2`)) %>%
  `rownames<-`(NULL) %>%
  column_to_rownames("Genes")%>% 
  t()

# set colours for coorplots
col<- colorRampPalette(c("blue", "white", "white", "red"))(50) 

# plot corrs for Hb
corrplot(corr_matrix_Hb_2, method="circle",  
         tl.cex = 0.8,
         #tl.pos='n', 
         col.lim = c(0, 1),
         cl.pos="n" ,
         type= "lower",
         order = 'original', col= col) # bg= "white")



# Figure 4B - Expand correlations from Platelet associated cluster
# plot dynamics of Platelet associated proteins
corr_matrix_Plt <- df_prot_correlation_all$r %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  separate(ID, into = c("UID", "Genes"), sep ="_") %>%
  mutate(Cluster = ifelse(Genes %in% clusterPlt, "keep", "remove")) %>%
  filter(Cluster == "keep") %>%
  mutate(ID = paste(UID, Genes, sep = "_"))  

# remove UID from column names
names(corr_matrix_Plt) = gsub(pattern = ".*._", replacement="", x=names(corr_matrix_Plt))

# select only data  from cluster 4
corr_matrix_Plt_2 <- corr_matrix_Plt %>%
  select(c(Genes, PPBP, THBS1, `PF4;PF4V1`)) %>%
  `rownames<-`(NULL) %>%
  column_to_rownames("Genes")%>%
  t()

# plot correlation data
corrplot(corr_matrix_Plt_2, method="circle",  
         tl.cex = 0.8,
         #  tl.pos='n', 
         cl.pos="n" ,
         col.lim = c(0, 1),
         type="lower", 
         order = 'original', col= col) 


# Figure 4D longitudinal dynamics of hemoglobin associated proteins in AML patients
# define colours for plot
palette_Hb<- brewer.pal(n=4,"Set1")
names(palette_Hb) <- c("HBB", "HBA1", "HBD", "HBB;HBD;HBG2")

# plot
Hb_cluster_plot <-
  df_prot_valid %>%
  dplyr::select(UID, pheno_without_qc$Samples) %>%
  tidyr::gather("Samples", "Intensity", -UID) %>%
  merge(.,pheno_without_qc , by ="Samples") %>%
  ungroup() %>%
  left_join(df_uid_information, "UID") %>% 
  filter(., Genes %in% clusterHb) %>% 
  filter(., Group_fac =="AML") %>%
  arrange(ID, Samples, Date_Trial) %>% 
  ggplot() +
  geom_point(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes), size=2) +
  geom_line(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes, group=Genes),size=0.3)+
  geom_hline(data = df_prot_valid %>%
               dplyr::select(UID, pheno_without_qc$Samples) %>%
               tidyr::gather("Samples", "Intensity", -UID) %>%
               filter(., !is.na(Intensity)) %>%
               merge(.,pheno_without_qc , by ="Samples") %>%
               ungroup() %>%
               left_join(df_uid_information, "UID") %>% 
               filter(., Genes %in% clusterHb) %>% 
               filter(., Group_fac !="AML") %>%
               group_by(Genes) %>% 
               summarise(mean = mean(Intensity, na.rm =T)),
             linetype = "dotted",
             aes(yintercept = mean, col = Genes)) +
  scale_color_manual(values=palette_Hb ) +
  scale_y_continuous(limits=c(18, 28), breaks=c(18,20,22,24,26,28)) +
  facet_grid(~ID, scales ="free", space="free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size=9),
    strip.text.x = element_text(size=9, angle = 0, 
                                margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))) +
  labs(
    y = expression('log'[2]*' LFQ intensity'), 
    x = "Trial day",
    color = ""
  )

# Figure 4E  longitudinal dynamics of platelet associated proteins in AML patients
# define colours for plot
palette_Plt<- c("#F1A340", "purple4", "#998EC3") 
names(palette_Plt) <- c("PPBP", "THBS1", "PF4;PF4V1")

#plot 
Plt_cluster_plot <-df_prot_valid %>%
  dplyr::select(UID, pheno_without_qc$Samples) %>%
  tidyr::gather("Samples", "Intensity", -UID) %>%
  merge(.,pheno_without_qc , by ="Samples") %>%
  ungroup() %>%
  left_join(df_uid_information, "UID") %>% 
  filter(., Genes %in% clusterPlt) %>% 
  filter(., Group_fac =="AML") %>%
  arrange(ID, Samples, Date_Trial) %>% 
  ggplot() +
  geom_point(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes), size=2) +
  geom_line(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes, group=Genes),size=0.3)+
  geom_hline(data = df_prot_valid %>%
               dplyr::select(UID, pheno_without_qc$Samples) %>%
               tidyr::gather("Samples", "Intensity", -UID) %>%
               filter(., !is.na(Intensity)) %>%
               merge(.,pheno_without_qc , by ="Samples") %>%
               ungroup() %>%
               left_join(df_uid_information, "UID") %>% 
               filter(., Genes %in% clusterPlt) %>% 
               filter(., Group_fac !="AML") %>%
               group_by(Genes) %>% 
               summarise(mean = mean(Intensity, na.rm =T)),
             linetype = "dotted",
             aes(yintercept = mean, col = Genes)) +
  scale_color_manual(values=palette_Plt ) +
  scale_y_continuous(limits=c(14, 20), breaks=c(14,16,18,20)) +
  facet_grid(~ID, scales ="free", space="free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size=9),
    strip.text.x = element_text(size=9, angle = 0, 
                                margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))) +
  labs(
    y = expression('log'[2]*' LFQ intensity'), 
    x = "Trial day",
    color = ""
  )

### plot all Plt and Hb dynamics together
ggpubr::ggarrange(Hb_cluster_plot 
                  , Plt_cluster_plot 
                  , ncol=1)

######################################################################
######### LONGITUDINAL TRENDS OF PROTEINS IN MAJOR CLUSTERS ########## 
######################################################################

# plot longitudinal plots for inflammation cluster ####
# define colours for plot
palette_Infl_short<- c( "#67a9cf","lightgrey", "#1c9099")
names(palette_Infl_short) <- c("VCAM1",  "PROS1", "C6")

# plot dynamics
Inflam_examples <- df_prot_valid %>%
  dplyr::select(UID, pheno_without_qc$Samples) %>%
  tidyr::gather("Samples", "Intensity", -UID) %>%
  merge(.,pheno_without_qc , by ="Samples") %>%
  ungroup() %>%
  left_join(df_uid_information, "UID") %>% 
  filter(., Genes %in% c("VCAM1",  "PROS1", "C6")) %>% 
  filter(., Group_fac =="AML") %>%
  arrange(ID, Date_Trial) %>% 
  ggplot() +
  geom_point(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes), size=2) +
  geom_line(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes, group=Genes),size=0.3)+
  geom_hline(data = df_prot_valid %>%
               dplyr::select(UID, pheno_without_qc$Samples) %>%
               tidyr::gather("Samples", "Intensity", -UID) %>%
               filter(., !is.na(Intensity)) %>%
               merge(.,pheno_without_qc , by ="Samples") %>%
               ungroup() %>%
               left_join(df_uid_information, "UID") %>% 
               filter(., Genes %in% c("VCAM1",  "PROS1", "C6")) %>% 
               filter(., Group_fac !="AML") %>%
               group_by(Genes) %>% 
               summarise(mean = mean(Intensity, na.rm =T)),
             linetype = "dotted",
             aes(yintercept = mean, col = Genes)) +
  scale_color_manual(values=palette_Infl_short ) +
  facet_grid(~ID, scales ="free", space="free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size=9),
    strip.text.x = element_text(size=9, angle = 0, 
                                margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))) +
  labs(
    y = expression('log'[2]*' LFQ intensity'), 
    x = "Trial day",
    color = ""
  )


# plot longitudinal plots for endopeptidase cluster
# define colours for example proteins
palette_Endop_short<- c("lightsalmon", "bisque3",  "grey45")
names(palette_Endop_short) <- c("SERPIND1", "C3", "APOA1" )

Endopep_examples <- df_prot_valid %>%
  dplyr::select(UID, pheno_without_qc$Samples) %>%
  tidyr::gather("Samples", "Intensity", -UID) %>%
  merge(.,pheno_without_qc , by ="Samples") %>%
  ungroup() %>%
  left_join(df_uid_information, "UID") %>% 
  filter(., Genes %in% c("SERPIND1", "C3", "APOA1")) %>% 
  filter(., Group_fac =="AML") %>%
  arrange(ID, Date_Trial) %>% 
  ggplot() +
  geom_point(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes), size=2) +
  geom_line(aes(x= as.factor(Date_Trial), y =Intensity, col = Genes, group=Genes),size=0.3)+
  geom_hline(data = df_prot_valid %>%
               dplyr::select(UID, pheno_without_qc$Samples) %>%
               tidyr::gather("Samples", "Intensity", -UID) %>%
               filter(., !is.na(Intensity)) %>%
               merge(.,pheno_without_qc , by ="Samples") %>%
               ungroup() %>%
               left_join(df_uid_information, "UID") %>% 
               filter(., Genes %in% c("SERPIND1", "C3", "APOA1")) %>% 
               filter(., Group_fac !="AML") %>%
               group_by(Genes) %>% 
               summarise(mean = mean(Intensity, na.rm =T)),
             linetype = "dotted",
             aes(yintercept = mean, col = Genes)) +
  scale_color_manual(values=palette_Endop_short ) +
  facet_grid(~ID, scales ="free", space="free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size=9),
    strip.text.x = element_text(size=9, angle = 0, 
                                margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"))) +
  labs(
    y = expression('log'[2]*' LFQ intensity'), 
    x = "Trial day",
    color = ""
  )


# plot both panels together
ggarrange (Inflam_examples , Endopep_examples, 
            nrow=2, ncol=1) 


# Tissue enrichment of proteins in major clusters with proteins highly correlated
Clust_1_2_names <- c( cluster1, cluster2) ## genes in clusters 1 and 2 - endopeptidase and inflammation

inputGenes_Clust_1_2 <- Clust_1_2_names ##  specific proteins in highest correlated clusters
inputGenes_background <- df_uid_information$Genes %>% na.omit() ## all proteins

gs_Clust_1_2 <-GeneSet(geneIds=inputGenes_Clust_1_2,organism="Homo Sapiens",geneIdType=SymbolIdentifier()) ## Define input gene set of interest with Gene Symbols
gs_background <-GeneSet(geneIds=unique(inputGenes_background),organism="Homo Sapiens",geneIdType=SymbolIdentifier()) ## Define background input gene set with Gene Symbols

output_Clust_1_2 <-teEnrichment(inputGenes = gs_Clust_1_2, backgroundGenes = gs_background) # calculate enrichment

seEnrichmentOutput_Clust_1_2 <-output_Clust_1_2[[1]] # Extract tissue enrichment summary results
enrichmentOutput_Clust_1_2 <-setNames(data.frame(assay(seEnrichmentOutput_Clust_1_2),row.names = rowData(seEnrichmentOutput_Clust_1_2)[,1]), colData(seEnrichmentOutput_Clust_1_2)[,1])
enrichmentOutput_Clust_1_2$Tissue<-row.names(enrichmentOutput_Clust_1_2)
