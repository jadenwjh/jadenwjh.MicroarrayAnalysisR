#PRE-PROCESSING


library(GEOquery); library(affy); library(limma); library(oligo); library(readr); library(ggplot2); library(ggpubr); library(ReactomePA); library(formattable); library(reactome.db); library(huex10sttranscriptcluster.db)

#----------------------------
#uncomment to unpack CEL files
#set directory to file location

#untar("GSE143150_RAW.tar",list=TRUE)  ## check contents
#untar("GSE143150_RAW.tar")
#list.files(pattern="*.CEL.gz")
#file.rename(list.files(pattern="*.CEL.gz"), paste0("GSM4250",986:997,".CEL.gz"))
#list.files()
#-----------------------------
#-------- Data extraction from files
#-------- Functions for annotating IDs and Graph plotting
#-----------------------------

Annot <- data.frame(SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=","),
                    GENENAME=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=","),
                    ENSEMBLID=sapply(contents(huex10sttranscriptclusterENSEMBL), paste, collapse=","),
                    ENTREZID=sapply(contents(huex10sttranscriptclusterENTREZID), paste,collapse=","))

gse <- getGEO('GSE143150',GSEMatrix = F)
genotype <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2] #culture treatment is 3
}
treatment <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][3] #EMT promoted
}

annotate_id <- function(x) {
  y <-rownames(x)
  z <-as.character(as.factor(x$logFC))
  item <- matrix(data=NA, nrow=length(y), ncol=5)
  i <- 1
  for (ids in y) {
    item[i,1] <- ids
    temp <- Annot[grep(ids,rownames(Annot)),]
    item[i,2] <- as.character(temp[[1]])
    item[i,3] <- as.character(temp[[2]])
    item[i,4] <- as.character(temp[[3]])
    item[i,4] <- as.character(temp[[3]])
    p <- as.character(temp[[4]])
    item[i,5] <- p
    i <- i+1
  }
  colnames(item) <- c('transcript_cluster_id', 'SYMBOL','GENENAME','ENSEMBLID','ENTREZID')
  item <- as.data.frame(item)
  return(item)
}

findplot <- function(goi) {  #insert gene of interest to function in string
  goi <- goi
  df <-data.frame(Treatment=as.factor(sapply(GSMList(gse),treatment)),
                  Genotype=as.factor(sapply(GSMList(gse),genotype)),
                  Expression_Level = as.factor(eset[goi,])) #put GOI inside [] 
  df$Genotype<- as.factor(df$Genotype)
  levels(df$Genotype) <- c("KD","WT")
  df$Treatment<- as.factor(df$Treatment)
  levels(df$Treatment) <- c("DIFF","FULL")
  
  wtdiff <- mean(as.numeric(as.character(subset(df$Expression_Level, df$Treatment == 'DIFF' & df$Genotype == 'WT'))))
  wtfull <- mean(as.numeric(as.character(subset(df$Expression_Level, df$Treatment == 'FULL' & df$Genotype == 'WT'))))
  mdiff <- mean(as.numeric(as.character(subset(df$Expression_Level, df$Treatment == 'DIFF' & df$Genotype == 'KD'))))
  mfull <- mean(as.numeric(as.character(subset(df$Expression_Level, df$Treatment == 'FULL' & df$Genotype == 'KD'))))
  
  shell <-data.frame(Treatment=as.factor(c('DIFF','FULL','DIFF','FULL')),
                     Genotype=as.factor(c('WT','WT','KD','KD')),
                     Expression_Level = as.factor(c(wtdiff,wtfull,mdiff,mfull)))
  shell$Expression_Level <- as.numeric(as.character(shell$Expression_Level))
  shell$Expression_Level <- round(shell$Expression_Level ,digit=2)
  shell$Expression_Level<- as.factor(shell$Expression_Level)
  
  g <- ggplot(data=shell,
              aes(x=Genotype,y=Expression_Level,group=Treatment)) +
    geom_line(aes(color=Treatment)) +
    geom_point(aes(color=Treatment)) +
    ggtitle(paste(as.character(Annot[grep(goi,rownames(Annot)),][[2]])))
  g
}
#------------------------
#------------------------ Processing expression Data and applying eBayes
#------------------------

apd <- data.frame(treatment=as.factor(sapply(GSMList(gse),treatment)),genotype=as.factor(sapply(GSMList(gse),genotype)))
apd$cond <- as.factor(paste(apd$treatment,apd$genotype,sep="_"))
levels(apd$cond) <- c("DIFF_KD","DIFF_WT","FULL_KD","FULL_WT")
acelfiles <- paste0(rownames(apd),'.CEL.gz')
data <- read.celfiles(acelfiles,phenoData = new("AnnotatedDataFrame",as.data.frame(apd)))


expression_data <- oligo::rma(data) # Background correction, Normalistation using rma() on dataset
eset <- exprs(expression_data)
model <- model.matrix( ~ 0 + expression_data$cond) #linear model, with intercept and the coefficient for all conditions ("DIFF_KD","DIFF_WT","FULL_KD","FULL_WT")
colnames(model) <- levels(expression_data$cond)
contrasts <- makeContrasts(DIFF_KD - DIFF_WT, #Contrast between genotypes(KD and WT) in DIFF medium
                           FULL_KD - FULL_WT, #Contrast between genotypes(KD and WT) in FULL medium
                           FULL_KD - DIFF_KD, #Contrast between media(DIFF and FULL) in KD genotype
                           FULL_WT - DIFF_WT, #Contrast between media(DIFF and FULL) in WT genotype
                           interaction=(DIFF_KD-DIFF_WT) - (FULL_KD - FULL_WT), #Contrast between different genotype with different media, also known as interaction
                           (DIFF_KD - DIFF_WT) + ( FULL_KD - FULL_WT), #Contrast between genotypes across media
                           (FULL_KD - DIFF_KD) + ( FULL_WT - DIFF_WT), #Contrast between media across genotypes
                           levels = model)

expdata_fitted_contrasts <- lmFit(expression_data,model) #Expression data undergoes Empirical Bayes method w.r.t linear model
fitted.contrasts <- contrasts.fit(expdata_fitted_contrasts,contrasts) #Subsequently fitted with the 7 differnt contrasts above
fitted.aebayes <- eBayes(fitted.contrasts) #Dataset with Empirical Bayes method applied

true_gendiff <- topTable(fitted.aebayes,coef = 1,number=Inf,p.value = 0.05,lfc=1)
true_genfull <- topTable(fitted.aebayes,coef = 2,number=Inf,p.value = 0.05,lfc=1)
true_mediakd <- topTable(fitted.aebayes,coef = 3,number=Inf,p.value = 0.05,lfc=1)
true_mediawt <- topTable(fitted.aebayes,coef = 4,number=Inf,p.value = 0.05,lfc=1)
true_intxn <- topTable(fitted.aebayes,coef = 5,number=Inf,p.value = 0.05,lfc=1)
gen_in_allmedia <- topTable(fitted.aebayes,coef = 6,number=Inf,p.value = 0.05,lfc=1)
media_in_allgen <- topTable(fitted.aebayes,coef = 7,number=Inf,p.value = 0.05,lfc=1)


#--------------------------- #table 2
GOIS <- data.frame(Genotype_in_DIFF=nrow(true_gendiff),
                   Genotype_in_FULL=nrow(true_genfull),
                   Media_in_KD=nrow(true_mediakd),
                   Media_in_WT=nrow(true_mediawt),
                   Assuming_Interaction=nrow(true_intxn),
                   Genotypes_across_media=nrow(gen_in_allmedia),
                   Media_across_genotypes=nrow(media_in_allgen))
row.names(GOIS) <- c("Transcript clusters with p-values<0.05 and logFC>1")

#--- Mapped clusters that pass cutoff to their gene names and gene IDs using "huex10sttranscriptcluster.db", annotation file for "Affymetrix Human Exon 1.0 ST Array"
#-------------- table 4a and 4b
# identify the differentially expressed genes found in KD across media
# identify the differentially expressed genes found in WT across media
#--------------

media_in_allgen1= as.data.frame(annotate_id(media_in_allgen)) #MC1
true_mediakd1= as.data.frame(annotate_id(true_mediakd)) #MC2
true_mediawt1= as.data.frame(annotate_id(true_mediawt)) #MC3

unique_genes_kd <-subset(true_mediakd1, !(GENENAME %in% media_in_allgen1$GENENAME)) #Genes_KD
unique_genes_wt <-subset(true_mediawt1, !(GENENAME %in% media_in_allgen1$GENENAME)) #Genes_WT

#--------- Parsing Reactome
# find the pathways affected during KD

KD_after_R <- enrichPathway(unique_genes_kd$ENTREZID,organism = "human", pvalueCutoff = 1, readable = T)
gene_media_kd = summary(KD_after_R)
gene_media_kd <- subset(gene_media_kd, select=c('geneID','Description','GeneRatio','BgRatio','pvalue','Count')) #5
unique_pathways_kd<-as.data.frame(gene_media_kd$Description)
colnames(unique_pathways_kd) <- "pathways affected in KD cells"

#-----------
# plotting of line plots for genes in X
#-----------

cluster_X <- ggarrange(findplot("2555490"),findplot("2639054"),findplot("2813060"),
                       ncol = 3, nrow = 1)


#-------------------------------------------- POST-PROCESSING
#-------------------------------------------- POST-PROCESSING
#-------------------------------------------- POST-PROCESSING

# number of transcripts meeting cut off for each contrast (2)
# Vplots  (3a,3b,3c) --> looking into 3 contrasts (1:Varying media in KD, 2: Varying media in WT, 3: Varying media across genotype)
# Mapped clusters that pass cutoff to their gene names and gene IDs using "huex10sttranscriptcluster.db", annotation file for "Affymetrix Human Exon 1.0 ST Array"
# identify the differentially expressed genes found in KD across media (4a) <20>
# identify the differentially expressed genes found in WT across media (4b) <2>
# Parsed REACTOME: Pathways affected in KD across media , WT across media, WT and KD across media (5a,5b,5c)
# Line plots of classified genes X (6)

formattable(GOIS) #table 2

#194 GOIs <MEDIA_KD> #table 3a
volcanoplot(fitted.aebayes,coef = 3, main=sprintf("%d features (Between media in KD) pass cutoff [LOG FOLD CHANGE >1 , P-VALUE<0.05]",nrow(true_mediakd))); points(true_mediakd[['logFC']],-log10(true_mediakd[['P.Value']]),col='red')

#89 GOIs <MEDIA_WT> #table 3b
volcanoplot(fitted.aebayes,coef = 4, main=sprintf("%d features (Between media in WT) pass cutoff [LOG FOLD CHANGE >1 , P-VALUE<0.05]",nrow(true_mediawt))); points(true_mediawt[['logFC']],-log10(true_mediawt[['P.Value']]),col='red')

#576 GOIs  #table 3c
volcanoplot(fitted.aebayes,coef = 7, main=sprintf("%d features (Between media across genotypes) pass cutoff [LOG FOLD CHANGE >1 , P-VALUE<0.05]",nrow(media_in_allgen))); points(media_in_allgen[['logFC']],-log10(media_in_allgen[['P.Value']]),col='red')

formattable(unique_genes_kd)  #table 4a <20 genes>
formattable(unique_genes_wt) #table 4b <2 genes>

formattable(head(unique_pathways_kd,n=10)) #table 5 <162 pathways>
write.csv(unique_pathways_kd,"fig5.csv")

cluster_X #table 6
