##### Libraries ####
library(BSgenome.Hsapiens.UCSC.hg19)
library(SomaticSignatures)
{
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(mclust)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(survival)
  library(ranger)
  library(ggfortify)
  library(reshape2) 
}

##### Data Import & cleaning ####

#### extracting meta data ####
all_data=fread("Desktop/pancreatic cancer/Data/simple_somatic_mutation.open.tsv")
dss <- fread("Desktop/pancreatic cancer/Data/dss.tsv")
d <- fread("Desktop/pancreatic cancer/Data/d.tsv")


meta.spec=fread("/Users/Amin-PC/Desktop/work/data/specimen.tsv")
meta.sample=fread("/Users/Amin-PC/Desktop/work/data/sample.tsv")
meta.donor=fread("/Users/Amin-PC/Desktop/work/data/donor.tsv")

dss=unique(all_data[,c("icgc_donor_id","icgc_specimen_id","icgc_sample_id")])
dss=dss[which(nchar(as.character(dss$icgc_sample_id)) != 0 &
                nchar(as.character(dss$icgc_specimen_id)) != 0 &
                nchar(as.character(dss$icgc_donor_id)) != 0),]

a=meta.donor[which(meta.donor$donor_vital_status%in%c("deceased","alive") &
                     nchar(meta.donor$donor_survival_time)!=0 &
                     meta.donor$donor_survival_time!=0 &
                     is.na(meta.donor$donor_survival_time)!=T &
                     nchar(meta.donor$donor_age_at_diagnosis)!=0 &
                     meta.donor$donor_age_at_diagnosis!=0 &
                     is.na(meta.donor$donor_age_at_diagnosis)!=T 
),
c("icgc_donor_id","donor_vital_status","donor_survival_time","donor_age_at_diagnosis")]


b=subset(meta.spec,meta.spec$specimen_type%in%c("Primary tumour - solid tissue","Primary tumour - other"
                                                ,"Metastatic tumour - other","Metastatic tumour - NOS"))
b=b[,c("icgc_specimen_id","specimen_type","tumour_grade","tumour_stage","tumour_histological_type")]

c=merge(dss,a)

d=merge(dss,b)

meta.all=merge(c,b,"icgc_specimen_id")

# fwrite(meta.all,"../../../All meta data.tsv",sep = "\t")

####

idata=all_data[,c("icgc_sample_id","icgc_specimen_id","icgc_donor_id","project_code","chromosome"
                  ,"chromosome_start","chromosome_end","chromosome_strand","reference_genome_allele"
                  ,"mutated_to_allele","gene_affected","transcript_affected","sequencing_strategy"
                  ,"consequence_type")]
remove(all_data)
#meta.all=fread("../../../All meta data.tsv")

idata <- readRDS("Desktop/pancreatic cancer/Data/idata.rds")
d <- readRDS("Desktop/pancreatic cancer/Data/d.rds")
dss <- readRDS("Desktop/pancreatic cancer/Data/dss.rds")

idata=idata[which(as.character(idata$reference_genome_allele) %in% c('A','C','G','T') 
                  & as.character(idata$mutated_to_allele) %in% c('A','C','G','T') 
                  & nchar(as.character(idata$icgc_donor_id)) != 0 
                  & as.character(idata$chromosome) %in% c(1:22, 'X','Y','M','MT')
                  & nchar(as.character(idata$gene_affected)) != 0 
                  & nchar(as.character(idata$project_code)) != 0 
                  & nchar(as.character(idata$icgc_sample_id)) != 0 
                  & nchar(as.character(idata$icgc_specimen_id)) != 0 
                  & idata$icgc_sample_id %in% d$icgc_sample_id
), ]

idata$chromosome <- paste0('chr',idata$chromosome)
idata$chromosome[which(idata$chromosome == 'chrMT')] <- 'chrM'

# protein

ens.prot= fread(file="Desktop/pancreatic cancer/file_971130/ensemble_protein_coding_hg19.tsv")
protein = subset(idata, idata$gene_affected %in% ens.prot$gene_id)

a=as.data.table(mutationContext(
  makeVRangesFromGRanges(
    makeGRangesFromDataFrame(
      protein, keep.extra.columns = T,
      seqnames.field = 'chromosome',
      start.field = 'chromosome_end',
      end.field = 'chromosome_end'),
    ref.field = 'reference_genome_allele',
    alt.field = 'mutated_to_allele'),
  ref = Hsapiens, k = 3))

protein$motif=paste(a$alteration,sep = "-",a$context)
#fwrite(protein,"c://Users/Amin-PC/Desktop/all protein coding.csv")
protein2=protein[,c("icgc_sample_id","icgc_specimen_id","icgc_donor_id","gene_affected","motif","chromosome_end","chromosome")]
protein2=protein2[which(as.character(protein2$chromosome) %in% c("chr16","chr10","chr1","chr6","chr3","chr2","chr5","chr7","chr8","chr4","chr12",
                                                                 "chr14","chr11","chr20","chr13","chr9","chrX","chr21","chr19","chr17","chr15",
                                                                 "chr18","chr22","chrY")),]

protein2=unique(protein2)
fwrite(protein2,"Desktop/pancreatic cancer/Data/protein2.tsv",sep = "\t")

remove(a)
remove(protein)

#lncRNA

ens_lncRNA=fread("Desktop/pancreatic cancer/Data/ensemble_lncRNA_hg19.tsv")
lncRNA <- subset(idata, idata$gene_affected %in% ens_lncRNA$gene_id)

a=as.data.table(mutationContext(
  makeVRangesFromGRanges(
    makeGRangesFromDataFrame(lncRNA, keep.extra.columns = T,
                             seqnames.field = 'chromosome',
                             start.field = 'chromosome_end',
                             end.field = 'chromosome_end'),
    ref.field = 'reference_genome_allele',
    alt.field = 'mutated_to_allele'),
  ref = Hsapiens, k = 3))

lncRNA$motif=paste(a$alteration,sep = "-",a$context)
#fwrite(lncRNA,"c://Users/Amin-PC/Desktop/all lncRNA coding.tsv",sep = "\t")
lncRNA2=lncRNA[,c("icgc_sample_id","icgc_specimen_id","icgc_donor_id","gene_affected","motif","chromosome_end","chromosome")]
lncRNA2=lncRNA2[which(as.character(lncRNA2$chromosome) %in% c("chr16","chr10","chr1","chr6","chr3","chr2","chr5","chr7","chr8","chr4","chr12",
                                                              "chr14","chr11","chr20","chr13","chr9","chrX","chr21","chr19","chr17","chr15",
                                                              "chr18","chr22","chrY")),]
lncRNA2=unique(lncRNA2)
lncRNA2$gm=paste0(lncRNA2$gene_affected,"-",lncRNA2$motif)
fwrite(lncRNA2,"Desktop/pancreatic cancer/Data/lncRNA2.tsv",sep = "\t")

remove(a)
remove(lncRNA)

##### Fitting & Feature Extraction####

#table of genes
gt=as.data.table(table(unique(protein2[,c("icgc_sample_id","gene_affected")])$gene_affected))

#table of gene-motifs
protein2$gm=paste0(protein2$gene_affected,"-",protein2$motif)
gmt=as.data.table(table(unique(protein2[,c("icgc_sample_id","gm")])$gm))

#Empirical distribution gene extraction
k=ecdf(x = gt$N)
mini=1000
for (i in 1:max(gt$N)) {
  if (k(i)>=0.99) {
    if(i < mini){
      mini <- i
    }
  }
}

sig_g_01=gt[which(gt$N>=mini),]

mini=1000
for (i in 1:max(gt$N)) {
  if (k(i)>=0.999) {
    if(i < mini){
      mini <- i
    }
  }
}

sig_g_001=gt[which(gt$N>=mini),]

#Empirical distribution gene-motif extraction
j=ecdf(x = gmt$N)
mini=1000
for (i in 1:max(gmt$N)) {
  if (j(i)>=0.99) {
    if(i < mini){
      mini <- i
    }
  }
}

sig_gm_01=gmt[which(gmt$N>=mini),]

mini=1000
for (i in 1:max(gmt$N)) {
  if (j(i)>=0.999) {
    if(i < mini){
      mini <- i
    }
  }
}

sig_gm_001=gmt[which(gmt$N>=mini),]

#sig. gene-motif sig. gene (Features)

sig_gm_01$gene=substr(sig_gm_01$V1,1,nchar(as.character(sig_gm_01$V1))-7)
sig_g_sig_gm_01=subset(sig_gm_01,sig_gm_01$gene%in%sig_g_01$V1)

sig_gm_001$gene=substr(sig_gm_001$V1,1,nchar(as.character(sig_gm_001$V1))-7)
sig_g_sig_gm_001=subset(sig_gm_001,sig_gm_001$gene%in%sig_g_001$V1)

##### Clustering ####

#Preparing data matrix for clustering

gm_data=unique(subset(protein2,protein2$gm%in%sig_g_sig_gm_01$V1)[,c("icgc_sample_id","gm")])
colnames(gm_data)<-c('Sample','motif')
gene_sample <- gm_data  %>% mutate(value=1) %>% complete(motif,Sample,fill=list(value=0)) %>%
  mutate(key=paste0(motif)) %>%
  group_by(Sample,key) %>%
  summarize(value = sum(value)) %>%
  spread(key,value) %>% 
  as.data.frame
rownames(gene_sample) <- gene_sample[,1] 
gene_sample <- gene_sample[,-1]

#Finding the excluded patients during Feature extraction
sample.ids=as.data.frame(unique(protein2$icgc_sample_id))
excluded=matrix()
if(nrow(gene_sample)!=nrow(sample.ids)){
  j=0
  for(i in 1:nrow(sample.ids)){
    if(is.na(match(sample.ids[i,],rownames(gene_sample)))){
      j=j+1
      excluded[j]=as.matrix(sample.ids[i,])
    }
  }
}

a=matrix(data = 0,nrow = length(excluded),ncol = ncol(gene_sample))
rownames(a)=excluded
colnames(a)=colnames(gene_sample)

gene_sample=rbind(gene_sample,a)

remove(a)
# Clustering

cluster=function(sample.motif.data){
  y=list()
  a=Mclust(sample.motif.data)
  assignments = as.numeric(a$classification)
  for (i in 1:max(assignments)) {
    y[[i]] = sample.motif.data[which(assignments == i),]
  }
  return(y)
}

cookie=function(input, s){
  hh <<- hh + 1 #counter for iteration
  print(hh)     #counter for iteration
  for(i in 1:length(input)){
    
    if(length(input) == 1){   #start for indexing 
      ss <- s
    } else {
      if(s == ''){
        ss <- as.character(i)
      } else {
        ss <- paste0(s,'-',i)
      }
      
    }
    if(dim(input[[i]])[1] == 1){
      L[[length(L)+1]] <<- input[[i]]
      index <<- c(index, ss)
    } else {
      out=cluster(input[[i]])
      if(length(out)==1){
        L[[length(L)+1]]<<-out[[1]]
        index <<- c(index, ss)
        
      } else{
        cookie(out, ss)
      }
    }
  }
}

L=list()
index = c()
hh <- 0
cookie(list(gene_sample),'')

for (i in 1:length(L)) {
  print(nrow(L[[i]]))
}

index
# saveRDS(L,"Desktop/pancreatic cancer/Data/L.rds")


##### Data Import ####

{
  L <- readRDS("Desktop/pancreatic cancer/Data/L.rds")
  lncRNA2 <- fread("Desktop/pancreatic cancer/Data/lncRNA2.tsv", sep="\t")
  protein2 <- fread("Desktop/pancreatic cancer/Data/protein2.tsv", sep="\t")
  ens.prot= fread(file="Desktop/pancreatic cancer/Data/ensemble_gene_types_19.tsv", sep="\t")
  ens_lncRNA = fread(file="Desktop/pancreatic cancer/Data/ensemble_lncRNA_hg19.tsv", sep="\t")
  sig_g_01 <- readRDS("Desktop/pancreatic cancer/Data/sig_g_01.rds")
  sig_gm_01 <- readRDS("Desktop/pancreatic cancer/Data/sig_gm_01.rds")
  sig_g_sig_gm_01 <- readRDS("Desktop/pancreatic cancer/Data/sig_g_sig_gm_01.rds")
  boxplot_colors_fill = c("lightblue","burlywood","gray","green","deepskyblue")
  boxplot_colors_colors = c("darkblue","burlywood4","gray27","green4","deepskyblue4")
  boxplot_subtypes_colors = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
  # empirical clusters 
  EC <- lapply(L, function(o){
    return(as.data.table(subset(protein2,protein2$icgc_sample_id %in% row.names(o))))
  })
  for(i in 1:length(EC)){
    EC[[i]]$gm = paste0(EC[[i]]$gene_affected,"-",EC[[i]]$motif)
  }
}


#### significant gene frequency plot ####
{
  sgEC=lapply(EC, function(o){
    return(as.data.table(subset(o,o$gene_affected%in%as.character(sig_g_01$V1))))
  })
  u=data.table()
  grtEC=lapply(1:length(sgEC), function(i){
    u=as.data.table(table(unique(sgEC[[i]][,c("icgc_sample_id","gene_affected")])$gene_affected))
    u=cbind(u,as.data.frame(as.data.frame((as.data.frame(table(unique(sgEC[[i]][,c("icgc_sample_id","gene_affected")])$gene_affected)))$Freq/length(unique(EC[[i]]$icgc_sample_id)))))
    colnames(u)=c("gene","Freq","rel.freq")
    return(u)
  })
  
  sga=list()
  exsga=list()
  for(i in 1:length(EC)){
    sga[[i]]=as.data.frame(grtEC[[i]][,c("gene","rel.freq")])
    exsga[[i]] = data.table()
    exsga[[i]]$V1=sig_g_01[-which(sig_g_01$V1 %in% sga[[i]]$gene),"V1"]
    exsga[[i]]$freq=0
    colnames(exsga[[i]])=c("gene","rel.freq")
    sga[[i]]=rbind(sga[[i]],exsga[[i]])
    sga[[i]]=sga[[i]][order(as.character(sga[[i]]$gene)),]
  }
  
  ang_dis <- function(a,b)
  {
    if(((a%*%b)/(sqrt((a%*%a)*(b%*%b)))) > 1){return(0)}
    return((2*acos((a%*%b)/(sqrt((a%*%a)*(b%*%b)))))/pi)
  }
  ang_sim <- function(a,b)
  {
    return(1-ang_dis(a,b))
  }
  # significant gene overall level of mutation and patterns comparison
  sga_comp=matrix(ncol = 3,nrow = 20)
  kk=1
  # run ang_sim <- and and_dis <- first
  for (i in c(1,2,3,4,10)) {
    for(j in c(1,2,3,4,10)) {
      if(i!=j){
        sga_comp[kk,1]=paste0("S",i," vs. ","S",j)
        sga_comp[kk,2]=t.test(sga[[i]]$rel.freq,sga[[j]]$rel.freq)$p.value
        sga_comp[kk,3]=ang_sim(sga[[i]]$rel.freq,sga[[j]]$rel.freq)
        kk=kk+1
      }
    }
  }
  sga_comp=as.data.table(sga_comp)
}

#### set threshold
# threshold_line = vector()
# for (i in c(1:14)) {
#   threshold_line[i] = 0.5
#   if (i %in% c(3,4,10)){
#     threshold_line[i] = 0.95
#   }
# }
# for (i in c(1:length(sga))){
#   sga[[i]][which(sga[[i]]$rel.freq < threshold_line[i]),"rel.freq"] <- 0
# }





plot = list(data)
gg=1
boxplot_colors_fill = c("lightblue","burlywood","gray","green","deepskyblue")
for (o in c(1,2,3,4,10)){
  data = data.frame(Significant_Genes = c(1:length(sga[[o]]$rel.freq)), Frequency = sga[[o]]$rel.freq)
  plot[[o]] = ggplot(data=data , aes(x=Significant_Genes, y=Frequency)) +
    geom_bar(stat="identity",width = 0.5, color=boxplot_colors_colors[[gg]])+
    ylim(0, 1)+ggtitle(paste0("PCS ",gg))+ylab("Mutational frequency")+xlab("Significant genes")+
    theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
          # Remove panel grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18.0),
          axis.title = element_text(size = 22.0)
    ) 
  
  # +geom_hline(yintercept=threshold_line[o], col="red", linetype="dotted")
  gg=gg+1
}

multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[10]], ncol = 2, nrow = 3)
annotate_figure(multi, top = text_grob("Figure 2 - Mutational load in significant genes"
                                       , color = "black", size = 36))

#### all gene frequency plot ####

u=data.frame()
agtEC=lapply(1:length(EC), function(i){
  u=as.data.frame(table(unique(EC[[i]][,c("icgc_sample_id","gene_affected")])$gene_affected))
  u=cbind(u,as.data.frame(as.data.frame((as.data.frame(table(unique(EC[[i]][,c("icgc_sample_id","gene_affected")])$gene_affected)))$Freq/length(unique(EC[[i]]$icgc_sample_id)))))
  colnames(u)=c("gene","Freq","rel.freq")
  return(u)
})

ga=list()
exga=list()
for(i in 1:length(EC)){
  ga[[i]]=agtEC[[i]][,c("gene","rel.freq")]
  exga[[i]]=ens.prot[-which(ens.prot$gene_id%in%ga[[i]]$gene),"gene_id"]
  exga[[i]]$freq=0
  colnames(exga[[i]])=c("gene","rel.freq")
  ga[[i]]=rbind(ga[[i]],exga[[i]])
  ga[[i]]=ga[[i]][order(as.character(ga[[i]]$gene)),]
}

# all gene overall level of mutation and patterns comparison
ga_comp=matrix(ncol = 3,nrow = 20)
kk=1
for (i in c(1,2,3,4,10)) {
  for(j in c(1,2,3,4,10)) {
    if(i!=j){
      ga_comp[kk,1]=paste0("S",i," vs. ","S",j)
      ga_comp[kk,2]=t.test(ga[[i]]$rel.freq,ga[[j]]$rel.freq)$p.value
      ga_comp[kk,3]=ang_sim(ga[[i]]$rel.freq,ga[[j]]$rel.freq)
      kk=kk+1
    }
  }
}
ga_comp=as.data.table(ga_comp)
ga <- readRDS("Desktop/pancreatic cancer/Data/ga.rds")

plot = list()
gg=1
for (o in c(1,2,3,4,10)){
  data = data.frame(Significant_Genes = c(1:length(ga[[o]]$rel.freq)), Frequency = ga[[o]]$rel.freq)
  plot[[o]] = ggplot(data=data , aes(x=Significant_Genes, y=Frequency)) +
    geom_bar(stat="identity",width = 0.5, color=boxplot_subtypes_colors[[gg]])+
    ylim(0, 1)+ggtitle(paste0("PCS",gg))+ylab("Mutational frequency")+xlab("Protein coding genes")+
    theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
          # Remove panel grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18.0),
          axis.title = element_text(size = 22.0)
    ) 
  
  # +geom_hline(yintercept=threshold_line[o], col="red", linetype="dotted")
  gg=gg+1
  
}
multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[10]], ncol = 1, nrow = 5)
annotate_figure(multi, top = text_grob("Figure 3 - Mutational load in all protein coding genes"
                                       , color = "black", face = "bold", size = 36))




plot[[1]]
plot[[2]]
plot[[3]]
plot[[4]]
plot[[10]]


##### Frequency significant gene-motif plot ####
sgmEC=lapply(EC,function(o){
  return(subset(o,o$gm%in%sig_gm_01$V1))
})

v=data.frame()
sgmtEC=lapply(1:length(sgmEC), function(i){
  u=as.data.frame(table(unique(sgmEC[[i]][,c("icgc_sample_id","gm")])$gm))
  u=cbind(u,as.data.frame(as.data.frame((as.data.frame(table(unique(sgmEC[[i]][,c("icgc_sample_id","gm")])$gm)))$Freq/length(unique(EC[[i]]$icgc_sample_id)))))
  colnames(u)=c("gm","Freq","rel.freq")
  return(u)
})

sgma=list()
exsgma=list()
for(i in 1:length(EC)){
  sgma[[i]]=sgmtEC[[i]][,c("gm","rel.freq")]
  if (nrow(sig_gm_01[-which(sig_gm_01$V1%in%sgma[[i]]$gm),"V1"]) == 0) {
    exsgma[[i]]$V1=""
  }else {
    exsgma[[i]]=sig_gm_01[-which(sig_gm_01$V1%in%sgma[[i]]$gm),"V1"]
  }
  # exsgma[[i]]=sig_gm_01[-which(sig_gm_01$V1%in%sgma[[i]]$gm),"V1"]
  exsgma[[i]]$freq=0
  colnames(exsgma[[i]])=c("gm","rel.freq")
  sgma[[i]]=rbind(sgma[[i]],exsgma[[i]])
  sgma[[i]]=sgma[[i]][order(as.character(sgma[[i]]$gm)),]
}

mrtEC=lapply(1:length(sgmEC), function(i){
  return((as.data.frame(table(unique(sgmEC[[i]][,c("icgc_sample_id","gm")])$gm)))$Freq/length(unique(EC[[i]]$icgc_sample_id)))
})
# significant gene-motif overall level of mutation and patterns comparison
sgm_comp=matrix(ncol = 3,nrow = 20)
kk=1
for (i in c(1,2,3,4,10)) {
  for(j in c(1,2,3,4,10)) {
    if(i!=j){
      sgm_comp[kk,1]=paste0("S",i," vs. ","S",j)
      sgm_comp[kk,2]=t.test(mrtEC[[i]],mrtEC[[j]])$p.value
      kk=kk+1
    }
  }
}
sgm_comp=as.data.table(sgm_comp)

plot = list()
gg=1
for (o in c(1,2,3,4,10)){
  data = data.frame(Significant_Genes = c(1:length(sgma[[o]]$rel.freq)), Frequency = sgma[[o]]$rel.freq)
  plot[[o]] = ggplot(data=data , aes(x=Significant_Genes, y=Frequency)) +
    geom_bar(stat="identity",color="#5DADE2", width = 0.001) +
    ylim(0, 1)+ggtitle(paste0("PCS ",gg))+ylab("Mutational frequency")+xlab("Significant gene-motifs")+
    theme(plot.title = element_text(hjust = 0.5,size = 6),panel.border = element_blank(),  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 4.0),
          axis.title = element_text(size = 5.0)
    )
  gg=gg+1
}
multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[10]], ncol = 2, nrow = 3)
annotate_figure(multi, top = text_grob("Figure 4 - Mutational load in significant Gene-motifs"
                                       , color = "black", face = "bold", size = 8))

##### sig gene sig gene-motif (Features) plot  ####

sgsmEC=lapply(EC,function(o){
  return(subset(o,o$gm%in%sig_g_sig_gm_01$V1))
})

u=data.frame()
sgmtEC=lapply(1:length(sgsmEC), function(i){
  u=as.data.frame(table(unique(sgsmEC[[i]][,c("icgc_sample_id","gm")])$gm))
  u=cbind(u,as.data.frame(as.data.frame((as.data.frame(table(unique(sgsmEC[[i]][,c("icgc_sample_id","gm")])$gm)))$Freq/length(unique(EC[[i]]$icgc_sample_id)))))
  colnames(u)=c("gm","Freq","rel.freq")
  return(u)
})

sgsgma=list()
exsgsgma=list()
for(i in 1:length(EC)){
  sgsgma[[i]]=sgmtEC[[i]][,c("gm","rel.freq")]
  exsgsgma[[i]]=sig_g_sig_gm_01[-which(sig_g_sig_gm_01$V1%in%sgsgma[[i]]$gm),"V1"]
  exsgsgma[[i]]$freq=0
  colnames(exsgsgma[[i]])=c("gm","rel.freq")
  sgsgma[[i]]=rbind(sgsgma[[i]],exsgsgma[[i]])
  sgsgma[[i]]=sgsgma[[i]][order(as.character(sgsgma[[i]]$gm)),]
}

# significant gene significant gene-motif overall level of mutation and patterns comparison
sgsgm_comp=matrix(ncol = 3,nrow = 20)
kk=1
for (i in c(1,2,3,4,10)) {
  for(j in c(1,2,3,4,10)) {
    if(i!=j){
      sgsgm_comp[kk,1]=paste0("S",i," vs. ","S",j)
      sgsgm_comp[kk,2]=t.test(sgsgma[[i]],sgsgma[[j]])$p.value
      kk=kk+1
    }
  }
}
sgsgm_comp=as.data.table(sgsgm_comp)

plot = list()
gg=1
for (o in c(1,2,3,4,10)){
  data = data.frame(Significant_Genes = c(1:length(sgsgma[[o]]$rel.freq)), Frequency = sgsgma[[o]]$rel.freq)
  plot[[o]] = ggplot(data=data , aes(x=Significant_Genes, y=Frequency)) +
    geom_bar(stat="identity",width = 0.01, color=boxplot_subtypes_colors[[gg]])+
    ylim(0, 1)+ggtitle(paste0("PCS",gg))+ylab("Mutational frequency")+xlab("Significant features")+
    theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
          # Remove panel grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18.0),
          axis.title = element_text(size = 22.0)
    ) 
  gg=gg+1
}
multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[10]], ncol = 1, nrow = 5)
annotate_figure(multi, top = text_grob("Figure 4 - Mutational load in significant features"
                                       , color = "black", face = "bold", size = 36))
plot[[1]]
plot[[2]]
plot[[3]]
plot[[4]]
plot[[10]]

##### Mutational load in All lncRNAs ####

gnmtf<-lncRNA2
lncEC <- lapply(L, function(o){
  return((subset(gnmtf,gnmtf$icgc_sample_id %in% row.names(o))))
})


u=data.frame()
lnc_grtEC=lapply(1:length(lncEC), function(i){
  u=as.data.frame(table(unique(lncEC[[i]][,c("icgc_sample_id","gene_affected")])$gene_affected))
  u=cbind(u,as.data.frame(as.data.frame((as.data.frame(table(unique(lncEC[[i]][,c("icgc_sample_id","gene_affected")])$gene_affected)))$Freq/length(unique(lncEC[[i]]$icgc_sample_id)))))
  colnames(u)=c("gene","Freq","rel.freq")
  return(u)
})

lncg=list()
exlncg=list()
for(i in 1:length(lncEC)){
  lncg[[i]]=lnc_grtEC[[i]][,c("gene","rel.freq")]
  exlncg[[i]]=ens_lncRNA[-which(ens_lncRNA$gene_id%in%lncg[[i]]$gene),"gene_id"]
  exlncg[[i]]$freq=0
  colnames(exlncg[[i]])=c("gene","rel.freq")
  lncg[[i]]=rbind(lncg[[i]],exlncg[[i]])
  lncg[[i]]=lncg[[i]][order(as.character(lncg[[i]]$gene)),]
}

# lncRNA overall level of mutation and patterns comparison
lncRNA_comp=matrix(ncol = 3,nrow = 20)
kk=1
for (i in c(1,2,3,4,10)) {
  for(j in c(1,2,3,4,10)) {
    if(i!=j){
      lncRNA_comp[kk,1]=paste0("S",i," vs. ","S",j)
      lncRNA_comp[kk,2]=t.test(lncg[[i]],lncg[[j]])$p.value
      #lncRNA_comp[kk,3]=ang_sim(lnc_grtEC[[i]],lnc_grtEC[[j]])
      kk=kk+1
    }
  }
}
lncRNA_comp=as.data.table(lncRNA_comp)

plot = list()
gg=1
for (o in c(1,2,3,4,10)){
  data = data.frame(Significant_Genes = c(1:length(lncg[[o]]$rel.freq)), Frequency = lncg[[o]]$rel.freq)
  plot[[o]] = ggplot(data=data , aes(x=Significant_Genes, y=Frequency)) +
    geom_bar(stat="identity",width = 0.01, color=boxplot_colors_colors[[gg]])+
    ylim(0, 1)+ggtitle(paste0("PCS ",gg))+ylab("Mutational frequency")+xlab("Significant genes")+
    theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
          # Remove panel grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18.0),
          axis.title = element_text(size = 22.0)
    ) 
  gg=gg+1
}
multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[10]], ncol = 1, nrow = 5)
annotate_figure(multi, top = text_grob("Figure 12 - Mutational load in all lncRNAs"
                                       , color = "black", face = "bold", size = 36))

##### Mutational load in Significant gene-motif lncRNAs ####

lnc_mrtEC=lapply(1:length(lncEC), function(i){
  return((as.data.frame(table(unique(lncEC[[i]][,c("icgc_sample_id","gm")])$gm)))$Freq/length(unique(lncEC[[i]]$icgc_sample_id)))
})

v=data.frame()
lnc_mrtEC=lapply(1:length(lncEC), function(i){
  u=as.data.frame(table(unique(lncEC[[i]][,c("icgc_sample_id","gm")])$gm))
  u=cbind(u,as.data.frame(as.data.frame((as.data.frame(table(unique(lncEC[[i]][,c("icgc_sample_id","gm")])$gm)))$Freq/length(unique(EC[[i]]$icgc_sample_id)))))
  colnames(u)=c("gm","Freq","rel.freq")
  return(u)
})

m=0
xx=data.frame()
for (i in 1:length(lnc_mrtEC)) {
  a=nrow(lnc_mrtEC[[i]])
  if(a>m){
    m=a
    s=i
  }
  xx=as.data.frame(rbind(xx,as.data.frame(lnc_mrtEC[[i]]$gm)))
}
xx=unique(xx)


lncgm=list()
exlncgm=list()
for(i in c(1:(s-1),(s+1):14)){
  lncgm[[i]]=lnc_mrtEC[[i]][,c("gm","rel.freq")]
  exlncgm[[i]]=as.data.frame(lnc_mrtEC[[s]][-which(xx[,1]%in%lncgm[[i]]$gm),"gm"])
  exlncgm[[i]]$freq=0
  colnames(exlncgm[[i]])=c("gm","rel.freq")
  lncgm[[i]]=rbind(lncgm[[i]],exlncgm[[i]])
  lncgm[[i]]=lncgm[[i]][order(as.character(lncgm[[i]]$gm)),]
}
lncgm[[s]]=lnc_mrtEC[[s]][,c("gm","rel.freq")]

# lncRNA gene-motif overall level of mutation and patterns comparison
lncRNA_m_comp=matrix(ncol = 3,nrow = 20)
kk=1
for (i in c(1,2,3,4,10)) {
  for(j in c(1,2,3,4,10)) {
    if(i!=j){
      lncRNA_m_comp[kk,1]=paste0("S",i," vs. ","S",j)
      lncRNA_m_comp[kk,2]=t.test(lncgm[[i]],lncgm[[j]])$p.value
      #lncRNA_m_comp[kk,3]=ang_sim(lnc_mrtEC[[i]],lnc_mrtEC[[j]])
      kk=kk+1
    }
  }
}
lncRNA_m_comp=as.data.table(lncRNA_m_comp)

plot = list()
gg=1
for (o in c(1,2,3,4,10)){
  data = data.frame(Significant_Genes = c(1:length(lncgm[[o]]$rel.freq)), Frequency = lncgm[[o]]$rel.freq)
  plot[[o]] = ggplot(data=data , aes(x=Significant_Genes, y=Frequency)) +
    geom_bar(stat="identity",width = 0.01, color=boxplot_colors_colors[[gg]])+
    ylim(0, 1)+ggtitle(paste0("PCS ",gg))+ylab("Mutational frequency")+xlab("Significant genes")+
    theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
          # Remove panel grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18.0),
          axis.title = element_text(size = 22.0)
    ) 
  gg=gg+1
}
multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[10]], ncol = 2, nrow = 3)
annotate_figure(multi, top = text_grob("Figure 5 - Mutational load in lncRNA-motifs"
                                       , color = "black", face = "bold", size = 20))


##### common top genes ####

topgEC=list()
x=matrix(data = 0,nrow = 5,ncol = 5)
ii=1
for (o in c(1,2,3,4,10)) {
  topgEC[[ii]]=as.data.frame(table(unique(sgEC[[o]][,c("icgc_sample_id","gene_affected")])$gene_affected))
  topgEC[[ii]]=head(topgEC[[ii]][order(topgEC[[ii]]$Freq,decreasing = T),],n = 100)
  ii=ii+1
}
for (i in c(1,2,3,4,5)){
  for (o in c(1,2,3,4,5)){
    x[i,o]=nrow(subset(topgEC[[o]],topgEC[[o]]$Var1%in%topgEC[[i]]$Var1))
  }
}
# for (i in c(1,2,3,4,10)) {
#   write.csv(topgEC[[i]],paste0("C:/Users/Amin-PC/Desktop/topg_cluster",i,".csv") )
# }
colnames(x)=c("Subtype 1","Subtype 2","Subtype 3","Subtype 4","Subtype 5")
rownames(x)=c("Subtype 1","Subtype 2","Subtype 3","Subtype 4","Subtype 5")
diag(x)=0
pheatmap(x, cluster_rows = F, cluster_cols = F,display_numbers = T,
         main = "Common top Genes (Most mutated)",show_rownames = T,show_colnames = T,number_format = "%.0f")

##### region analysis ####

region=fread("C:/Users/Amin-PC/Desktop/work/DATA/side/sample.tsv",sep="\t")
region=region[,c("icgc_sample_id","project_code")]
x=data.table()
raEC=list()
for(i in 1:length(EC)){
  raEC[[i]]=unique(subset(region,region$icgc_sample_id%in%EC[[i]]$icgc_sample_id)[,c("icgc_sample_id","project_code")])
  x=rbind(x,as.data.frame(table(raEC[[i]]$project_code)))
}

##### data preparation for signature analysis ####

#region=fread("C:/Users/Amin-PC/Desktop/work/DATA/side/donor.tsv",sep="\t")
masroor.all=as.data.table(unique(idata[,c("icgc_sample_id","chromosome","chromosome_end","reference_genome_allele","mutated_to_allele")]))
colnames(masroor.all)=c("sample_id","chromosome","position","reference","mutated_to")
fwrite(masroor.all,"../../../all samples.tsv",sep = "\t")
masroor=list()
for (i in 1:length(EC)) {
  masroor[[i]]=subset(masroor.all,masroor.all$sample_id%in%EC[[i]]$icgc_sample_id)
}

for (i in 1:length(masroor)) {
  fwrite(masroor[[i]],paste0("C:/Users/Amin-PC/Desktop/cluster ",i,".tsv"),sep = "\t")
}

##### Gene association ####

for (i in 1:length(grtEC)) {
  l=subset(grtEC[[i]],grtEC[[i]]$gene%in%protein2$gene_affected)
  assign(paste0("gnmtf",i),as.data.table(l[,1:2]))
}

# for (i in 1:length(EC)){
#   
# }

gnmtf1$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf1$fisher<-0
t<-as.matrix(gnmtf1[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 68-t[i,1], 649-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf1$fisher<-t[,3]

gnmtf2$non<-gnmtf1$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf2$fisher<-0
t<-as.matrix(gnmtf2[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 259-t[i,1], 458-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf2$fisher<-t[,3]

gnmtf3$non<-gnmtf2$Freq+gnmtf1$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf3$fisher<-0
t<-as.matrix(gnmtf3[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 151-t[i,1], 566-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf3$fisher<-t[,3]

gnmtf4$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf1$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf4$fisher<-0
t<-as.matrix(gnmtf4[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 125-t[i,1], 592-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf4$fisher<-t[,3]

gnmtf5$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf1$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf5$fisher<-0
t<-as.matrix(gnmtf5[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 107-t[i,1], 610-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf5$fisher<-t[,3]

gnmtf6$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf1$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf6$fisher<-0
t<-as.matrix(gnmtf6[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 1-t[i,1], 716-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf6$fisher<-t[,3]

gnmtf7$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf1$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf7$fisher<-0
t<-as.matrix(gnmtf7[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 1-t[i,1], 716-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf7$fisher<-t[,3]

gnmtf8$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf1$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf8$fisher<-0
t<-as.matrix(gnmtf8[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 1-t[i,1], 716-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf8$fisher<-t[,3]

gnmtf9$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf1$Freq+gnmtf10$Freq+gnmtf11$Freq
gnmtf9$fisher<-0
t<-as.matrix(gnmtf7[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 2-t[i,1], 715-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf9$fisher<-t[,3]

gnmtf10$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf1$Freq+gnmtf11$Freq
gnmtf10$fisher<-0
t<-as.matrix(gnmtf7[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 1-t[i,1], 716-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf10$fisher<-t[,3]

gnmtf11$non<-gnmtf2$Freq+gnmtf3$Freq+gnmtf4$Freq+gnmtf5$Freq+gnmtf6$Freq+gnmtf7$Freq+gnmtf8$Freq+gnmtf9$Freq+gnmtf10$Freq+gnmtf1$Freq
gnmtf11$fisher<-0
t<-as.matrix(gnmtf7[,2:4])
for (i in 1:18967) {
  tt<-matrix(c(t[i,1], t[i,2], 1-t[i,1], 716-t[i,2]), nrow = 2)
  t[i,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
}
gnmtf11$fisher<-t[,3]

# Common associated genes by fisher

fisher=list()
for(i in 1:5){
  fisher[[i]]=fread(paste0("C:/Users/Amin-PC/Desktop/top_g_fisher_sc",i,".csv"))
}
y=matrix(nrow = 5,ncol = 5)
for (i in 1:5){
  for (o in 1:5){
    y[i,o]=nrow(subset(fisher[[i]],fisher[[i]][2:101,2]%in%fisher[[o]][2:101,2]))
  }
}

##### signature comparison ####


signature_all=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_all-samples/3mer_Signatures_(N=8).tsv")
signature_c1=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster1/3mer_Signatures_(N=3).tsv")
signature_c2=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster2/3mer_Signatures_(N=2).tsv")
signature_c3=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster3/3mer_Signatures_(N=4).tsv")
signature_c4=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster4/3mer_Signatures_(N=8).tsv")
signature_c10=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster10/3mer_Signatures_(N=5).tsv")

exposure_all=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_all-samples/Exposures_to_3mer_Signatures_(N=8).tsv")
exposure_c1=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster1/Exposures_to_3mer_Signatures_(N=3).tsv")
exposure_c2=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster2/Exposures_to_3mer_Signatures_(N=2).tsv")
exposure_c3=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster3/Exposures_to_3mer_Signatures_(N=4).tsv")
exposure_c4=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster4/Exposures_to_3mer_Signatures_(N=8).tsv")
exposure_c10=fread("Desktop/pancreatic cancer/Data/output_(incomplete)/results_for_cluster10/Exposures_to_3mer_Signatures_(N=5).tsv")

Alexandrov <- fread("Desktop/pancreatic cancer/Data/alex.tsv", sep="\t")
Alexandrov = Alexandrov[,c("Signature.1","Signature.2","Signature.3","Signature.5","Signature.6","Signature.13")]

ang_dis <- function(a,b)
{
  if(((a%*%b)/(sqrt((a%*%a)*(b%*%b)))) > 1){return(0)}
  return((2*acos((a%*%b)/(sqrt((a%*%a)*(b%*%b)))))/pi)
}
ang_sim <- function(a,b)
{
  return(1-ang_dis(a,b))
}

aa=signature_c1
bb=signature_c3

sigcomp = matrix(nrow = ncol(aa),ncol = ncol(bb))
rownames(sigcomp) <- paste0('set1_signature(',1:nrow(sigcomp),')')
colnames(sigcomp) <- paste0('set2_signature(',1:ncol(sigcomp),')')

for(i in 1:dim(aa)[2]){
  for(j in 1:dim(bb)[2]){
    sigcomp[i,j] <- ang_sim(aa[,i,with = F][[1]], bb[,j,with = F][[1]])
  }
}


pheatmap(sigcomp, cluster_rows = F, cluster_cols = F,display_numbers = T,
         main = "Cluster 1 VS Cluster 3")

# exposures

exposure_all_norm <- cbind(exposure_all[,1], exposure_all[,-1]/rowSums(exposure_all[,-1]))
exposure_c1_norm = cbind(exposure_c1[,1], exposure_c1[,-1]/rowSums(exposure_c1[,-1]))
exposure_c2_norm = cbind(exposure_c2[,1], exposure_c2[,-1]/rowSums(exposure_c2[,-1]))
exposure_c3_norm = cbind(exposure_c3[,1], exposure_c3[,-1]/rowSums(exposure_c3[,-1]))
exposure_c4_norm = cbind(exposure_c4[,1], exposure_c4[,-1]/rowSums(exposure_c4[,-1]))
exposure_c10_norm = cbind(exposure_c10[,1], exposure_c10[,-1]/rowSums(exposure_c10[,-1]))

all_exposures = list()
all_exposures[[1]] = exposure_c1_norm
all_exposures[[2]] = exposure_c2_norm
all_exposures[[3]] = exposure_c3_norm
all_exposures[[4]] = exposure_c4_norm
all_exposures[[10]] = exposure_c10_norm




{
  plot = list()
  gg = 1
  for (i in c(1,2,3,4,10)){
    m <- melt(all_exposures[[i]])
    m$variable <- factor(m$variable, levels = unique(m$variable))
    plot[[i]] =
      ggplot()+geom_boxplot(data = m, aes(x= variable, y= value), fill =boxplot_colors_fill[[gg]], color=boxplot_colors_colors[[gg]]
                            , width = 0.5, outlier.size = 0.1)+
      ylab("")+xlab(paste0("PCS ",gg))+
      theme(legend.position = "none",panel.border = element_blank(),  
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # Remove panel background
            panel.background = element_blank(),
            # Add axis line
            axis.line = element_line(colour = "grey"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 6.0),
            axis.title = element_text(size = 5.0))
    gg = gg+1
  }
  multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[10]], ncol = 5, nrow = 1)
  annotate_figure(multi, top = text_grob("Figure 13 - Exposure to signatures"
                                         , color = "black", face = "bold", size = 8))
}


##### difference analysis ####

par(mfrow=c(5,5))
barplot(grtEC[[1]] - t2$Freq,main = 'SC1-2 ',ylim=c(-1, 1) )
barplot(t2$Freq - t2$Freq,main = 'SC3-4 ',ylim=c(-1 ,1) )
barplot(t3$Freq - t2$Freq,main = 'SC3-5 ',ylim=c(-1, 1) )
barplot(t4$Freq - t2$Freq,main = 'SC4-5 ',ylim=c(-1, 1) )


##### Gender Analysis ####

gender=merge(dss,meta.donor[,c("icgc_donor_id","donor_sex")])

sexEC=list()
for(i in 1:14){
  sexEC[[i]]=subset(gender,gender$icgc_sample_id%in%allclusters[which(allclusters$cluster==i),]$icgc_sample_id)
  sexEC[[i]]=table(sexEC[[i]]$donor_sex)
  print(sexEC[[i]])
}





##### Transcripts analysis ####

tftable=list()
j=1
for(i in 1:14){
  tftable[[j]]=as.data.table(table(unique(idata[which(idata$icgc_sample_id%in%clusters[which(clusters[,"cluster"]==i),"icgc_sample_id"] &
                                                        nchar(idata$transcript_affected)!=0),c("icgc_sample_id","transcript_affected")])$transcript_affected))
  tftable[[j]]=tftable[[j]][order(tftable[[j]]$N,decreasing = T),]
  j=j+1
}

tfnames=as.data.table(unique(idata[which(idata$transcript_affected!=""),]$transcript_affected))
tf=list()
for(i in 1:14){
  tf[[i]]=merge(tfnames,tftable[[i]],all=T)
  tf[[i]][which(is.na(tf[[i]][,2])),2]=0
  tf[[i]]=tf[[i]][order(tf[[i]][,1]),]
  tf[[i]][,"rel.freq"]=tf[[i]][,2]/length(which(clusters[,2]==i))
}

a=unique(idata[,c("gene_affected","transcript_affected")])
colnames(a)=c("gene_affected","V1")
a=unique(a[which(a$gene_affected%in%ens.prot$gene_id),])
b=list()
c=data.frame()
for (i in 1:length(tf)){
  b[[i]]=merge(tf[[i]],a)
  b[[i]]=cbind(b[[i]],i)
  c=rbind(c,b[[i]])
}
colnames(c)=c("Transcript","N","Rel. Freq.","Gene Affected","Subtype")

plot1=ggplot(data=c[which(c$Subtype%in%c(1,2,3,4,10)&c$`Gene Affected`=="ENSG00000078328"),], aes(x=Transcript, y=N, fill=as.factor(Subtype))) +
  # geom_bar(stat="identity", position=position_dodge())
  geom_bar(stat="identity",position = "fill")+
  scale_fill_discrete(name = "",labels = c("Subtype 1", "Subtype 2", "Subtype 3","Subtype 4","Subtype 5"))+
  xlab("Transcripts of gene RBFOX1")+ylab("Relative frequency of mutations in subtypes")+ggtitle("RBFOX1 (ENSG00000078328)")+
  theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))
# geom_bar(stat="identity",position = "stack")

plot2=ggplot(data=c[which(c$Subtype%in%c(1,2,3,4,10)&c$`Gene Affected`=="ENSG00000130226"),], aes(x=Transcript, y=N, fill=as.factor(Subtype))) +
  geom_bar(stat="identity",position = "fill")+
  scale_fill_discrete(name = "",labels = c("Subtype 1", "Subtype 2", "Subtype 3","Subtype 4","Subtype 5"))+
  xlab("Transcripts of gene DPP6")+ylab("Relative frequency of mutations in subtypes")+ggtitle("DPP6 (ENSG00000130226)")+
  theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))

plot3=ggplot(data=c[which(c$Subtype%in%c(1,2,3,4,10)&c$`Gene Affected`=="ENSG00000075884"),], aes(x=Transcript, y=N, fill=as.factor(Subtype))) +
  geom_bar(stat="identity",position = "fill")+
  scale_fill_discrete(name = "",labels = c("Subtype 1", "Subtype 2", "Subtype 3","Subtype 4","Subtype 5"))+
  xlab("Transcripts of gene ARHGAP15")+ylab("Relative frequency of mutations in subtypes")+ggtitle("ARHGAP15 (ENSG00000075884)")+
  theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))

plot4=ggplot(data=c[which(c$Subtype%in%c(1,2,3,4,10)&c$`Gene Affected`=="ENSG00000153707"),], aes(x=Transcript, y=N, fill=as.factor(Subtype))) +
  geom_bar(stat="identity",position = "fill")+
  scale_fill_discrete(name = "",labels = c("Subtype 1", "Subtype 2", "Subtype 3","Subtype 4","Subtype 5"))+
  xlab("Transcripts of gene PTPRD")+ylab("Relative frequency of mutations in subtypes")+ggtitle("PTPRD (ENSG00000153707)")+
  theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))

plot5=ggplot(data=c[which(c$Subtype%in%c(1,2,3,4,10)&c$`Gene Affected`=="ENSG00000155657"),], aes(x=Transcript, y=N, fill=as.factor(Subtype))) +
  geom_bar(stat="identity",position = "fill")+
  scale_fill_discrete(name = "",labels = c("Subtype 1", "Subtype 2", "Subtype 3","Subtype 4","Subtype 5"))+
  xlab("Transcripts of gene TTN")+ylab("Relative frequency of mutations in subtypes")+ggtitle("TTN (ENSG00000155657)")+
  theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))




##### Top motif, then mutational load analysis ####

top_motifs <- lapply(EC, function(o){
  u <- as.data.frame(unique(o[,c("icgc_sample_id","gene_affected","motif")]))
  u <- as.data.frame(u[,c("icgc_sample_id","motif")])
  colnames(u)<-c('sample_id','gm')
  u <- as.data.frame(table(as.character(u$gm)))
  u <- u[order(-u$Freq),]
  u <- u[1:5,]
  return(u)
})

top_motif_genes = list()
for(i in 1:length(top_motifs)){
  top_motif_genes[[i]] <- subset(EC[[i]], EC[[i]]$motif %in% top_motifs[[i]]$Var1)
  top_motif_genes[[i]] <- unique(top_motif_genes[[i]][,c("icgc_sample_id","gene_affected")])
  top_motif_genes[[i]] <- as.data.frame(table(top_motif_genes[[i]]$gene_affected))
  top_motif_genes[[i]]$Freq <- top_motif_genes[[i]]$Freq/nrow(L[[i]])
  colnames(top_motif_genes[[i]]) <- c('gene','rel.freq')
}



plot = list()
gg=1
for (o in c(1,4)){
  data = data.frame(all_Genes = c(1:length(top_motif_genes[[o]]$gene)), Frequency = top_motif_genes[[o]]$rel.freq)
  plot[[o]] = ggplot(data=data , aes(x=all_Genes, y=Frequency)) +
    geom_bar(stat="identity",width = 0.0001, color=boxplot_colors_colors[[gg]])+
    ylim(0, 1)+ggtitle(paste0("PCS ",gg))+ylab("Mutational frequency")+xlab("Protein coding gene")+
    theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
          # Remove panel grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18.0),
          axis.title = element_text(size = 22.0)
    ) 
  gg=gg+1
}
multi = ggarrange(plot[[1]],plot[[4]], ncol = 2, nrow = 1)
annotate_figure(multi, top = text_grob("Figure 11 - Mutational load in genes with top motifs"
                                       , color = "black", face = "bold", size = 20))






##### motif rate for main subtypes ##### 
EC[[5]] <- EC[[10]]
for (i in 6:14){
  EC[[6]] <- NULL
}

selected <- list()
for (i in 1:length(EC)){
  selected[[i]] <- EC[[i]][,c('icgc_sample_id','motif')]
  selected[[i]] <- as.data.frame(table(selected[[i]]$motif))
  total_freq <- sum(selected[[i]]$Freq)
  selected[[i]]$Freq <- selected[[i]]$Freq/total_freq
  selected[[i]] <- as.data.frame(selected[[i]]$Freq)
  names(selected[[i]]) <- paste0('Signature_No_1')
  write.csv(selected[[i]], paste0('Desktop/pancreatic cancer/Motif rate input/3mer_Signatures_(N=',i,').tsv')
            ,row.names=FALSE)
}

plot <- list()
for(i in 1:length(EC)){
  plot[[i]] <- plot_signatures_3mer('Desktop/pancreatic cancer/Motif rate input/',i)
}
multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[5]],
                  ncol = 2, nrow = 3)
annotate_figure(multi, top = text_grob("Figure 20 - Motif rate in main subtypes"
                                       , color = "black", size = 36))
# print inch: 12,33,landscape


##### motif Rate for outlier subtypes##### 
all_motif = unique(protein2$motif)
selected <- list()
for (i in c(5,6,7,8,9,11,12,13,14)){
  selected[[i]] <- EC[[i]][,c('icgc_sample_id','motif')]
  selected[[i]] <- as.data.frame(table(selected[[i]]$motif))
  # adding missing motif to list for accurate plotting 
  if (length(selected[[i]]$Var1) != 96){
    missing_motif <- subset(all_motif, !(all_motif %in% selected[[i]]$Var1))
    t_df <- data.frame("Var1" = missing_motif, "Freq" = rep(0, length(missing_motif)))
    # colnames(t_df) <- c("Var1", "Freq")
    selected[[i]] <- rbind(selected[[i]], t_df)
    selected[[i]]$Var1 <- as.character(selected[[i]]$Var1)
    selected[[i]] <- selected[[i]][order(selected[[i]]$Var1),]
    remove(t_df)
  }
  total_freq <- sum(selected[[i]]$Freq)
  selected[[i]]$Freq <- selected[[i]]$Freq/total_freq
  selected[[i]] <- as.data.frame(selected[[i]]$Freq)
  names(selected[[i]]) <- paste0('Signature_No_1')
  write.csv(selected[[i]], paste0('Desktop/pancreatic cancer/Motif rate input/3mer_Signatures_(N=',i,').tsv')
            ,row.names=FALSE)
}

plot <- list()
for(i in c(5,6,7,8,9,11,12,13,14)){
  plot[[i]] <- plot_signatures_3mer('Desktop/pancreatic cancer/Motif rate input/',i)
}

multi = ggarrange(plot[[5]], plot[[6]], plot[[7]], plot[[8]], plot[[9]],
                  plot[[11]], plot[[12]], plot[[13]], plot[[14]],
                  ncol = 2, nrow = 5)
annotate_figure(multi, top = text_grob("Figure 20 - Motif rate in main subtypes"
                                       , color = "black", size = 36))


##### consequence types ####
idata <- readRDS("Desktop/pancreatic cancer/Data/idata.rds")
EC_detail <- lapply(L, function(o){
  return(as.data.table(subset(idata,idata$icgc_sample_id %in% row.names(o))))
})
for(i in 1:length(EC_detail)){
  EC_detail[[i]]$gm = paste0(EC_detail[[i]]$gene_affected,"-",EC_detail[[i]]$motif)
}
EC_detail[[5]] <- EC_detail[[10]]
for (i in 6:14){
  EC_detail[[6]] <- NULL
}

selected <- list()
cons_freq <- list()
for (i in 1:length(EC_detail)){
  selected[[i]] <- unique(EC_detail[[i]][,c('icgc_sample_id','gene_affected','consequence_type')])
  cons_freq[[i]] <- as.data.frame(table(selected[[i]]$consequence_type))
  colnames(cons_freq[[i]]) <- c("Consequence_Type", "Freq")
  cons_freq[[i]]$Freq <- cons_freq[[i]]$Freq/dim(unique(selected[[i]]))[[1]]
  cons_freq[[i]]$clust <- paste0(i)
}

cons_df <- rbind(cons_freq[[1]],cons_freq[[2]],cons_freq[[3]],cons_freq[[4]],cons_freq[[5]])
cons_df$clust <- as.factor(cons_df$clust)



plot = ggplot(data=cons_df , aes(x=Consequence_Type, y=Freq, fill = paste0('PCS',clust))) +
  geom_bar(stat="identity",position=position_dodge(), width = 0.9)+
  ylim(0, 0.6)+ylab("Relative frequency")+xlab("Consequence type")+
  theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
        # Remove panel grid lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Remove panel background
        panel.background = element_blank(),
        # Add axis line
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 90,size = 10),
        axis.text.y = element_text(size = 18.0),
        axis.title = element_text(size = 22.0)
  ) 

plot






##### Differential analysis of all protein coding genes ####
EC[[5]] <- EC[[10]]
for (i in 6:14){
  EC[[6]] <- NULL
}
selected <- list()
for (i in 1:length(EC)){
  selected[[i]] <- EC[[i]][,c('icgc_sample_id','gene_affected')]
  selected[[i]] <- unique(selected[[i]])
  colnames(selected[[i]]) <- c('sample_id','gene_affected')
  selected[[i]] <- as.data.frame(table(as.character(selected[[i]]$gene_affected)))
}

# differences in all protein-coding genes
gene_list <- subset(ens.prot$gene_id, ens.prot$gene_biotype == "protein_coding")
gene_list <- unique(gene_list)



samples <- list()
for (i in 1:length(selected)){
  samples[[i]] <- nrow(unique(as.data.frame(EC[[i]]$icgc_sample_id)))
  additional <- as.data.frame(subset(gene_list , !(gene_list %in% selected[[i]]$Var1)))
  additional$Freq <- 0
  names(additional)[1]<-paste("Var1")
  selected[[i]] <- rbind(selected[[i]],additional)
  selected[[i]]<-selected[[i]][order(as.character(selected[[i]]$Var1)),]
  selected[[i]]$Freq <- selected[[i]]$Freq/samples[[i]]
}

differences <- list()
index = 1
for (i in 1:5){
  k = i + 1
  for (j in k:5){
    if (index == 11) {
      next
    }
    differences[[index]] <- selected[[1]]
    differences[[index]] <- selected[[i]]$Freq - selected[[j]]$Freq
    index = index + 1
  }
}

differ <- list()
for (i in c(1:length(differences))){
  differ[[i]] <- data.frame(
    freq = differences[[i]],
    color = 0
  )
}

k = 1
t = length(EC)-1
for (i in c(1:t)){
  tt = i + 1
  for (j in c(tt:length(EC))){
    for (l in c(1:nrow(differ[[k]]))){
      if (differ[[k]][l,1] > 0){
        differ[[k]][l,2] <- boxplot_subtypes_colors[i]
      }
      else {
        differ[[k]][l,2] <- boxplot_subtypes_colors[j]
      }
    }
    k = k+1
    print(paste0("i is: ", i , " and j is: ", j))
  }
}


plot = list()
gg=1



for (o in c(1:length(differ))){
  data = data.frame(Significant_Genes = c(1:nrow(selected[[1]])), Frequency = differ[[o]]$freq)
  plot[[o]] = ggplot(data=data , aes(x=Significant_Genes, y=Frequency)) +
    geom_bar(stat="identity",width = 0.001, color=differ[[o]]$color)+
    ylim(-1, 1)+ylab("Differences of Mutation Rate")+xlab("Protein coding gene")+
    theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
          # Remove panel grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18.0),
          axis.title = element_text(size = 22.0)
    ) 
  
  # +geom_hline(yintercept=threshold_line[o], col="red", linetype="dotted")
  gg=gg+1
}
plot[[1]]
plot[[2]]
plot[[3]]
plot[[4]]
plot[[5]]
plot[[6]]
plot[[7]]
plot[[8]]
plot[[9]]
plot[[10]]
#20356
multi = ggarrange(plot[[2]], plot[[4]], plot[[7]], ncol = 1, nrow = 3)
annotate_figure(multi, top = text_grob("Figure 19 - Differences in relative frequency of mutations in all protein-coding gene"
                                       , color = "black", size = 36))
# print inch: 24,17,landscape



##### Differential analysis of significant features ####
selected <- list()
for (i in 1:length(EC)){
  selected[[i]] <- EC[[i]][,c('icgc_sample_id','gm')]
  selected[[i]] <- unique(selected[[i]])
  colnames(selected[[i]]) <- c('sample_id','gm')
  selected[[i]] <- as.data.frame(table(as.character(selected[[i]]$gm)))
}

# differences in all gene-motifs
gene_list <-  read.csv('Desktop/pancreatic cancer/Data/significant_features.csv')
gene_list <- unique(gene_list$Var1)

samples <- list()
for (i in 1:length(selected)){
  samples[[i]] <- nrow(unique(as.data.frame(EC[[i]]$icgc_sample_id)))
  selected[[i]] <- subset(selected[[i]] , selected[[i]]$Var1 %in% gene_list)
  additional <- as.data.frame(subset(gene_list , !(gene_list %in% selected[[i]]$Var1)))
  if (nrow(additional) == 0){
    next
  }
  additional$Freq <- 0
  names(additional)[1]<-paste("Var1")
  selected[[i]] <- rbind(selected[[i]],additional)
  selected[[i]]<- selected[[i]][order(as.character(selected[[i]]$Var1)),]
  selected[[i]]$Freq <- selected[[i]]$Freq/samples[[i]]
}

# samples <- list()
# for (i in 1:length(selected)){
#   samples[[i]] <- nrow(unique(as.data.frame(EC[[i]]$icgc_sample_id)))
#   additional <- as.data.frame(subset(gene_list , !(gene_list %in% selected[[i]]$Var1)))
#   if (nrow(additional) == 0){
#     next
#   }
#   additional$Freq <- 0
#   names(additional)[1]<-paste("Var1")
#   selected[[i]] <- rbind(selected[[i]],additional)
#   selected[[i]]<- selected[[i]][order(as.character(selected[[i]]$Var1)),]
#   selected[[i]]$Freq <- selected[[i]]$Freq/samples[[i]]
# }

differences <- list()
index = 1
for (i in 1:5){
  k = i + 1
  for (j in k:5){
    if (index == 11) {
      next
    }
    differences[[index]] <- selected[[1]]
    differences[[index]] <- selected[[i]]$Freq - selected[[j]]$Freq
    index = index + 1
  }
}

differ <- list()
for (i in c(1:length(differences))){
  differ[[i]] <- data.frame(
    freq = differences[[i]],
    color = 0
  )
}

k = 1
t = length(EC)-1
for (i in c(1:t)){
  tt = i + 1
  for (j in c(tt:length(EC))){
    for (l in c(1:nrow(differ[[k]]))){
      if (differ[[k]][l,1] > 0){
        differ[[k]][l,2] <- boxplot_subtypes_colors[i]
      }
      else {
        differ[[k]][l,2] <- boxplot_subtypes_colors[j]
      }
    }
    k = k+1
    print(paste0("i is: ", i , " and j is: ", j))
  }
}


plot = list(data)
gg=1

for (o in c(1:length(differ))){
  data = data.frame(Significant_Genes = c(1:nrow(selected[[1]])), Frequency = differ[[o]]$freq)
  plot[[o]] = ggplot(data=data , aes(x=Significant_Genes, y=Frequency)) +
    geom_bar(stat="identity",width = 0.001, color=differ[[o]]$color)+
    ylim(-1, 1)+ylab("Differences of Mutation Rate")+xlab("Significant features")+
    theme(plot.title = element_text(hjust = 0.01,size = 36),panel.border = element_blank(),  
          # Remove panel grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Remove panel background
          panel.background = element_blank(),
          # Add axis line
          axis.line = element_line(colour = "grey"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 18.0),
          axis.title = element_text(size = 22.0)
    ) 
  
  # +geom_hline(yintercept=threshold_line[o], col="red", linetype="dotted")
  gg=gg+1
}
plot[[1]]
plot[[2]]
plot[[3]]
plot[[4]]
plot[[5]]
plot[[6]]
plot[[7]]
plot[[8]]
plot[[9]]
plot[[10]]

multi = ggarrange(plot[[2]], plot[[4]], plot[[7]], ncol = 1, nrow = 3)

# multi = ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], plot[[5]], 
#                   plot[[6]], plot[[7]], plot[[8]], plot[[9]], plot[[10]], ncol = 1, nrow = 10)

annotate_figure(multi, top = text_grob("Figure 20 - Differences in relative frequency of mutations in gene-motifs"
                                       , color = "black", size = 36))
#482724










##### Fisher Exact test for significant features ####
EC[[5]] <- EC[[10]]
for (i in 6:14){
  EC[[6]] <- NULL
}


mcl <- list()
for (i in 1:5) {
  mcl[[i]] <- EC[[i]][,c("icgc_sample_id","gm")]
  mcl[[i]] <- subset(mcl[[i]], mcl[[i]]$gm %in% sig_g_sig_gm_01$V1)
  mcl[[i]] <- unique(mcl[[i]])
  colnames(mcl[[i]])<-c('sample_id','gm')
  mcl[[i]] <- as.data.frame(table(as.character(mcl[[i]]$gm)))
  additional <- as.data.frame(subset(sig_g_sig_gm_01$V1 , !(sig_g_sig_gm_01$V1 %in% mcl[[i]]$Var1)))
  if (nrow(additional) == 0){
    next
  }
  additional$Freq <- 0
  names(additional)[1]<-paste("Var1")
  mcl[[i]] <- rbind(mcl[[i]],additional)
  mcl[[i]] <- mcl[[i]][order(as.character(mcl[[i]]$Var1)),]
}



new_selected <- mcl
samples <- list()

for (i in 1:5) {
  samples[[i]] <- unique(as.data.frame(EC[[i]]$icgc_sample_id))
  all_ind <- c(1,2,3,4,5)
  new_selected[[i]]$non <- new_selected[[all_ind[-c(i)][1]]]$Freq + new_selected[[all_ind[-c(i)][2]]]$Freq + 
    new_selected[[all_ind[-c(i)][3]]]$Freq + new_selected[[all_ind[-c(i)][4]]]$Freq
  new_selected[[i]]$fisher<-0
  t<-as.matrix(new_selected[[i]][,2:4])
  this_cluster_sample_number = nrow(samples[[i]])
  others_cluster_sample_number = 758 - this_cluster_sample_number
  for (k in 1:length(sig_g_sig_gm_01$V1)) {
    # 70 is #cluster1 patients
    # 758 is all_patients - #cluster1 patients
    tt<-matrix(c(t[k,1], t[k,2], this_cluster_sample_number-t[k,1], 
                 others_cluster_sample_number-t[k,2]), nrow = 2)
    t[k,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
  }
  new_selected[[i]]$fisher <- t[,3]
}
saveRDS(new_selected, "Desktop/pancreatic cancer/significant_feature_fisher.rds")

##### Fisher Exact test for associated gene then gene-motif ####
EC[[5]] <- EC[[10]]
for (i in 6:14){
  EC[[6]] <- NULL
}


associated_genes <- readRDS("Desktop/pancreatic cancer/Data/associated_genes.rds")


mcl <- list()
for (i in 1:5) {
  mcl[[i]] <- EC[[i]][,c("icgc_sample_id","gene_affected","gm")]
  mcl[[i]] <- subset(mcl[[i]], mcl[[i]]$gene_affected %in% associated_genes[[i]]$gene)
  mcl[[i]] <- mcl[[i]][,c("icgc_sample_id","gm")]
  mcl[[i]] <- unique(mcl[[i]])
}
associated_gm <- rbind(mcl[[1]],mcl[[2]],mcl[[3]],mcl[[4]],mcl[[5]])
associated_gm <- as.data.frame(unique(associated_gm$gm))
colnames(associated_gm)<-c('gm')

mcl <- list()
for (i in 1:5) {
  mcl[[i]] <- EC[[i]][,c("icgc_sample_id","gene_affected","gm")]
  mcl[[i]] <- subset(mcl[[i]], mcl[[i]]$gene_affected %in% associated_genes[[i]]$gene)
  mcl[[i]] <- mcl[[i]][,c("icgc_sample_id","gm")]
  mcl[[i]] <- unique(mcl[[i]])
  colnames(mcl[[i]])<-c('sample_id','gm')
  mcl[[i]] <- as.data.frame(table(as.character(mcl[[i]]$gm)))
  additional <- as.data.frame(subset(associated_gm$gm , !(associated_gm$gm %in% mcl[[i]]$Var1)))
  if (nrow(additional) == 0){
    next
  }
  additional$Freq <- 0
  names(additional)[1]<-paste("Var1")
  mcl[[i]] <- rbind(mcl[[i]],additional)
  mcl[[i]] <- mcl[[i]][order(as.character(mcl[[i]]$Var1)),]
}



new_selected <- mcl
samples <- list()

for (i in 1:5) {
  samples[[i]] <- unique(as.data.frame(EC[[i]]$icgc_sample_id))
  all_ind <- c(1,2,3,4,5)
  new_selected[[i]]$non <- new_selected[[all_ind[-c(i)][1]]]$Freq + new_selected[[all_ind[-c(i)][2]]]$Freq + 
    new_selected[[all_ind[-c(i)][3]]]$Freq + new_selected[[all_ind[-c(i)][4]]]$Freq
  new_selected[[i]]$fisher<-0
  t<-as.matrix(new_selected[[i]][,2:4])
  this_cluster_sample_number = nrow(samples[[i]])
  others_cluster_sample_number = 758 - this_cluster_sample_number
  for (k in 1:length(sig_g_sig_gm_01$V1)) {
    # 70 is #cluster1 patients
    # 758 is all_patients - #cluster1 patients
    tt<-matrix(c(t[k,1], t[k,2], this_cluster_sample_number-t[k,1], 
                 others_cluster_sample_number-t[k,2]), nrow = 2)
    t[k,3] <- fisher.test(tt,alternative = "greater", conf.level = 0.99)$p.value
  }
  new_selected[[i]]$fisher <- t[,3]
}
saveRDS(new_selected, "Desktop/pancreatic cancer/associated_gm_fisher.rds")


##### significant difference in common gene for motif rate #### 

EC[[5]] <- EC[[10]]
for (i in 6:14){
  EC[[6]] <- NULL
}

associated_genes <- readRDS("Desktop/pancreatic cancer/Data/fisher_0.5.rds")
cosmic <- readRDS("Desktop/pancreatic cancer/Data/cosmic.rds")
for (i in c(1:length(associated_genes))){
  associated_genes[[i]] <- subset(associated_genes[[i]], associated_genes[[i]] %in% cosmic$converted_alias)
}
common_gene <- list()
index = 1
for (i in 1:5){
  k = i + 1
  for (j in k:5){
    if (index == 11) {
      next
    }
    # differences[[index]] <- selected[[1]]
    common_gene[[index]] <-subset(associated_genes[[i]], associated_genes[[i]] %in% associated_genes[[j]])
    index = index + 1
  }
}

# 2 -> (1,3)
# 3 -> (1,4)
# 4 -> (1,5)
# 8 -> (3,4)
# 9 -> (3,5)
# 10 -> (4,5)
all_motif = unique(protein2$motif)
{
  g_index = 1
  fisher_list <- list()
  z = 10
  for (gene in common_gene[[z]]) {
    i = 4
    j = 5
    
    
    
    selected <- list()
    
    selected[[i]] <- subset(EC[[i]], EC[[i]]$gene_affected == as.character(gene))
    selected[[i]] <- unique(selected[[i]][,c('icgc_sample_id','motif')])
    selected[[i]] <- as.data.frame(table(selected[[i]]$motif))
    if (length(selected[[i]]$Var1) != 96){
      missing_motif <- subset(all_motif, !(all_motif %in% selected[[i]]$Var1))
      t_df <- data.frame("Var1" = missing_motif, "Freq" = rep(0, length(missing_motif)))
      # colnames(t_df) <- c("Var1", "Freq")
      selected[[i]] <- rbind(selected[[i]], t_df)
      selected[[i]]$Var1 <- as.character(selected[[i]]$Var1)
      selected[[i]] <- selected[[i]][order(selected[[i]]$Var1),]
      remove(t_df)
    }
    i_number_of_samples = nrow(unique(as.data.frame(EC[[i]]$icgc_sample_id)))
    
    
    selected[[j]] <- subset(EC[[j]], EC[[j]]$gene_affected == as.character(gene))
    selected[[j]] <- unique(selected[[j]][,c('icgc_sample_id','motif')])
    selected[[j]] <- as.data.frame(table(selected[[j]]$motif))
    if (length(selected[[j]]$Var1) != 96){
      missing_motif <- subset(all_motif, !(all_motif %in% selected[[j]]$Var1))
      t_df <- data.frame("Var1" = missing_motif, "Freq" = rep(0, length(missing_motif)))
      # colnames(t_df) <- c("Var1", "Freq")
      selected[[j]] <- rbind(selected[[j]], t_df)
      selected[[j]]$Var1 <- as.character(selected[[j]]$Var1)
      selected[[j]] <- selected[[j]][order(selected[[j]]$Var1),]
      remove(t_df)
    }
    j_number_of_samples = nrow(unique(as.data.frame(EC[[j]]$icgc_sample_id)))
    
    
    fisher_list[[g_index]] <- selected[[i]]
    for (k in 1:length(all_motif)) {
      tt <- matrix(c(selected[[i]][k,2], selected[[j]][k,2], i_number_of_samples-selected[[i]][k,2], 
                     j_number_of_samples-selected[[j]][k,2]), nrow = 2)
      fisher_list[[g_index]][k,2] <- fisher.test(tt,alternative = "two.sided", conf.level = 0.99)$p.value
    }
    colnames(fisher_list[[g_index]]) <- c("motif","fisher")
    fisher_list[[g_index]]$gene <- gene
    fisher_list[[g_index]]$compare <- paste0('PCS',i,' vs PCS',j)
    
    g_index = g_index + 1
  }
  
  t <- list()
  for (k in 1:length(fisher_list)){
    t[[k]] <- subset(fisher_list[[k]], fisher_list[[k]]$fisher <= 0.05)
  }
  saveRDS(t,paste0('Desktop/pancreatic cancer/common_gene_motif_rate/PCS',i,' vs PCS',j,'.rds'))
}

{
  mf <- list()
  mf[[1]] <- readRDS('Desktop/pancreatic cancer/common_gene_motif_rate/PCS1 vs PCS3.rds')
  mf[[2]] <- readRDS('Desktop/pancreatic cancer/common_gene_motif_rate/PCS1 vs PCS4.rds')
  mf[[3]] <- readRDS('Desktop/pancreatic cancer/common_gene_motif_rate/PCS1 vs PCS5.rds')
  mf[[4]] <- readRDS('Desktop/pancreatic cancer/common_gene_motif_rate/PCS3 vs PCS4.rds')
  mf[[5]] <- readRDS('Desktop/pancreatic cancer/common_gene_motif_rate/PCS3 vs PCS5.rds')
  mf[[6]] <- readRDS('Desktop/pancreatic cancer/common_gene_motif_rate/PCS4 vs PCS5.rds')
  
  t <- mf[[1]][[1]]
  for(i in 1:length(mf)){
    t <- t[0,]
    for(j in 1:length(mf[[i]])){
      t <- rbind(t, mf[[i]][[j]])
    }
    write.csv(t, paste0('Desktop/pancreatic cancer/common_gene_motif_rate/',i,'.csv'))
  }
}





##### Permutation test ####

{
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(reshape2) 
}
L <- readRDS("Desktop/pancreatic cancer/Data/L.rds")
protein2 <- fread("Desktop/pancreatic cancer/Data/protein2.tsv", sep="\t")
ens.prot= fread(file="Desktop/pancreatic cancer/Data/ensemble_gene_types_19.tsv", sep="\t")
EC <- lapply(L, function(o){
  return(as.data.table(subset(protein2,protein2$icgc_sample_id %in% row.names(o))))
})
for(i in 1:length(EC)){
  EC[[i]]$gm = paste0(EC[[i]]$gene_affected,"-",EC[[i]]$motif)
}
EC[[5]] <- EC[[10]]
for (i in 6:14){
  EC[[6]] <- NULL
}


allpo <- list()
tmp <- data.table()
samples <- list()
tmp$gene <- subset(ens.prot$gene_id, ens.prot$gene_biotype == "protein_coding")
for (i in 1:5) {
  samples[[i]] <- unique(as.data.frame(EC[[i]]$icgc_sample_id))
  colnames(samples[[i]]) <- c("sample_id")
  allpo[[i]] <- data.table()
  for (j in c(1:length(samples[[i]]$sample_id))) {
    tmp$sample <- samples[[i]][j,1]
    allpo[[i]] <- rbind(allpo[[i]], tmp)
  }
}


mcl <- list()
sum_of_frequency <- list() 
for (i in 1:5) {
  mcl[[i]] <- EC[[i]][,c("icgc_sample_id","gene_affected")]
  mcl[[i]] <- unique(mcl[[i]])
  colnames(mcl[[i]])<-c('sample_id','gm')
  mcl[[i]] <- as.data.frame(table(as.character(mcl[[i]]$gm)))
  sum_of_frequency[[i]] <- sum(mcl[[i]]$Freq)
}


getFisherColumn <- function(allpo,sum_of_frequency,gene_list,i){
  
  new_selected <- lapply(1:length(allpo), function(i){
    indexes <- sample(1:nrow(allpo[[i]]), sum_of_frequency[[i]])
    u <- unique(allpo[[i]][indexes,])
    colnames(u)<-c('gm','sample_id')
    u <- as.data.frame(table(as.character(u$gm)))
    additional <- as.data.frame(subset(gene_list , !(gene_list %in% u$Var1)))
    if (nrow(additional) == 0){
      return(u)
    }
    additional$Freq <- 0
    names(additional)[1]<-paste("Var1")
    u <- rbind(u,additional)
    u <- u[order(as.character(u$Var1)),]
    return(u)
  })
  
  all_ind <- c(1,2,3,4,5)
  
  new_selected[[i]]$non <- new_selected[[all_ind[-c(i)][1]]]$Freq + new_selected[[all_ind[-c(i)][2]]]$Freq + 
    new_selected[[all_ind[-c(i)][3]]]$Freq + new_selected[[all_ind[-c(i)][4]]]$Freq
  new_selected[[i]]$fisher<-0
  t<-as.matrix(new_selected[[i]][,2:4])
  this_cluster_sample_number = nrow(samples[[i]])
  others_cluster_sample_number = 758 - this_cluster_sample_number
  for (k in 1:length(gene_list)) {
    # 70 is #cluster1 patients
    # 704 is all_patients - #cluster1 patients
    tt<-matrix(c(t[k,1], t[k,2], this_cluster_sample_number-t[k,1], 
                 others_cluster_sample_number-t[k,2]), nrow = 2)
    t[k,3] <- fisher.test(tt,alternative = "two.sided", conf.level = 0.99)$p.value
    # print(paste0('in gene list: item ',k,' and pcs ', i))
  }
  
  new_selected[[i]]$fisher <- t[,3]
  return(new_selected[[i]]$fisher)
  
}


library(doParallel)  
library(foreach)

#setup parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores[1]-1, outfile = "") #not to overload your computer
registerDoParallel(cl)



gene_list <- subset(ens.prot$gene_id, ens.prot$gene_biotype == "protein_coding")
all_fisher <- foreach(m=1:2, .combine='cbind') %:% 
  foreach (i=1:5) %dopar% {
    getFisherColumn(allpo,sum_of_frequency,gene_list,i)
    
  }
#stop cluster

stopCluster(cl)




for (i in c(6:length(all_fisher))){
  if (i%%5 == 0){
    all_fisher[[5]] <- cbind(all_fisher[[i]],all_fisher[[5]])
  }else{
    all_fisher[[i%%5]] <- cbind(all_fisher[[i]],all_fisher[[i%%5]])
  }
}
new_all_fisher <- list()
for (i in 1:5){
  new_all_fisher[[i]] <- all_fisher[[i]]
}

remove(all_fisher)




new_all_fisher <- readRDS("Desktop/pancreatic cancer/Data/new_all_fisher.rds")



for (i in 1:5){
  for (j in 1:nrow(all_fisher[[i]])){
    u <- sort(new_all_fisher[[i]][j,])
    if (all_fisher[[i]][j,4] > u[100]){
      all_fisher[[i]][j,4] <- 1
    }
  }
}

all_fisher_backup <- all_fisher

for(i in 1:5){
  all_fisher[[i]] <- subset(all_fisher[[i]], all_fisher[[i]]$fisher != 1)
}
all_fisher_0.01 <- all_fisher
all_fisher_0.001 <- all_fisher
saveRDS(all_fisher, "Desktop/pancreatic cancer/Data/fisher.rds") 

for (i in 1:5) {
  all_fisher[[i]] <- subset(all_fisher[[i]], all_fisher[[i]]$fisher < 0.05)
}


fisher <- readRDS("Desktop/pancreatic cancer/Data/fisher.rds")
samples <- list()
for (i in 1:5) {
  samples[[i]] <- nrow(unique(as.data.frame(EC[[i]]$icgc_sample_id)))
  fisher[[i]]$sample_mutation <- c(1:nrow(fisher[[i]]))
  for (j in 1:nrow(fisher[[i]])){
    mutation_rate <- unique(subset(EC[[i]]$icgc_sample_id, EC[[i]]$gene_affected == fisher[[i]][j,1]))
    fisher[[i]][j,5] <- length(mutation_rate)/samples[[i]]
  }
}
saveRDS(fisher, "Desktop/pancreatic cancer/Data/fisher.rds") 

fisher <- readRDS("Desktop/pancreatic cancer/Data/fisher.rds")
fisher_0.1 <- list()
fisher_0.2 <- list()
fisher_0.3 <- list()
fisher_0.4 <- list()
fisher_0.5 <- list()
fisher_0.6 <- list()
fisher_0.7 <- list()
fisher_0.8 <- list()
fisher_0.9 <- list()

for (i in 1:5){
  fisher_0.1[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.1)
  fisher_0.2[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.2)
  fisher_0.3[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.3)
  fisher_0.4[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.4)
  fisher_0.5[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.5)
  fisher_0.6[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.6)
  fisher_0.7[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.7)
  fisher_0.8[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.8)
  fisher_0.9[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.9)
}
saveRDS(fisher_0.5, "Desktop/pancreatic cancer/Data/fisher_0.5.rds")

fisher_0.5 <- readRDS("Desktop/pancreatic cancer/Data/fisher_0.5.rds")
library(venn)
venn(fisher_0.5, zcolor = "style")


for (i in 1:5){
  fisher_0.5[[i]] <- subset(fisher[[i]]$Var1, fisher[[i]]$sample_mutation >= 0.5)
}
saveRDS(fisher_0.5, "Desktop/pancreatic cancer/Data/fisher_0.5.rds")


##### reading expression data, cleaning and preparing meta data####

# pancreatic cancer samples input and cleaning
seq1=fread("../../data/expression/exp_seq.tsv.gz")
seq1=seq1[,-c("fold_change","submitted_sample_id","analysis_id","assembly_version",
              "platform","experimental_protocol","alignment_algorithm",
              "other_analysis_algorithm","reference_sample_type","raw_data_accession",
              "raw_data_repository","sequencing_strategy","gene_model","normalized_read_count",
              "normalization_algorithm","total_read_count","icgc_specimen_id")]

dss=merge(dss,clusters,"icgc_sample_id")
don=as.data.table(table(dss[,2]))[which(as.data.table(table(dss[,2]))[,2]==1)]
don=don[,-2]
colnames(don)=c("icgc_donor_id")
don=merge(don,dss[,c("icgc_donor_id","cluster")],by="icgc_donor_id")
###
seq2=subset(seq1,seq1$icgc_donor_id%in%don$icgc_donor_id | seq1$icgc_sample_id%in%dss$icgc_sample_id)
###
remove(seq1)
seq2=unique(seq2)
seq3=seq2[,c("icgc_sample_id","gene_id","raw_read_count")]

# removing duplicated values and genes
dup=duplicated(seq3[,1:2]) | duplicated(seq3[,1:2], fromLast = TRUE)  # finding duplicated rows
dup2=seq3[dup]  # finding duplicated rows in data
seq3=seq3[-which(seq3$gene_id%in%dup2$gene_id),] # removing duplicated genes in data
pseq=seq3[which(seq3$gene_id%in%ens.prot$gene_id)]
pseq2=dcast(pseq,formula = icgc_sample_id~gene_id,value.var = "raw_read_count")

pseq3=t(pseq2)
colnames(pseq3)=pseq3[1,]
pseq3=pseq3[-1,]

pseq4=apply(pseq3, 2, as.numeric)
rownames(pseq4)=rownames(pseq3)
pseq4=pseq4[,order(colnames(pseq4))]

pseq4=pseq4[-c(unique(which(is.na(pseq4), arr.ind=TRUE)[,1])),]

## meta data matrix
x=unique(seq2[,c("icgc_donor_id","icgc_sample_id")])
x1=merge(x,don,by="icgc_donor_id")
x2=merge(x[which(!(x$icgc_sample_id%in%x1$icgc_sample_id)),]
         ,dss[,c("icgc_sample_id","cluster")]
         ,by="icgc_sample_id")
x3=rbind(x1,x2)
x4=merge(x3,unique(seq2[,c("icgc_donor_id","project_code")]))
x4[,4]=as.factor(x4$project_code)

colD=x4[,c(2:4)]
colD=as.data.frame(colD)
rownames(colD)=as.character(colD$icgc_sample_id)
colD=colD[which(colD$cluster%in%c(1,2,3,4,10)),]
colD[which(colD$cluster==10),"cluster"]=5
colD$cluster=as.factor(as.character(colD$cluster))
colD=colD[,-1]
colnames(colD)=c("Subtypes","projects")

colD=colD[order(rownames(colD)),]
pseq4=pseq4[,order(colnames(pseq4))]
rownames(colD)=NULL

##### expression: one group VS remaning groups ####
# sorting samples in both matrices based on subtypes
colD=colD[order(rownames(colD)),]
pseq4=pseq4[,order(colnames(pseq4))]
colD$ind=seq(1,nrow(colD))
colD=colD[order(colD$Subtypes),]
pseq4=pseq4[,colD$ind]
colD=colD[,-3]

# data input and preprocessing
rna.o.ova=DESeqDataSetFromMatrix(pseq4,colD,design = ~0+Subtypes+projects)
rna.d.ova=DESeq(rna.o.ova,parallel = T)

# 1 vs all
res1=results(rna.d.ova,contrast=list(c("Subtypes1"), c("Subtypes2","Subtypes3","Subtypes4","Subtypes5")),
             listValues=c(1, -1/4),parallel = T,pAdjustMethod = "BH")
# res1$padj=p.adjust(res1$pvalue,method = "BH")
res1=res1[order(res1$padj,decreasing = F),]

# plotCounts(rna.d, gene=rownames(res1[11,]), intgroup=c("Subtypes"))
# plotCounts(rna.d, gene=rownames(res1[1,]), intgroup=c("projects","Subtypes"))

# 2 vs all
res2=results(rna.d.ova,contrast=list(c("Subtypes2"), c("Subtypes1","Subtypes3","Subtypes4","Subtypes5")),
             listValues=c(1, -1/4),parallel = T,pAdjustMethod = "BH")
# res2$padj=p.adjust(res2$pvalue,method = "BH")
res2=res2[order(res2$padj,decreasing = F),]

# plotCounts(rna.d, gene=rownames(res2[4,]), intgroup=c("Subtypes"))
# plotCounts(rna.d, gene=rownames(res2[1,]), intgroup=c("projects","Subtypes"))

# 3 vs all
res3=results(rna.d.ova,contrast=list(c("Subtypes3"), c("Subtypes1","Subtypes2","Subtypes4","Subtypes5")),
             listValues=c(1, -1/4),parallel = T,pAdjustMethod = "BH")
# res3$padj=p.adjust(res3$pvalue,method = "BH")
res3=res3[order(res3$padj,decreasing = F),]

# plotCounts(rna.d, gene=rownames(res3[11,]), intgroup=c("Subtypes"))
# plotCounts(rna.d, gene=rownames(res3[1,]), intgroup=c("projects","Subtypes"))

# 4 vs all
res4=results(rna.d.ova,contrast=list(c("Subtypes4"), c("Subtypes1","Subtypes2","Subtypes3","Subtypes5")),
             listValues=c(1, -1/4),parallel = T,pAdjustMethod = "BH")
# res4$padj=p.adjust(res4$pvalue,method = "BH")
res4=res4[order(res4$padj,decreasing = F),]

# plotCounts(rna.d, gene=rownames(res4[11,]), intgroup=c("Subtypes"))
# plotCounts(rna.d, gene=rownames(res4[1,]), intgroup=c("projects","Subtypes"))

# 5 vs all
res5=results(rna.d.ova,contrast=list(c("Subtypes5"), c("Subtypes1","Subtypes2","Subtypes3","Subtypes4")),
             listValues=c(1, -1/4),parallel = T,pAdjustMethod = "BH")
# res5$padj=p.adjust(res5$pvalue,method = "BH")
res5=res5[order(res5$padj,decreasing = F),]

# plotCounts(rna.d, gene=rownames(res5[1,]), intgroup=c("Subtypes"))
# plotCounts(rna.d.ova, gene=rownames(res5[1,]), intgroup=c("projects","Subtypes"))

###### unique genes for subtypes

spec.gene=list()
spec.gene[[1]]=rownames(res1[which(res1$padj<10^-2),])
spec.gene[[2]]=rownames(res2[which(res2$padj<10^-2),])
spec.gene[[3]]=rownames(res3[which(res3$padj<10^-2),])
spec.gene[[4]]=rownames(res4[which(res4$padj<10^-2),])
spec.gene[[5]]=rownames(res5[which(res5$padj<10^-2),])

venn(spec.gene,zcolor = "style",snames = c("PCS1","PCS2","PCS3","PCS4","PCS5"))

int=intersect(spec.gene[[1]],intersect(spec.gene[[2]],intersect(spec.gene[[3]],intersect(spec.gene[[4]],spec.gene[[5]]))))

uniwo1=union(spec.gene[[2]],union(spec.gene[[3]],union(spec.gene[[4]],spec.gene[[5]])))
uniwo2=union(spec.gene[[1]],union(spec.gene[[3]],union(spec.gene[[4]],spec.gene[[5]])))
uniwo3=union(spec.gene[[1]],union(spec.gene[[2]],union(spec.gene[[4]],spec.gene[[5]])))
uniwo4=union(spec.gene[[1]],union(spec.gene[[2]],union(spec.gene[[3]],spec.gene[[5]])))
uniwo5=union(spec.gene[[1]],union(spec.gene[[2]],union(spec.gene[[3]],spec.gene[[4]])))

s1=setdiff(spec.gene[[1]],uniwo1)
s2=setdiff(spec.gene[[2]],uniwo2)
s3=setdiff(spec.gene[[3]],uniwo3)
s4=setdiff(spec.gene[[4]],uniwo4)
s5=setdiff(spec.gene[[5]],uniwo5)

# UDEG with p-values
pv1 = subset(res1,rownames(res1)%in%s1)
pv2 = subset(res2,rownames(res2)%in%s2)
pv3 = subset(res3,rownames(res3)%in%s3)
pv4 = subset(res4,rownames(res4)%in%s4)
pv5 = subset(res5,rownames(res5)%in%s5)

pv1 = data.frame(x = rownames(pv1), y = pv1$padj)
pv2 = data.frame(x = rownames(pv2), y = pv2$padj)
pv3 = data.frame(x = rownames(pv3), y = pv3$padj)
pv4 = data.frame(x = rownames(pv4), y = pv4$padj)
pv5 = data.frame(x = rownames(pv5), y = pv5$padj)

# DEG with p-values
pv1 = res1[which(res1$padj<10^-2),]
pv2 = res2[which(res2$padj<10^-2),]
pv3 = res3[which(res3$padj<10^-2),]
pv4 = res4[which(res4$padj<10^-2),]
pv5 = res5[which(res5$padj<10^-2),]

pv1 = data.frame(x = rownames(pv1), y = pv1$padj)
pv2 = data.frame(x = rownames(pv2), y = pv2$padj)
pv3 = data.frame(x = rownames(pv3), y = pv3$padj)
pv4 = data.frame(x = rownames(pv4), y = pv4$padj)
pv5 = data.frame(x = rownames(pv5), y = pv5$padj)


##### Survival ####
aaa=meta.all
aaa[which(aaa$donor_vital_status=="deceased"),"donor_vital_status"]=as.numeric(2)
aaa[which(aaa$donor_vital_status=="alive"),"donor_vital_status"]=as.numeric(1)
aaa=as.data.table(aaa)

bbb=merge(aaa,clusters,by="icgc_sample_id")
bbb=bbb[which(bbb$cluster%in%c(1,2,3,4,10)&
                nchar(bbb$tumour_grade)!=0&
                nchar(bbb$tumour_stage)!=0),]

bbb[which(bbb$tumour_grade=="1 - Well differentiated"),"tumour_grade"]="Well differentiated"
bbb[which(bbb$tumour_grade=="2 - Moderately differentiated"),"tumour_grade"]="Moderately differentiated"
bbb[which(bbb$tumour_grade=="3 - Poorly differentiated"),"tumour_grade"]="Poorly differentiated"
bbb[which(bbb$tumour_grade=="4 - Undifferentiated"),"tumour_grade"]="Undifferentiated"
bbb[which(bbb$tumour_grade=="G1"),"tumour_grade"]="Well differentiated"
bbb[which(bbb$tumour_grade=="G2"),"tumour_grade"]="Moderately differentiated"
bbb[which(bbb$tumour_grade=="G3"),"tumour_grade"]="Poorly differentiated"
bbb[which(bbb$tumour_grade=="G4"),"tumour_grade"]="Undifferentiated"
bbb[which(bbb$tumour_grade=="G2 and G3"),"tumour_grade"]="Moderate to Poor"
bbb[which(bbb$tumour_grade=="Moderate to Poor"),"tumour_grade"]="Poorly differentiated"

bbb=bbb[-which(bbb$tumour_grade%in%c("G1 to G3","X - Cannot be assessed","Status Post Therapy"))]

table(bbb$tumour_grade[which(nchar(bbb$tumour_grade)!=0)])

bbb[which(bbb$tumour_stage=="T1N0M0"),"tumour_stage"]="IA"
bbb[which(bbb$tumour_stage=="T2N0M0"),"tumour_stage"]="IB"
bbb[which(bbb$tumour_stage=="T3N0M0"),"tumour_stage"]="IIA"
bbb[which(bbb$tumour_stage=="T2N1bM0"),"tumour_stage"]="IIB"
bbb[which(bbb$tumour_stage=="T2N1M0"),"tumour_stage"]="IIB"
bbb[which(bbb$tumour_stage=="T3N1M0"),"tumour_stage"]="IIB"
bbb[which(bbb$tumour_stage=="T3N1bM0"),"tumour_stage"]="IIB"
bbb[which(bbb$tumour_stage=="T3N1aM0"),"tumour_stage"]="IIB"
bbb[which(bbb$tumour_stage=="T3N3M0"),"tumour_stage"]="III"
bbb[which(bbb$tumour_stage=="T4N1M0"),"tumour_stage"]="III"
bbb[which(bbb$tumour_stage=="T4N0M0"),"tumour_stage"]="III"
bbb[which(bbb$tumour_stage=="T1N1M1"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="T4N1bM1"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="T2N0M1"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="T3N1bM1"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="T3N0M1"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="T2N1M1"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="T3N1M1"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="T3N1M2"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="T4N1M1"),"tumour_stage"]="IV"
bbb[which(bbb$tumour_stage=="TXNXM1"),"tumour_stage"]="IV"

bbb[-which(bbb$tumour_stage%in%c("IA","IB","IIA","IIB","III","IV")),"tumour_stage"]="X"

bbb[which(bbb$tumour_stage%in%c("IA","IB")),"tumour_stage"]="I"
bbb[which(bbb$tumour_stage%in%c("IIA","IIB")),"tumour_stage"]="II"

table(bbb$tumour_stage[which(nchar(bbb$tumour_stage)!=0)])


lrt.surv=function(mod.reduced,mod.full) {
  lrts=round(((-2)*(mod.reduced$loglik[2]- mod.full$loglik[2])),10)
  df=abs(summary(mod.reduced)$logtest[2]-summary(mod.full)$logtest[2])
  pvalue=round((1-pchisq(lrts,df)),20)
  return(pvalue)
}


s1=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$donor_age_at_diagnosis,data=bbb,na.action = "na.omit")

s2=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$tumour_stage,data=bbb,na.action = "na.omit")

s3=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$tumour_grade,data=bbb,na.action = "na.omit")

s4=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$clusters,data=bbb,na.action = "na.omit")

s5=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$donor_age_at_diagnosis+
           bbb$tumour_stage+
           bbb$tumour_grade+
           bbb$cluster,data=bbb,na.action = "na.omit")

s6=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$clusters+
           bbb$tumour_stage+
           bbb$cluster*bbb$tumour_stage
         ,data=bbb,na.action = "na.omit")

s7=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$cluster+
           bbb$tumour_grade+
           bbb$cluster*bbb$tumour_grade
         ,data=bbb,na.action = "na.omit")

s8=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$cluster+
           bbb$donor_age_at_diagnosis+
           bbb$cluster*bbb$donor_age_at_diagnosis
         ,data=bbb,na.action = "na.omit")

s9=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
           bbb$cluster+
           bbb$tumour_stage
         ,data=bbb,na.action = "na.omit")

s10=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
            bbb$cluster+
            bbb$tumour_grade
          ,data=bbb,na.action = "na.omit")

s11=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
            bbb$cluster+
            bbb$donor_age_at_diagnosis
          ,data=bbb,na.action = "na.omit")

s12=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
            bbb$donor_age_at_diagnosis+
            bbb$tumour_stage+
            bbb$tumour_grade,data=bbb,na.action = "na.omit")

s13=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
            bbb$donor_age_at_diagnosis+
            bbb$tumour_stage+
            bbb$cluster,data=bbb,na.action = "na.omit")

s14=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
            bbb$donor_age_at_diagnosis+
            bbb$tumour_grade+
            bbb$cluster,data=bbb,na.action = "na.omit")

s15=coxph(Surv(bbb$donor_survival_time,as.numeric(bbb$donor_vital_status==2))~
            bbb$tumour_grade+
            bbb$tumour_stage+
            bbb$cluster,data=bbb,na.action = "na.omit")



colnames(bbb)[11]="Subtype"
bbb=bbb[which(bbb$Subtype%in%c(1,2,3,4,10))]

fit <- survfit(Surv(time = (as.numeric(aaa$donor_survival_time/365)), 
                    event = (as.numeric(aaa$donor_vital_status))==2) ~ aaa$Subtype, data = aaa)
summary(fit, times = 5)

ggsurvplot(
  km_trt_fit,
  data = surv.data,
  risk.table = TRUE,
  pval = TRUE,
  xlab = "Time (Years)",
  break.time.by = 1,
  ggtheme = theme_classic(),
  risk.table.y.text.col = T,
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  pval.method = T,
  conf.int.style = "step",
  surv.median.line = "hv")