#required libraries
library(edgeR)
library(ggplot2)
library(colorspace)
library(limma)
library(forcats)


setwd("C:/Users/iwt/Dropbox (Duke Bio_Ea)/Benfey Lab/Current Lab Members/Isaiah_Taylor/manuscripts/KL_circumnutation/scripts_for_review/data/")
#rice gene annotations
#from:
#http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/

annotation=read.csv("all.locus_brief_info.7.csv", header = T)

#cytokinin response genes from rice roots. From Raines, Tracy, et al. "Characterization of the cytokinin-responsive transcriptome in rice." BMC plant biology 16.1 (2016): 260.
CRGs=read.csv("root_ctyo_response_genes.csv", header=T)

#Each sample is read in here. Samples starting with "WT" are wildtype. They are associated with the mutant
#box that was grown closest to them in the growth chamber (which is what is indicated by the number next
#in the names). This was tracked in case there was a position effect which would necessitate paired analysis.
#We didn't observe a strong position effect so analyzed the data as a randomized design.  The number "2"
#or "4" represents the 1-2mm section or 3-4 mm section, respectively.
x1 <- read.table("WT_179_2_counts", header = F)
x2 <- read.table("WT_179_4_counts", header = F)
x3 <- read.table("FN179_2_counts", header = F)
x4 <- read.table("FN179_4_counts", header = F)
x5 <- read.table("WT_287_2_counts", header = F)
x6 <- read.table("WT_287_4_counts", header = F)
x7 <- read.table("FN287_2_counts", header = F)
x8 <- read.table("FN287_4_counts", header = F)
x9 <- read.table("WT_790_2_counts", header = F)
x10 <- read.table("WT_790_4_counts", header = F)
x11 <- read.table("FN790_2_counts", header = F)
x12 <- read.table("FN790_4_counts", header = F)

#create single data frame of counts
raw_counts=cbind(x1$V2, x2$V2, x3$V2, x4$V2, x5$V2, x6$V2, x7$V2, x8$V2, x9$V2, x10$V2, x11$V2, x12$V2)
raw_counts=data.frame(raw_counts)

#label rows and columns
rownames(raw_counts) = x1$V1
colnames(raw_counts) = c("WT_179_2", "WT_179_4" 
                        ,"FN179_2", "FN179_4" 
                        ,"WT_287_2", "WT_287_4" 
                        ,"FN287_2", "FN287_4"
                        ,"WT_790_2", "WT_790_4" 
                        ,"FN790_2", "FN790_4") 

#remove ambiguously mapped or unmapped reads 
raw_counts=raw_counts[1:58061,]
raw_counts=raw_counts[!rownames(raw_counts) %in% c("no_feature"),]
raw_counts[is.na(raw_counts)]=0

#make DGElist
x <- DGEList(counts = raw_counts, genes = rownames(raw_counts))

#reads per library uniquely mapped to a gene:
#21201419 18474498 19783688 21343921 23122302 18485086 22868584 22729138 16519412 17247120 16504482 16459380

#make cpm and lcpm
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

#keep only genes that are expressed. "Expressed" here means counts observed in at least 3 samples
dim(x)
keep.exprs <- rowSums(cpm>0)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x) #compare to dim(x) above

#normalize data after removing low expressed genes
x <- calcNormFactors(x)

#cpm, lcpm of normalized values
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)


#account for factors in experiment
treatment=as.factor(rep(c("WT_2","WT_4","MUT_2","MUT_4"),3))

#assign factors to DGElist
x$samples$treatment=treatment

#make design matrix of experimental factors
design <- model.matrix(~0+treatment)
colnames(design) <- levels(treatment)
design #check matrix

#fit normalized dge list to model matrix
x <- estimateDisp(x, design)
fit <- glmQLFit(x, design)

#tr <- glmTreat(fit, coef=2, lfc=1)
#topTags(tr)

#making contrast matrix for tests of interest
my.contrasts <- makeContrasts(WT2_v_MUT2=WT_2-MUT_2, WT4_v_MUT4=WT_4-MUT_4, levels=design)

#QLF tests 
WT2_v_MUT2 <- glmQLFTest(fit, contrast=my.contrasts[,"WT2_v_MUT2"]) 
WT4_v_MUT4 <- glmQLFTest(fit, contrast=my.contrasts[,"WT4_v_MUT4"]) 

#extract results. toptags adds FDR column.
#section 2
sec2 <- topTags(WT2_v_MUT2, n = "Inf")$table
sec2 <- subset(sec2, sec2$FDR<0.05)
sec2 <- subset(sec2, abs(sec2$logFC)>0.4)
sec2$WT_av <- apply(cpm[sec2$genes,c(1,5,9)], 1,mean)
sec2$MUT_av <- apply(cpm[sec2$genes,c(3,7,11)], 1,mean)
sec2 = cbind(sec2, cpm[sec2$genes,c(1,5,9,3,7,11)])
WT2_low_MUT2_high=sec2[sec2$logFC<0,]
WT2_high_MUT2_low=sec2[sec2$logFC>0,]
rownames(WT2_high_MUT2_low)=WT2_high_MUT2_low$genes
rownames(WT2_low_MUT2_high)=WT2_low_MUT2_high$genes


#section 4
sec4 <- topTags(WT4_v_MUT4, n = "Inf")$table
sec4 <- subset(sec4, sec4$FDR<0.05)
sec4 <- subset(sec4, abs(sec4$logFC)>0.4)
sec4$WT_av <- apply(cpm[sec4$genes,c(2,6,10)], 1,mean)
sec4$MUT_av <- apply(cpm[sec4$genes,c(4,8,12)], 1,mean)
sec4 = cbind(sec4, cpm[sec4$genes,c(2,6,10,4,8,12)])
WT4_low_MUT4_high=sec4[sec4$logFC<0,]
WT4_high_MUT4_low=sec4[sec4$logFC>0,]
rownames(WT4_high_MUT4_low)=WT4_high_MUT4_low$genes
rownames(WT4_low_MUT4_high)=WT4_low_MUT4_high$genes

all_genes=cpm(x)

#write to file
write.table(WT2_high_MUT2_low, file="WT2_high_MUT2_low")
write.table(WT2_low_MUT2_high, file="WT2_low_MUT2_high")
write.table(WT4_high_MUT4_low, file="WT4_high_MUT4_low")
write.table(WT4_low_MUT4_high, file="WT4_low_MUT4_high")
write.table(all_genes, file="all_genes")



#enrichment analysis

#keep only genes that are expressed in the 1-2mm sectoin
keep.exprs <- rowSums(cpm[,c(1,3,5,7,9,11)]>0)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x) #compare to dim(x) above

rownames(CRGs)=CRGs$gene

#up_genes are genes with positive logFC from rice root cytokinin treatment data
up_genes=data.frame(CRGs[CRGs$log2FC>0,]$gene)

#down_genes are genes with negative logFC from rice root cytokinin treatment data
down_genes=data.frame(CRGs[CRGs$log2FC<0,]$gene)

#we filter out genes which are counted as "non-expressed" during our filtering above
up_genes_expressed=intersect(up_genes[,1],x$genes$genes)
down_genes_expressed=intersect(down_genes[,1],x$genes$genes)

#here are the genes which are in the two tested lists
up_in_cyto_down_in_mut=intersect(up_genes_expressed, WT2_high_MUT2_low$genes)
down_in_cyto_up_in_mut=intersect(down_genes_expressed, WT2_low_MUT2_high$genes)
down_in_cyto_down_in_mut=intersect(down_genes_expressed, WT2_high_MUT2_low$genes)
up_in_cyto_up_in_mut=intersect(up_genes_expressed, WT2_low_MUT2_high$genes)

#enrichment of cytokinin up genes in mutant down
1-phyper(length(up_in_cyto_down_in_mut),length(up_genes_expressed),dim(x)[1]-length(up_genes_expressed),nrow(WT2_high_MUT2_low))

#enrichment of cytokinin down genes in mutant up
1-phyper(length(down_in_cyto_up_in_mut),length(down_genes_expressed),dim(x)[1]-length(down_genes_expressed),nrow(WT2_low_MUT2_high))

#enrichment of cytokinin up genes in mutant up
1-phyper(length(up_in_cyto_up_in_mut),length(up_genes_expressed),dim(x)[1]-length(up_genes_expressed),nrow(WT2_low_MUT2_high))

#enrichment of cytokinin down genes in mutant down
1-phyper(length(down_in_cyto_down_in_mut),length(down_genes_expressed),dim(x)[1]-length(down_genes_expressed),nrow(WT2_high_MUT2_low))


#plot observed versus expected for the two comparisons of interest
enrichment=data.frame(matrix(ncol=2, nrow=2))
colnames(enrichment)=c("expected","observed")
rownames(enrichment)=c("Cytokinin induced/Higher in wildtype than in mutant", "Cytokinin repressed/Lower in wildtype than in mutant")
enrichment[] = c(nrow(WT2_high_MUT2_low)*length(up_genes_expressed)/dim(x)[1], nrow(WT2_low_MUT2_high)*length(down_genes_expressed)/dim(x)[1],length(up_in_cyto_down_in_mut),length(down_in_cyto_up_in_mut))
barplot(height=unlist(enrichment[1,]), ylim=c(0,200), main=wrapper(paste(rownames(enrichment)[1],", p-value = 0.0"),25))
barplot(height=unlist(enrichment[2,]), ylim=c(0,27.5), main=wrapper(paste(rownames(enrichment)[2],", p-value = 1.925917e-08"),25))

