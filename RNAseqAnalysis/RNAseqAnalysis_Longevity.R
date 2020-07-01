### RNAseq analysis for Longevity project

library("dplyr")
library("DESeq2")
library("apeglm")
library("sjPlot")
library("gridExtra")

###
### Differential expression analysis in DESEQ2
###


#load raw gene counts
  #10915 expressed genes with mean CPM>1
  #162 samples
  counts <- read.table("../CleanCode/ReadCounts_forDeseq.txt",h=T)

#load metadata
  metadata <- read.table("../CleanCode/Metadata_forDeseq.txt",h=T, row.names = 1)
  metadata$Plate <- as.factor(metadata$Plate)
  metadata$Age  <- as.factor(metadata$Age)
  
#create DESEQ object
  long <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = metadata,
                                 design =  ~ Plate + Age)
  
#variance stability the counts
  vs.long <- vst(long,blind = F)
  
#Set "young" flies as reference level
  long$Age = relevel( long$Age, "young")
  
#Differential expression analysis
  dea <- DESeq(long, test = "Wald")
  
  #Age effect on Control condition (null lfc=0, FDR 5%)
  res.control <- results(dea, name = "Age_oldC_vs_young", alpha = 0.05, lfcThreshold = 0)
  res.control <- lfcShrink(dea, coef = "Age_oldC_vs_young", res = res.control, type = "apeglm")
  res.control.sig <- subset(res.control, res.control$padj<0.05)

  #Age effect on High-sugar condition (null lfc=0, FDR 5%)
  #extracting young vs old HIGH SUGAR results
  res.hs <- results(dea, name = "Age_oldHS_vs_young", alpha = 0.05, lfcThreshold = 0)
  res.hs <- lfcShrink(dea, coef = "Age_oldHS_vs_young", res = res.hs, type = "apeglm")
  res.hs.sig <- subset(res.hs, res.hs$padj<0.05)

#Classify genes according to the significance group (hs_ony, ctrl_only, both)
  res.hs <- res.hs[order(match(rownames(res.hs), rownames(res.control))),]
  all_equal(rownames(res.hs), rownames(res.control))
  res.all <- cbind(res.control[,c(2,4,5)], res.hs[c(2,4,5)])
  colnames(res.all) <- c("lfc_ctrl", "p_ctrl", "FDR_ctrl", "lfc_hs", "p_hs", "FDR_hs" )
  res.all <- as.data.frame(res.all)
  
  #Assigment to signficance groups follows the criteria used in AlleleFrequencyChanges analysis
  siggroup <- function(q){
    if(is.na(q[3]) | is.na(q[6])) {
      "NA"
    } else {
      if(q[3]<0.05 & q[5]>0.05){
        "ctrl"  
      } else {
        if(q[2]>0.05 & q[6]<0.05){
          "hs"
        } else {
          if(q[3]<0.05 & q[6]<0.05){
            "both"
          } else { 
            "none" 
            }}}}}          
  res.all$diffexp_in <- apply(res.all,1, function(x) siggroup(x))
  
  #Add information about lfc>1
  res.all$lfc_1 <-apply(res.all, 1, function(x) if(abs(as.numeric(x[1]))>1 | abs(as.numeric(x[4]))>1) {"yes"} else {"no"})
  
  #write.table(res.all, "Results_RNAseqAge_longevity.txt", quote=F, sep="\t", col.names = T, row.names = T)
  
###
### PLOTS 
###
  
#Volcano plot for High-sugar and Control results  
  volcano <- function(x){
      tab = data.frame(logFC = x$log2FoldChange, negLogPval = -log10(x$pvalue))
      par(mar = c(5, 4, 4, 4))
      plot(tab, pch = 1, cex = 0.7, xlab = "Log Fold Change (LFC)", ylab = expression(~-log[10] (pval)), col="darkgrey", cex.lab=.8,
           main="" )
      lfc = 0.0
      pval = 0.05
      signGenes = (abs(tab$logFC) > lfc & -log10(x$padj) > -log10(pval))
      points(tab[signGenes, ], pch = 1, cex = 0.7, col="darkorange") 
      }
    
  par(mfrow=c(1,2))  
  volcano(res.control) ; title(main="Expression changes due to age in Control diet")
  volcano(res.hs); title(main="Expression changes due to age in High-sugar diet")
  

#Comparison of effect size in each category (hs_only, ctrl_only, both)
  par(mfrow=c(2,2))  
  {
    res.all.both <- subset(res.all, res.all$diffexp_in=="both")
    res.all.hs <- subset(res.all, res.all$diffexp_in=="hs")
    res.all.c <- subset(res.all, res.all$diffexp_in=="ctrl")
   
    plot(res.all$lfc_ctrl~res.all$lfc_hs, xlab="LFC High-sugar", ylab="LFC Control",
         main="", type="n")
    points(res.all.hs$lfc_ctrl ~ res.all.hs$lfc_hs, col="salmon")
    points(res.all.c$lfc_ctrl ~ res.all.c$lfc_hs, col="orange")
    points(res.all.both$lfc_ctrl ~ res.all.both$lfc_hs, col="blue")
    abline(h=0, lty=2); abline(v=0, lty=2)
    legend(-5, 4, c("Both","HS_only", "Ctrl_only"), bty = "n", fill=c("blue", "salmon", "orange"), border = F)
  }
  
  #The effect size of genes differentially expressed in both diets tends to be larger than diet-specific effects
  #pval 2.2e-16, odds ratio 3.24
  # 1417 = all both (2098) - lfc1 both (681)
  # 3375 = all non both (3875) - lfc1 non both (500)
  fisher.test(matrix(c(681, 1417, 500,  3375 ), nrow = 2))
  
#Gene expression of the genes used in the Validation Experiment
  {  
  #Using normalized counts 
  #The bimodal distribution is due to the plate the samples were processed in ("Plate" effect)
  par(mfrow=c(3,1))
  validation <- c("FBgn0266347", "FBgn0267429", "FBgn0004797")
  for (i in 1:length(validation)){
    plotCounts(long, gene=validation[i], intgroup="Age", main=validation[i], xlab="Age group") 
  }
  
  #Ploting the marginal effects after accounting for the batch effect
  par(mfrow=c(2,2))
  validation <- c("FBgn0266347", "FBgn0267429", "FBgn0004797")
  a=list()
  for (i in 1:length(validation)){
    tmp <- lm(counts(dea,normalized=T)[rownames(counts(dea))==validation[i]] ~ Plate + Age, data = colData(long))
    a[[i]]<- plot_model(tmp, type="pred", terms="Age", axis.title = c("Age group","Normalized Gene counts"), title=validation[i])
  }
  do.call(gridExtra::grid.arrange, a)
  }

#PCA of the RNAseq samples
  {
  #subset of most variable genes
  par(mfrow=c(2,2))
  topgenes <- plotPCA(vs.long, intgroup=c("Age"), ntop=500, returnData=T)
  plot(topgenes$PC1, topgenes$PC2, xlab="PC1 87%", ylab="PC2 3%", col="blue", pch=16, main="PCA with 500 most variable genes")
  points(topgenes$PC1[topgenes$group=="young"], topgenes$PC2[topgenes$group=="young"], col="salmon", pch=16)
  points(topgenes$PC1[topgenes$group=="oldC"], topgenes$PC2[topgenes$group=="oldC"], col="orange", pch=16)
  legend(-65,-15, fill=c("blue","salmon","orange"), c("Old Hs", "Young", "Old Ctrl"), bty="n", border = F)
  
  #using all expressed genes
  allgenes <- plotPCA(vs.long, intgroup=c("Age"), ntop=nrow(counts(long)), returnData=T)
  plot(allgenes$PC1, allgenes$PC2, xlab="PC1 87%", ylab="PC2 5%", col="blue", pch=16, main="PCA with 500 most variable genes")
  points(allgenes$PC1[allgenes$group=="young"], allgenes$PC2[allgenes$group=="young"], col="salmon", pch=16)
  points(allgenes$PC1[allgenes$group=="oldC"], allgenes$PC2[allgenes$group=="oldC"], col="orange", pch=16)
  legend(-30,-10, fill=c("blue","salmon","orange"), c("Old Hs", "Young", "Old Ctrl"), bty="n", border = F)
  }
  

###
### OVERLAP WITH AGE-RELATED SNP EFFECTS
###
  
#Load gene-level results from Allele Frequency Changes analysis
#The overlap between significant SNPs and genes was done by extending gene body +-1kb to include regulatory regions
  afc <- read.table("Results_AFC_genelevel.txt",h=T)

#Load Differential Expression results
  dge <- read.table("Results_RNAseqAge_longevity.txt",h=T)
  dge$gene <- rownames(dge)

#Overlap between DGE and AFC genes
  length(intersect(afc$gene[afc$sig==1], dge$gene[dge$diffexp_in!="none"]))
  l<- list(as.character(afc$gene[afc$sig==1]), as.character(dge$gene[dge$diffexp_in!="none"]))
  VennDiagram::venn.diagram(x=l,category.names = c("AFC","DEG"), print.mode = c("raw","percent"), filename = "Venn.AFC.DEG.png",fill=c("salmon","orange"))

#Significant AFC genes are enriched for DE genes (pval 3.5e-7, odds ratio 1.3)
  a = length(intersect(afc$gene[afc$sig==1], dge$gene[dge$diffexp_in!="none"])) #Number of significant AFC genes that are DE genes
  b = sum(afc$sig==1) - a #Number of significant AFC genes that are not DE genes
  c = length(intersect(afc$gene[afc$sig==0], dge$gene[dge$diffexp_in!="none"]))#Number of non significant AFC genes that are DE genes
  d = sum(afc$sig==0) - c #Number of non significant AFC that are not DE genes
  fisher.test(matrix(c(a,b,c,d), nrow = 2))
  
  
#Overlap between AFC genes and large effect size DGE (lfc>1)
  length(intersect(afc$gene[afc$sig==1], dge$gene[dge$diffexp_in!="none" & dge$lfc_1=="yes"]))
  l<- list(as.character(afc$gene[afc$sig==1]), as.character(dge$gene[dge$diffexp_in!="none" & dge$lfc_1=="yes"]))
  VennDiagram::venn.diagram(x=l,category.names = c("AFC","DEG lfc>1"), print.mode = c("raw","percent"), filename = "Venn.AFC.DEGlfc1.png",fill=c("salmon","orange"),
                            cat.pos=c(0,0))

#Significant AFC genes are not enriched for DE genes of large effect (lfc>1)
  a = length(intersect(afc$gene[afc$sig==1], dge$gene[dge$diffexp_in!="none" & dge$lfc_1=="yes"])) #Number of significant AFC genes that are lfc>1 DE genes
  b = length(intersect(afc$gene[afc$sig==1], dge$gene[dge$diffexp_in!="none" & dge$lfc_1=="no"]))  #Number of significant AFC genes that are lfc<1 DE genes
  c = length(intersect(afc$gene[afc$sig==0], dge$gene[dge$diffexp_in!="none" & dge$lfc_1=="yes"])) #Number of non significant AFC genes that are lfc>1 DE genes
  d = length(intersect(afc$gene[afc$sig==0], dge$gene[dge$diffexp_in!="none" & dge$lfc_1=="no"]))  #Number of non significant AFC that are lfc<1 DE genes
  fisher.test(matrix(c(a,b,c,d), nrow = 2))    
  

#34% of the overlap between AFC and DE genes have effects in the same group (both diets, hs only, ctrl only) 
  afc_dge <- merge(afc[afc$sig==1,], dge[dge$diffexp_in!="none",], by="gene")
  afc_dge$overlap<- apply(afc_dge, 1, function(x) if(x[13]=="ctrl" && x[4]==1 | 
                                                     x[13]=="hs"   && x[3]==1 |
                                                     x[13]=="both" && x[2]==1 ) {
    "samegroup" } else { "diffgroup" })
  table(afc_dge$overlap) 
 
   
#Non exonic Age-related SNPs are not enriched for Differentially Expressed genes compared to exonic SNPs
  #using genes that are tagged by just one type of SNP (either exonic or not)
  a = length(intersect(afc$gene[afc$sig==1 & afc$snptype=="noexon"], dge$gene[dge$diffexp_in!="none"])) # Number of DGE genes tagged by non exonic AFC sig SNPs
  b = length(intersect(afc$gene[afc$sig==1 & afc$snptype=="noexon"], dge$gene[dge$diffexp_in=="none"])) # Number of not DGE genes tagged by non exonic AFC sig SNPs
  c = length(intersect(afc$gene[afc$sig==1 & afc$snptype=="exon"], dge$gene[dge$diffexp_in!="none"]))   # Number of DGE genes tagged by exonic AFC sig SNPs
  d = length(intersect(afc$gene[afc$sig==1 & afc$snptype=="exon"], dge$gene[dge$diffexp_in=="none"]))   # Number of not DGE genes tagged by exonic AFC sig SNPs
  fisher.test(matrix(c(a,b,c,d), nrow = 2))    
  
  #including genes that are tagged by both SNP types
  a = length(intersect(afc$gene[afc$sig==1 & afc$signotexon==1], dge$gene[dge$diffexp_in!="none"])) # Number of DGE genes tagged by non exonic AFC sig SNPs
  b = length(intersect(afc$gene[afc$sig==1 & afc$signotexon==1], dge$gene[dge$diffexp_in=="none"])) # Number of not DGE genes tagged by non exonic AFC sig SNPs
  c = length(intersect(afc$gene[afc$sig==1 & afc$sigexon==1], dge$gene[dge$diffexp_in!="none"]))   # Number of DGE genes tagged by exonic AFC sig SNPs
  d = length(intersect(afc$gene[afc$sig==1 & afc$sigexon==1], dge$gene[dge$diffexp_in=="none"]))   # Number of not DGE genes tagged by exonic AFC sig SNPs
  fisher.test(matrix(c(a,b,c,d), nrow = 2))    
  
  #If AFC effects are restricted to GxE (diet specific effect), the effect is also not-signifcant 
  a = length(intersect(afc$gene[afc$sig==1 & afc$both==0 & afc$snptype=="noexon"], dge$gene[dge$diffexp_in!="none"])) # Number of DGE genes tagged by non exonic GxE AFC sig SNPs
  b = length(intersect(afc$gene[afc$sig==1 & afc$both==0 & afc$snptype=="noexon"], dge$gene[dge$diffexp_in=="none"])) # Number of not DGE genes tagged by non exonic GxE AFC sig SNPs
  c = length(intersect(afc$gene[afc$sig==1 & afc$both==0 & afc$snptype=="exon"], dge$gene[dge$diffexp_in!="none"]))   # Number of DGE genes tagged by exonic GxE AFC sig SNPs
  d = length(intersect(afc$gene[afc$sig==1 & afc$both==0 & afc$snptype=="exon"], dge$gene[dge$diffexp_in=="none"]))   # Number of not DGE genes tagged by exonic GxE AFC sig SNPs
  fisher.test(matrix(c(a,b,c,d), nrow = 2))    
  

  
  
