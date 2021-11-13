library(data.table)
library(plyr)
library(dplyr)
library(igraph)
library(lsa)
library(ggplot2)
library(tidyr)
library(corrr)
rm(list = ls())
# Learn where data files are stored.
source("code/dataConfig.R")
#load utility functions
source("code/utilityFunctions.R")
# We compute similarity of two pathways by using different gene similarity measures

simFiles  = c(proDomainSimFile,
              subcellularSimFile,molecularSimFile,bioProcessSimFile
)


# GeneSim computes similarity of two genes by analyzing the pathways they co-appear
# The similaritySourceFile has pathways and genes in them.
geneSim <- function(similaritySourceFile, testing = FALSE) {
  dat <- read.csv(
    similaritySourceFile,
    header = FALSE,
    sep = ",",
    stringsAsFactors = F
  )
  colnames(dat) <- c("path", "genelist")
  # check if the pathway file is clean.
  if(length(dat$path)!=length(unique(dat$path))){
    message("Some pathways appear multiple times in ",similaritySourceFile)
    return(-1);
  }
  # give ids to pathways
  numPathways <- nrow(dat)
  pathwayId <- transform(dat, id = seq_len(numPathways))
  # how many genes are there?
  geneNames <- sort(unique(unlist(as.list(strsplit(
    dat$genelist, "\t"
  )))))
  numberOfGenes = length(geneNames)
  geneId<-transform(data.frame(geneNames), id = seq_len(numberOfGenes))
  colnames(geneId) <- c("gene", "id")
  
  # how many pathways?
  if (testing) {
    message("There are ",numberOfGenes, " genes, ", numPathways, " pathways.")
  }
  
  # we will create a new 1-to-1 mapping gene->pathway
  # TRhis code snippet takes too much time.
  newDat <- data.frame(gene = NA, pathway = NA)
  pathwayColumn <- numeric(length=nrow(pathwayId))
  geneColumn <- numeric(length=nrow(pathwayId))
  i <- 0
  for (row in 1:nrow(pathwayId)) {
    thisPathway <- pathwayId[row, "id"]
    for (geneName in unlist(as.list(strsplit(pathwayId[row, "genelist"], "\t")))) {
      i = i + 1
      pathwayColumn[i] <- thisPathway
      geneColumn[i] <- geneId[geneId$gene==geneName,"id"]
    }
  }
  newDat <- data.frame(geneColumn, pathwayColumn, stringsAsFactors = FALSE)
  colnames(newDat) <- c("gene", "pathway")
  # we will cache all pathways of a gene in a map to speed up future comps. 
  pathwaysOfGenes <- list(length=numberOfGenes)
  for (currentGene1 in 1:numberOfGenes) {
    pathwaysOfGenes[[currentGene1]] <- newDat[newDat$gene == currentGene1,]$pathway
  }
  
  # We reduce the number of sim computations by choosing to compute the similarity of genes in a pathway.
  ma<-matrix(0,ncol=numberOfGenes,nrow=numberOfGenes)
  for(pathway1Id in unique(newDat$pathway)){
    genesOfPathway1= newDat[newDat$pathway==pathway1Id,"gene"]
    numGenesInPathway=length(genesOfPathway1)       
    for (g1index in 1:(numGenesInPathway-1)) {
      currentGene1=genesOfPathway1[[g1index]]
      vec1 = pathwaysOfGenes[[currentGene1]]
      for (g2index in (g1index + 1):numGenesInPathway) {
        currentGene2=genesOfPathway1[[g2index]]
        ma[currentGene1,currentGene2]=ma[currentGene2,currentGene1]=1
      }
    }
  }
  
  # we now start computing the pairwise similarities of genes
  geneColumn1 <- character()
  geneColumn2 <- character()
  simColumn <- numeric()
  index = 0
  for (currentGene1 in 1:(nrow(ma) - 1)) {
    vec1 = pathwaysOfGenes[[currentGene1]]
    gene1Name = (geneId[geneId$id == currentGene1, "gene"])
    for (currentGene2 in (currentGene1 + 1):nrow(ma)) {
      if (ma[currentGene1, currentGene2] == 1) {
        vec2 = pathwaysOfGenes[[currentGene2]]
        sim = cosineSim(vec1, vec2)
        if (sim > 0) {
          gene2Name = (geneId[geneId$id == currentGene2, "gene"])
          index = index + 1
          geneColumn1[index] <- gene1Name
          geneColumn2[index] <- gene2Name
          simColumn[index] <- sim
        }
      }
    }
  }
  return(data.frame(geneColumn1, geneColumn2, simColumn, stringsAsFactors = FALSE))
}


loadSim<-function(dataDir,similaritySourceFile){
  message(similaritySourceFile, " is being processed to create gene similarity values.")
  fname=gsub("[^[:alnum:][:space:]]","",similaritySourceFile)
  
  subDir <- "simcache"
  ifelse(!dir.exists(file.path(dataDir, subDir)), dir.create(file.path(dataDir, subDir)), FALSE)
  resFile <-  paste0(dataDir,"simcache/","geneSim",
                     fname,
                     ".rds"
  )
  result = ""
  tryCatch({
    result <-readRDS(file = resFile)
    message("Previously computed similarity results were found and loaded.")
  }, error = function(e) {
    message("There is no previous save. Gene similarity will be computed once and stored in the simcache subfolder for future use.")
    result = geneSim(similaritySourceFile)
    saveRDS(result,
            file = resFile)
  }, finally = {
    
  })
  return(result)
}


useIdenticalGenes = F
pathways = read.csv(pathwayFile)
load(bioPlanetsRdata)

for(nearestNeighborCount in c(5,10,20,50,100,200)){
  results<-data.frame()
  
  for(simFile in simFiles){
    geneSimMap = as.data.table(loadSim(dataDir,simFile))
    colnames(geneSimMap)<-c("gene1","gene2","similarity")
    
    message("Similarity reference file: ",simFile,", k-neighbor:",nearestNeighborCount)
    
    for(graphIndex in 1:length(pathway.net.list)){
      gr=graph_from_edgelist(as.matrix(pathway.net.list[[graphIndex]][,c("source","target")]))
      # graph is here. 
      # How do we know that this specific crosstalk pathway network is any good?
      # We can quantify edge similarity. If two pathways are connected, 
      # we can compute their similarity by taking the average similarity of their genes.
      # assumption: edges between highly similar pathways are good
      sumOfAllPathwaySim=0.0
      simBagAll<-c()
      
      for(pathway1Name in V(gr)$name){
        pathway1Genes <- pathways[pathways$PATHWAY_NAME==pathway1Name,]$GENE_SYMBOL
        sumOfPathway1Sim<-0.0
        neList <- neighbors(gr,pathway1Name)$name
        df10<-geneSimMap[geneSimMap$gene1%in%pathway1Genes,]
        df20<-geneSimMap[geneSimMap$gene2%in%pathway1Genes,]
        simBag<-c()
        for (pathway2Name in neList){
          # todo: Bioplanet has gene symbols that are 9818 Levels: 1-Dec 1-Sep 2-Sep 4-Sep 5-Sep 6-Mar 7-Sep 9-Sep A1BG A1CF A2M ... ZYX
          pathway2Genes = as.character(pathways[pathways$PATHWAY_NAME==pathway2Name,]$GENE_SYMBOL)
          df1<-df10[df10$gene2%in%pathway2Genes,]
          df2<-df20[df20$gene1%in%pathway2Genes,]
          both = length(intersect(pathway1Genes,pathway2Genes))
          
          
          if(useIdenticalGenes){
            N2 = nearestNeighborCount-both
            if(N2>0){
              sumOfSim = sum(sort(union(df1$similarity,df2$similarity),decreasing=T)[1:N2],rm.na=T)
              avgSimOfPathway1to2 = (both+sumOfSim)/(nearestNeighborCount+both)
            }else {
              avgSimOfPathway1to2 = 1.0
            }
          }
          else{
            avgSimOfPathway1to2 = sum(sort(union(df1$similarity,df2$similarity),decreasing=T)[1:nearestNeighborCount],rm.na=T)/nearestNeighborCount
          }
          if(is.na(avgSimOfPathway1to2))
            simBag<-c(simBag,0.0)
          else   simBag<-c(simBag,avgSimOfPathway1to2)
        }
        neighLen<-length(neList)
        if(is.finite(neighLen)&neighLen>0){
          meanOfPathwaySim = mean(simBag,rm.na=T)
          simBagAll<-c(simBagAll,meanOfPathwaySim)
        }
      }
      avgGraphSimBag = sum(simBagAll)/length(V(gr))
      x=c(graph=graphIndex,avgNSim=avgGraphSimBag,gorder=gorder(gr),gsize=gsize(gr),simFile=simFile)
      message("\tGraph id: ",graphIndex," Average Pathway Similarity: ",avgGraphSimBag)
      results= bind_rows(results,x)
    }   
  }
  
  fname=gsub("[^[:alnum:][:space:]]","",simFile)
  
  subDir <- "results"
  ifelse(!dir.exists(file.path(workingDir, subDir)), dir.create(file.path(workingDir, subDir)), FALSE)
  fname=paste0("results",nearestNeighborCount,fname,".rds")
  resFile <-  paste0(workingDir,subDir,"/",fname)
  saveRDS(results,
          file = resFile)
}

for(simFile in simFiles){
  geneSimMap = as.data.table(loadSim(dataDir,simFile));
  colnames(geneSimMap)<-c("gene1","gene2","similarity"); 
  plot2<-ggplot(data = geneSimMap,aes(x=similarity))+  
    geom_histogram (binwidth = 0.05, color = "white") +
    labs (title=paste(substring(simFile,first=79,99)), x="Gene similarity", y="Gene pairs")
  plot2
  fname=gsub("[^[:alnum:][:space:]]","",simFile)
  
  subDir <- "figures"
  ifelse(!dir.exists(file.path(workingDir, subDir)), dir.create(file.path(workingDir, subDir)), FALSE)
  fname=paste0(fname,".png")
  resFile <-  paste0(workingDir,subDir,"/",fname)
  ggsave(filename =resFile,plot=plot2,width=10,unit="cm")
}

for( nearestNeighborCount in c(5,10,20,50)){
  fname=gsub("[^[:alnum:][:space:]]","",simFile)
  results2<-readRDS(file=paste0(workingDir,"results/","results",nearestNeighborCount,fname,".RDS"))
  results2$graph<-as.numeric(results2$graph)
  results2$avgNSim<-round(as.numeric(results2$avgNSim),3)
  
  # a simnple hack for better file names in the legend
  results2$simFile=substr(results2$simFile,33,53)
  
  
  p1<-ggplot(data=results2,aes(x=graph,y=avgNSim,group=simFile,color=simFile))+geom_line(size=1.3)+
    ggtitle(paste0(nearestNeighborCount,"-NN average pathway-pathway similarity in graphs")) +
    scale_y_continuous(name="similarity")+ scale_x_continuous(name="graph id")+
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     panel.border = element_blank(),panel.background = element_blank())+
    theme(legend.position = c(0.2,0.9));p1
  ggsave(filename = paste(nearestNeighborCount,"sim.png",sep=""),plot=p1,width=20,height=15,unit="cm",path=paste0(workingDir,"figures/"))
  
  results2$gorder<-as.numeric(results2$gorder)
  results2$gsize<-as.numeric(results2$gsize)
  ggplot(data=results2,aes(x=graph,y=gorder))+geom_line(size=1) 
  ggplot(data=results2,aes(x=graph,y=gsize))+geom_line(size=1)+scale_y_log10() 
}

# August 2021 code
load("data/BioPlanetNetworks.RData")
simFiles  = c(proDomainSimFile,
              subcellularSimFile,molecularSimFile,bioProcessSimFile
)
for(simFile in simFiles){
  results<-data.frame()
  geneSimMap = as.data.table(loadSim(dataDir,simFile))
  
  # todo: compute pairwise sim of all
  gr=graph_from_edgelist(as.matrix(total.pathway.net[,c("source","target")]))
  # graph is here. 
  # How do we know that this specific crosstalk pathway network is any good?
  # We can quantify edge similarity. If two pathways are connected, 
  # we can compute their similarity by taking the average similarity of their genes.
  # assumption: edges between highly similar pathways are good
 
  nearestNeighborCount<-30
  for(pathway1Name in V(gr)$name){
    pathway1Genes <- pathways[pathways$PATHWAY_NAME==pathway1Name,]$GENE_SYMBOL
    sumOfPathway1Sim<-0.0
    neList <- neighbors(gr,pathway1Name)$name
    df10<-geneSimMap[geneSimMap$geneColumn1%in%pathway1Genes,]
    df20<-geneSimMap[geneSimMap$geneColumn2%in%pathway1Genes,]
    for (pathway2Name in neList){
      pathway2Genes = as.character(pathways[pathways$PATHWAY_NAME==pathway2Name,]$GENE_SYMBOL)
      df1<-df10[df10$geneColumn2%in%pathway2Genes,]
      df2<-df20[df20$geneColumn1%in%pathway2Genes,]
      both = length(intersect(pathway1Genes,pathway2Genes))
      
      
      if(useIdenticalGenes){
        N2 = nearestNeighborCount-both
        if(N2>0){
          sumOfSim = sum(sort(union(df1$simColumn,df2$simColumn),decreasing=T)[1:N2],rm.na=T)
          avgSimOfPathway1to2 = (both+sumOfSim)/(nearestNeighborCount+both)
        }else {
          avgSimOfPathway1to2 = 1.0
        }
      }
      else{
        avgSimOfPathway1to2 = sum(sort(union(df1$simColumn,df2$simColumn),decreasing=T)[1:nearestNeighborCount],rm.na=T)/nearestNeighborCount
      }
      if(is.na(avgSimOfPathway1to2))
        avgSimOfPathway1to2<-0.0
      if(avgSimOfPathway1to2>1)avgSimOfPathway1to2=1.0
      x=c(pw1=pathway1Name,pw2=pathway2Name,sim=avgSimOfPathway1to2,sourceFile=simFile)
      #message(pathway1Name," ",pathway2Name," ",avgSimOfPathway1to2)
      results= bind_rows(results,x)
    }
  }
   saveRDS(results,file=paste0(simFile ,"P2PSim.rds"))
}

#Sept 18 2021
simFiles  = c(proDomainSimFile,
              subcellularSimFile,molecularSimFile,bioProcessSimFile
)
results<-data.frame()
for(simFile in simFiles){
  results<-bind_rows(results,readRDS(file=paste0(simFile ,"P2PSim.rds")))
}
nr<-results
nr$sim<-as.numeric(nr$sim)
nr$sourceFile<-substr(nr$sourceFile,33,(nchar(nr$sourceFile)-4))
nr$sourceFile <- gsub('GO_Biological_Process_2018', 'GoBio', nr$sourceFile)
nr$sourceFile <- gsub('GO_Cellular_Component_2018', 'GoCel', nr$sourceFile)
nr$sourceFile <- gsub('InterPro_Domains_2019', 'InPro', nr$sourceFile)
nr$sourceFile <- gsub('GO_Molecular_Function_2018', 'GoMol', nr$sourceFile)
dataAll<-spread(nr, key = sourceFile, value = sim)
saveRDS(dataAll,file="P2PSim.rds")

load("data/BioPlanetNetworks.RData")
dataAll<-readRDS(file="P2PSim.rds")
total.pathway.net <- read.delim("C:/Code/IntegratedPTMOmics/data/total.pathway.net.txt")
pathway.crosstalk.network<- read.delim("C:/Code/IntegratedPTMOmics/data/pathway.crosstalk.network.txt")
####################################3
# Normalize: put two weights on same scale to make Weight.clust from 0 to 1
tpnn <- total.pathway.net[,c(1:4,6,5)]
tpnn$Weight.clust <- tpnn$Weight.clust/max(tpnn$Weight.clust)
tpnn$Weight.normalized <- tpnn$Weight.clust - tpnn$Weight.bp
tpnn$Combined.Weight <- tpnn$Weight.clust + tpnn$Weight.bp
tpnn.no.bp <- tpnn[which(tpnn$Weight.bp==0), c(1,2,5,6)]
names(tpnn.no.bp)[3:4] <- c("Weight", "interaction")
#
# Do pathway.crosstalk.network
pcnn <- pathway.crosstalk.network
dim(pcnn[which(pcnn$interaction=="cluster evidence"),]) #645709 okay, same as total.pathway.net, tpnn
pcnn.clust <- pcnn[which(pcnn$interaction=="cluster evidence"),]
pcnn.clust$Weight <- pcnn.clust$Weight/max(pcnn.clust$Weight)
pcnn.bp <- pcnn[which(pcnn$interaction=="pathway Jaccard similarity"),]
pcnn <- rbind(pcnn.clust, pcnn.bp)
#

####################################

total.pathway.net$CoWei<-tpnn$Combined.Weight



data2<-merge(dataAll,total.pathway.net, by.x = c("pw1","pw2"), 
      by.y = c("source","target"), all.x = TRUE, all.y = TRUE)

corRes<-cor(data2[,c("GoBio", "CoWei", 
                  "GoCel","InPro","GoMol")])
data <- as.matrix((corRes))

# Default Heatmap
heatmap(data) 
install.packages("gplots")
library("gplots")


# heatmap with the defaults parameters
colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))

my_palette <- colorRampPalette(c("lightblue", "blue", "red"))(n = 299)
heatmap.2(x=data,trace="none",col=my_palette)
