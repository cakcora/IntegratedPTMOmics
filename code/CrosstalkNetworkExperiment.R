library(data.table)
library(plyr)
library(dplyr)
library(igraph)
library(lsa)
library(ggplot2)

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

for( nearestNeighborCount in c(5,10,20,50,100,200)){
  fname=gsub("[^[:alnum:][:space:]]","",simFile)
  results2<-readRDS(file=paste0(workingDir,"results/","results",nearestNeighborCount,fname,".RDS"))
  results2$graph<-as.numeric(results2$graph)
  results2$simFile<-fname
  results2$avgSim<-round(as.numeric(results2$avgNSim),3)
  p1<-ggplot(data=results2,aes(x=graph,y=avgSim,group=simFile,color=simFile))+geom_line(size=1.3)+
    ggtitle(paste0(nearestNeighborCount,"-NN average graph similarity")) +theme(legend.position = c(0.2,0.9));p1
  ggsave(filename = paste(nearestNeighborCount,"sim.png",sep=""),plot=p1,width=20,height=15,unit="cm",path=paste0(workingDir,"results/"))
  
  results2$gorder<-as.numeric(results2$gorder)
  results2$gsize<-as.numeric(results2$gsize)
  ggplot(data=results2,aes(x=graph,y=gorder))+geom_line(size=1) 
  ggplot(data=results2,aes(x=graph,y=gsize))+geom_line(size=1)+scale_y_log10() 
}
