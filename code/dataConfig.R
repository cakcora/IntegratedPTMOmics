

dataDir = "C:/Code/IntegratedPTMOmics/data/"

#Pathway file comes from BioPlanet
pathwayFile = paste0(dataDir,"bioplanet_pathway.csv")
bioPlanetsRdata<-paste0(dataDir,"BioPlanetNetworks.RData")

#Clusters come from Mark's experiments
clusterFile = paste(dataDir,"sites_by_cluster.txt",sep="")

#Files to be used in gene similarity computations
#We needed to do the following modifications on the original files:
#1 - The line endings contained many repeating \t. We removed them. 
#2- Change line endings with "\t\t\r\n" to "\r\n"
#3- Change pathway name and remove commas in the name
#4- Change pathway name and gene list seperator from "\t\t" to ","
subcellularSimFile =paste0(dataDir,"GO_Cellular_Component_2018.txt")
molecularSimFile = paste0(dataDir,"GO_Molecular_Function_2018.txt")
proDomainSimFile = paste0(dataDir,"InterPro_Domains_2019.txt")
bioProcessSimFile =paste0(dataDir, "GO_Biological_Process_2018.txt")