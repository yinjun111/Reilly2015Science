#Code for WGCNA network

#Input File
#expression_matrix_cortex.csv	#RNA-seq gene expression matrix downloaded from Brainspan

#Out Files
#WC2_4.txt	#Filtered gene expression matrix for period 2 to 4
#WC2_4_power.txt	#Power selectionto reach Scale Free Topology
#WC2_4_module.txt	#Module assignment for genes
#WC2_4_all_con.txt	#Connectivity scores of genes
#WC2_4_cmd.txt	#Multidimensional Scaling of genes
#WC2_4_module3_cor0.95_network.txt	#Network of genes by correlation cutoff of 0.95 for module 3
#WC2_4_eigengene.txt	#Module eigengene expression
#WC2_4_cmd_formodules.txt	#Multidimensional Scaling of modules based on eigengene expression
#WC2_4_cor0.3_modulecor.txt	#Network of modules by correlation cutoff of 0.5

#----------------------------------------------
#Code begins

library(WGCNA);
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads =8)

#Step 1, filter the data
#----------------------------------------------

#data table from Brainspan
expMatrix= read.csv("expression_matrix_cortex.csv",check.names=F)

#function to filter the data
filter.good.data<-function(data) {
	#only write the good samples out
	gsgdata = goodSamplesGenes(data, verbose = 3);

	#remove 0 variance genes
	if (!gsgdata$allOK)
	{
	  # Optionally, print the gene and sample names that were removed:
	  if (sum(!gsgdata$goodGenes)>0)
		 printFlush(paste("Removing genes:", paste(names(data)[!gsgdata$goodGenes], collapse = ", ")));
	  if (sum(!gsgdata$goodSamples)>0)
		 printFlush(paste("Removing samples:", paste(rownames(data)[!gsgdata$goodSamples], collapse = ", ")));
	  # Remove the offending genes and samples from the data:
	  gooddata = data[gsgdata$goodSamples, gsgdata$goodGenes]
	}
	gooddata
}	


data<-as.data.frame(t(expMatrix[, 6:84]))
names(data) = expMatrix$ensembl_gene_id;
rownames(data) = names(expMatrix)[6:84];
#get good genes/samples
gooddata<-filter.good.data(data)

#write
write.table(gooddata,file="WC2_4.txt",sep="\t",quote=F)


#Step2, obtain power for each file
#----------------------------------------------
infile<-"WC2_4.txt" #filtered expression file

r2cutoff<-0.80

#outfiles
outfile1<-sub(".txt","_tree.pdf",infile) #tree figure
outfile2<-sub(".txt","_power.pdf",infile) #power figure
outfile3<-sub(".txt","_power.txt",infile) #power file
outfile4<-sub(".txt","_step2_power.rdata",infile) #rdata

data<-read.table(file=infile,header=T,row.names=1,sep="\t",as.is=T,check.names=F)

sampleTree = flashClust(dist(data), method = "average");
pdf(outfile1)
par(cex = 0.8);
par(mar = c(0,4,2,0))
plot(sampleTree, main = paste("Sample clustering to detect outliers for ",infile,sep=""), sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#define power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(data, powerVector = powers, verbose = 5)

pdf(outfile2,width = 9, height = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
r2<- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
plot(sft$fitIndices[,1], r2,
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence for ",infile,sep=""));
text(sft$fitIndices[,1], r2,
    labels=powers,cex=cex1,col="red");
abline(h=r2cutoff,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity for ",infile,sep=""))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#mark power to passing R2 cutoff
write.table(sft$fitIndices,file=outfile3,quote=F,sep="\t")

#save image for future use
save.image(outfile4)


#Step3, build network, and write module assignment and dendrogram
#----------------------------------------------

infile<-"WC2_4.txt" #filtered expression file
npower<-3 #decided from the above analysis

outfile1<-sub(".txt","_module.txt",infile) #genes to color/module
outfile2<-sub(".txt","_dendro.pdf",infile) #colored dendrogram figure
outfile3<-sub(".txt","_step3_network.rdata",infile) #rdata
outfile4<-sub(".txt","_TOM",infile) #rdata for TOM

#parameters to change
no.maxBlockSize=20000

data<-read.table(file=infile,header=T,row.names=1,sep="\t",as.is=T,check.names=F)

#using 8 threads to calculate the modules
bwnet = blockwiseModules(data, maxBlockSize = no.maxBlockSize,power=npower,
                     TOMType = "unsigned", minModuleSize = 30,
                     reassignThreshold = 0, mergeCutHeight = 0.25,
                     numericLabels = T, pamRespectsDendro = FALSE,
                     saveTOMs = T, saveTOMFileBase= outfile4,verbose = 3,nThreads=8)

#plot the softConnectivity #not implemented yet
#scaleFreePlot(TOM)

#save module assignment
# Convert labels to colors for plotting
bwModuleGroups<-bwnet$colors
bwModuleColors = labels2colors(bwnet$colors)

#results
results<-cbind(colnames(data),bwModuleGroups,bwModuleColors,bwnet$blocks)
colnames(results)<-c("Genes","ModuleAssignment","ModuleColor","Blocks")
write.table(results,file=outfile1,sep="\t",quote=F,col.names=NA)

#plot dendrogram
pdf(outfile2)
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]], "Module colors",
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

save.image(outfile3)


#Step 4, calculate correlation and intra module connectivity
#----------------------------------------------

infile<-"WC2_4.txt" #filtered expression file
npower<-3 #decided from the above analysis

matfile<-infile #signal mat
modulefile<-sub(".txt","_module.txt",infile) #module assignment

outfile1<-sub(".txt","_all_con.txt",infile) #genes to connectivity
outfile2<-sub(".txt","_all_con.rdata",infile) #rdata

data<-read.table(file=matfile,header=T,row.names=1,sep="\t",as.is=T,check.names=F)
modules<-read.table(file=modulefile,header=T,row.names=1,sep="\t",as.is=T,check.names=F)

module.cor<-list()
all.con<-c()

for (module in unique(modules[,2])) {
	ADJ=cor(data[,modules[,2]==module],use="p")
	ADJ1=abs(ADJ)^npower
	module.cor[[as.character(module)]]<-ADJ
	#weighted k
	Alldegrees1=intramodularConnectivity(ADJ1, modules[modules[,2]==module,2])
	#unweighted k
	Alldegrees2=intramodularConnectivity(abs(ADJ), modules[modules[,2]==module,2])
	all.con<-rbind(all.con,cbind(Alldegrees1,Alldegrees2))
}

all.con<-all.con[order(rownames(all.con)),]
colnames(all.con)<-c(paste(colnames(all.con)[1:4],"Signed",sep="-"),paste(colnames(all.con)[5:8],"Unsigned",sep="-"))

write.table(file=outfile1,all.con,sep="\t",quote=F)

save.image(outfile2)


#Step 5, export gene network
#----------------------------------------------
load("WC2_4_all_con.rdata")

cor.cutoff<-0.95
module.sel<-3 #selected module to export

module.sel.cor<-module.cor[[module.sel]]
genes.sel<-c()
results<-c()

for(num1 in 1:(ncol(module.sel.cor)-1)) {
	for (num2 in (num1+1):ncol(module.sel.cor)) {
		if(abs(module.sel.cor[num1,num2])>=cor.cutoff) {
			results<-rbind(results,c(colnames(module.sel.cor)[num1],colnames(module.sel.cor)[num2],module.sel.cor[num1,num2]))
			genes.sel<-c(genes.sel,colnames(module.sel.cor)[num1],colnames(module.sel.cor)[num2])
		}
	}
}

cat(paste(length(unique(genes.sel))," out of ",ncol(module.sel.cor) ,"genes in the network",sep=""))

write.table(results,file="WC2_4_module3_cor0.95_network.txt",sep="\t",quote=F,row.names =F,col.names =F)

#Multidimensional Scaling of genes in all the modules

module.sel.cor<-module.cor[[as.character(module.sel)]]
module.cmdscale<-cmdscale(as.dist(1-module.sel.cor),k=3)

colnames(module.cmdscale)<-c("Axis1","Axis2","Axis3")
 
write.table(module.cmdscale,file="WC2_4_module3_cmd.txt",sep="\t",quote=F,col.names=NA)


#Step 6, network of modules by module eigengene 
#----------------------------------------------

load("WC2_4_step3_network.rdata")

cor.cutoff<- 0.3 #cut off for the linkage of modules on network

datME=moduleEigengenes(data,bwnet$colors)$eigengenes

write.table(datME,file="WC2_4_eigengene.txt",sep="\t",quote=F,col.names=NA)

eigen<-read.table("WC2_4_eigengene.txt",header=T,row.names=1,sep="\t")
eigen.sel<-eigen[,2:ncol(eigen)] #remove grey module

#correlation
eigen.cor<-cor(eigen.sel,method="pearson")

#MultiDimensional Scaling
module.eigen.cmdscale<-cmdscale(as.dist(1-eigen.cor),k=3)

write.table(file="WC2_4_cmd_formodules.txt",module.eigen.cmdscale*500,sep="\t",quote=F,col.names=NA)

#export links

genes.sel<-c()
results<-c()

for(num1 in 1:(ncol(eigen.cor)-1)) {
	for (num2 in (num1+1):ncol(eigen.cor)) {
		if(abs(eigen.cor[num1,num2])>=cor.cutoff) {
			results<-rbind(results,c(colnames(eigen.cor)[num1],colnames(eigen.cor)[num2],eigen.cor[num1,num2]))
			genes.sel<-c(genes.sel,colnames(eigen.cor)[num1],colnames(eigen.cor)[num2])
		}
	}
}

#Write file for CytoScape
interaction.type<-results[,3]
interaction.type[interaction.type>0]=1
interaction.type[interaction.type<0]=-1

results<-cbind(results,interaction.type)

cat(paste(length(unique(genes.sel))," out of ",ncol(eigen.cor) ," modules in the network\n",sep=""))
#file that can be imported to Cytoscape
write.table(results,file="WC2_4_cor0.5_modulecor.txt",sep="\t",quote=F,row.names =F,col.names =F)

