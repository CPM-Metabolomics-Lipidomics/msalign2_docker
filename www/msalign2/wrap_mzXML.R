###
### This script assumes that XML package is installed
library(XML)
###
###
all.args=commandArgs()
#
###***************************************************************************************
### There are 3 variables you have to set for your own dataset: PLF file, two datasets
### 
plfFile = paste(all.args[6], ".alignment.plf", sep="") 
DataFile1 = all.args[6]
DataFile2 = all.args[7]
###
### Read the PLF table
plf_tab <- read.table(plfFile, header = FALSE)
new_plftab <- plf_tab
###
### Parse the first mzXML file
dataset1 <- xmlTreeParse(DataFile1, addAttributeNamespaces = TRUE, encoding = "ISO-8859-1", useInternalNodes = TRUE)
### Parse the second mzXML file
dataset2 <- xmlTreeParse(DataFile2, addAttributeNamespaces = TRUE, encoding = "ISO-8859-1", useInternalNodes = TRUE)
### Find the root of XML
root1 <- xmlRoot(dataset1)
### Find the first child
node1 <- xmlChildren(root1)
node2 <- xmlChildren(node1$msRun)
### Find the last scan
last_scan1 <- length(node2)
### Find relations between scan numbers in PLF and retention time in mzXML
for (i in 1:last_scan1) {
	if (names(node2[i])=="scan") {
	scanNum <- as.numeric(xmlAttrs(node2[[i]])["num"][[1]])
		for (j in 1:(nrow(plf_tab)-1)) {
		if (identical(floor(plf_tab[j,1]),scanNum)) {
			retTime <- xmlAttrs(node2[[i]])["retentionTime"][[1]]
			new_plftab[j,1] <- as.numeric(strsplit(strsplit(retTime,"T")[[1]][2],"S"))		
		}
	}
	}
	else { 
		start <- i
	}
}
### Add the last retention time to PLF table
retTime <- xmlAttrs(node2[[last_scan1]])["retentionTime"][[1]]
new_plftab[(nrow(plf_tab)),1] <- as.numeric(strsplit(strsplit(retTime,"T")[[1]][2],"S"))		
###*******************************************************************************************
### Follow the same procedure for the second dataset
root2 <- xmlRoot(dataset2)
### Find the first child
node1_2 <- xmlChildren(root2)
node2_2 <- xmlChildren(node1_2$msRun)
### Find the last scan
last_scan2 <- length(node2_2)
### Find relations between scan numbers in PLF and retention time in mzXML
for (k in 1:last_scan2) {
	if (names(node2_2[k])=="scan") {
	scanNum <- as.numeric(xmlAttrs(node2_2[[k]])["num"][[1]])
		for (l in 1:(nrow(plf_tab)-1)) {
		if (identical(floor(plf_tab[l,2]),scanNum)) {
			retTime <- xmlAttrs(node2_2[[k]])["retentionTime"][[1]]
			new_plftab[l,2] <- as.numeric(strsplit(strsplit(retTime,"T")[[1]][2],"S"))		
		}
	}
	}
}
### Add the last retention time to PLF table
retTime <- xmlAttrs(node2_2[[last_scan2]])["retentionTime"][[1]]
new_plftab[nrow(new_plftab),2] <- as.numeric(strsplit(strsplit(retTime,"T")[[1]][2],"S"))		

###################
for (i in (start+1):(last_scan1-1)){
	retTime <- xmlAttrs(node2[[i]])["retentionTime"][[1]]
	retTime <- as.numeric(strsplit(strsplit(retTime,"T")[[1]][2],"S"))
		for (k in 1:(nrow(new_plftab)-1)) {
			if ((retTime >= new_plftab[k,1]) && (retTime < new_plftab[k+1,1]))  
			{
			new_retTime <- new_plftab[k,2]+((retTime-new_plftab[k,1])/(new_plftab[k+1,1]-new_plftab[k,1]))*(new_plftab[k+1,2]-new_plftab[k,2])
                    	new_retTime <- paste("PT", new_retTime,"S", sep="")
			
			node2[[i]] <- addAttributes(node2[[i]],"retentionTime"=new_retTime)	
			}
		}		
}
last_retTime <- paste("PT", new_plftab[nrow(new_plftab),2],"S", sep="") 
node2[[last_scan1]] <- addAttributes(node2[[last_scan1]],"retentionTime"=last_retTime)	
xmlChildren(node1$msRun) <- node2
xmlChildren(root1) <- node1
saveXML(root1, file = paste(all.args[6], "_ALIGNED.mzxml", sep=""), prefix = "<?xml version=\"1.0\"?>\n", encoding = "ISO-8859-1")
########################################
