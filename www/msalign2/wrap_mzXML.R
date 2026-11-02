###
### This script assumes that XML package is installed
library(XML)
###
###
all.args=commandArgs()
#
###***************************************************************************************
### The PLF file is in retention time (seconds): column 1 is the dataset to align,
### column 2 is the master dataset.
###
plfFile = paste(all.args[6], ".alignment.plf", sep="")
DataFile1 = all.args[6]
###
### Read the PLF table
plf_tab <- read.table(plfFile, header = FALSE)
### Keep only breakpoints which are strictly increasing in both columns, this keeps the order of the scans
keep <- rep(FALSE, nrow(plf_tab))
last <- c(-Inf, -Inf)
for (j in 1:nrow(plf_tab)) {
	if ((plf_tab[j,1] > last[1]) && (plf_tab[j,2] > last[2])) {
		keep[j] <- TRUE
		last <- c(plf_tab[j,1], plf_tab[j,2])
	}
}
plf_tab <- plf_tab[keep, ]
###
### Parse the mzXML file
dataset1 <- xmlTreeParse(DataFile1, addAttributeNamespaces = TRUE, encoding = "ISO-8859-1", useInternalNodes = TRUE)
root1 <- xmlRoot(dataset1)
### Find all scans (MS1 and MS2, also when nested) in the order of the file
scans <- getNodeSet(dataset1, "//*[local-name()='scan']")
retTime <- sapply(scans, function(scan) as.numeric(sub("^PT(.*)S$", "\\1", xmlGetAttr(scan, "retentionTime"))))
###
### Map the retention time of every scan through the PLF
new_retTime <- approx(plf_tab[,1], plf_tab[,2], xout = retTime, rule = 2)$y
### Make sure the order of the scans did not change
if (!is.unsorted(retTime, strictly = TRUE) && is.unsorted(new_retTime, strictly = TRUE)) {
	stop("aligned retention times are not strictly increasing, the order of the scans would change. No aligned mzXML file written.")
}
for (s in seq_along(scans)) {
	addAttributes(scans[[s]], "retentionTime"=paste("PT", new_retTime[s], "S", sep=""))
}
### Update the start and end time of the run
msRun <- getNodeSet(dataset1, "//*[local-name()='msRun']")[[1]]
invisible(addAttributes(msRun, "startTime"=paste("PT", min(new_retTime), "S", sep=""), "endTime"=paste("PT", max(new_retTime), "S", sep="")))
saveXML(root1, file = paste(all.args[6], "_ALIGNED.mzxml", sep=""), prefix = "<?xml version=\"1.0\"?>\n", encoding = "ISO-8859-1")
########################################
