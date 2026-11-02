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
### Map retention times (seconds) through the PLF
alignTimes <- function(retTime) {
	new_retTime <- approx(plf_tab[,1], plf_tab[,2], xout = retTime, rule = 2)$y
	### Make sure the order of the scans did not change
	if (!is.unsorted(retTime, strictly = TRUE) && is.unsorted(new_retTime, strictly = TRUE)) {
		stop("aligned retention times are not strictly increasing, the order of the scans would change. No aligned file written.")
	}
	new_retTime
}
###
### Is it an mzML or mzXML file?
fileStart <- readChar(DataFile1, 4000, useBytes = TRUE)
isMzML <- grepl("<mzML", fileStart, fixed = TRUE)
###
if (!isMzML) {
	###
	### Parse the mzXML file
	dataset1 <- xmlTreeParse(DataFile1, addAttributeNamespaces = TRUE, encoding = "ISO-8859-1", useInternalNodes = TRUE)
	root1 <- xmlRoot(dataset1)
	### Find all scans (MS1 and MS2, also when nested) in the order of the file
	scans <- getNodeSet(dataset1, "//*[local-name()='scan']")
	retTime <- sapply(scans, function(scan) as.numeric(sub("^PT(.*)S$", "\\1", xmlGetAttr(scan, "retentionTime"))))
	###
	### Map the retention time of every scan through the PLF
	new_retTime <- alignTimes(retTime)
	for (s in seq_along(scans)) {
		addAttributes(scans[[s]], "retentionTime"=paste("PT", new_retTime[s], "S", sep=""))
	}
	### Update the start and end time of the run
	msRun <- getNodeSet(dataset1, "//*[local-name()='msRun']")[[1]]
	invisible(addAttributes(msRun, "startTime"=paste("PT", min(new_retTime), "S", sep=""), "endTime"=paste("PT", max(new_retTime), "S", sep="")))
	saveXML(root1, file = paste(all.args[6], "_ALIGNED.mzxml", sep=""), prefix = "<?xml version=\"1.0\"?>\n", encoding = "ISO-8859-1")
} else {
	###
	### Parse the mzML file
	dataset1 <- xmlTreeParse(DataFile1, addAttributeNamespaces = TRUE, encoding = "UTF-8", useInternalNodes = TRUE)
	### The retention time is the scan start time cvParam of each spectrum, in the order of the file
	### (a single step XPath, nested steps are very slow on large files with this libxml2)
	rtParams <- getNodeSet(dataset1, "//*[local-name()='cvParam' and @accession='MS:1000016']")
	### Factor to convert the unit of the scan start time to seconds
	unitFactor <- sapply(rtParams, function(cv) {
		unitAcc <- xmlGetAttr(cv, "unitAccession", default = "")
		unitName <- xmlGetAttr(cv, "unitName", default = "")
		if (unitAcc == "UO:0000031" || unitName == "minute") {
			60
		} else if (unitAcc == "UO:0000028" || unitName == "millisecond") {
			0.001
		} else {
			1
		}
	})
	retTime <- sapply(rtParams, function(cv) as.numeric(xmlGetAttr(cv, "value"))) * unitFactor
	###
	### Map the retention time of every spectrum through the PLF, write it back in the original unit
	new_retTime <- alignTimes(retTime)
	for (s in seq_along(rtParams)) {
		addAttributes(rtParams[[s]], "value"=as.character(new_retTime[s] / unitFactor[s]))
	}
	###
	### The time axis of the chromatograms (TIC, BPC) is not aligned, remove them.
	### Also remove the old index, it is rebuilt below.
	invisible(removeNodes(getNodeSet(dataset1, "//*[local-name()='chromatogramList' or local-name()='indexList' or local-name()='indexListOffset' or local-name()='fileChecksum']")))
	###
	### Write the mzML element wrapped in a new indexedmzML element
	mzmlNode <- getNodeSet(dataset1, "//*[local-name()='mzML']")[[1]]
	mzmlText <- saveXML(mzmlNode)
	### Make sure the namespaces are declared on the mzML element itself
	mzmlTag <- regmatches(mzmlText, regexpr("<mzML[^>]*>", mzmlText))
	if (!grepl("xmlns=", mzmlTag, fixed = TRUE)) {
		mzmlText <- sub("<mzML", "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\"", mzmlText, fixed = TRUE)
	}
	if (grepl("xsi:", mzmlTag, fixed = TRUE) && !grepl("xmlns:xsi=", mzmlTag, fixed = TRUE)) {
		mzmlText <- sub("<mzML", "<mzML xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", mzmlText, fixed = TRUE)
	}
	body <- paste0("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",
		"<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ",
		"xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd\">\n",
		mzmlText, "\n")
	###
	### Build the spectrum index: byte offset of every spectrum element
	spectrumIds <- sapply(getNodeSet(dataset1, "//*[local-name()='spectrum']"), function(spectrum) xmlGetAttr(spectrum, "id"))
	offsets <- gregexpr("<spectrum ", body, fixed = TRUE, useBytes = TRUE)[[1]] - 1
	if (length(spectrumIds) == 0) {
		offsets <- integer(0)
	}
	if (length(offsets) != length(spectrumIds)) {
		stop("could not determine the offsets of all spectra. No aligned mzML file written.")
	}
	spectrumIds <- gsub("&", "&amp;", spectrumIds, fixed = TRUE)
	spectrumIds <- gsub("<", "&lt;", spectrumIds, fixed = TRUE)
	spectrumIds <- gsub(">", "&gt;", spectrumIds, fixed = TRUE)
	spectrumIds <- gsub("\"", "&quot;", spectrumIds, fixed = TRUE)
	indexText <- paste0("<indexList count=\"1\">\n<index name=\"spectrum\">\n",
		paste0("<offset idRef=\"", spectrumIds, "\">", format(offsets, scientific = FALSE, trim = TRUE), "</offset>\n", collapse = ""),
		"</index>\n</indexList>\n",
		"<indexListOffset>", format(nchar(body, type = "bytes"), scientific = FALSE), "</indexListOffset>\n",
		"<fileChecksum>")
	###
	### The checksum is the SHA-1 of the file up to and including <fileChecksum>
	outFile <- paste(sub("\\.mzML$", "", DataFile1, ignore.case = TRUE), "_ALIGNED.mzML", sep="")
	con <- file(outFile, "wb")
	writeBin(charToRaw(body), con)
	writeBin(charToRaw(indexText), con)
	close(con)
	checksum <- strsplit(system2("sha1sum", shQuote(outFile), stdout = TRUE), " ")[[1]][1]
	cat(checksum, "</fileChecksum>\n</indexedmzML>\n", file = outFile, append = TRUE, sep = "")
}
########################################
