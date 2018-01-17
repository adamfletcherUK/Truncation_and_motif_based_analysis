################################################################################################################
# Script 3: Merging the truncating mutation data from the TCGA and CCLE truncating scripts (scripts 1 and 2)   #
#                                                                                                              #
# Automatic merge function to combine the TCGA and CCLE trucation mutation data outputted from scripts 1 and 2.#
################################################################################################################ 

# Intro text.
cat("A. Hudson, N. Stephenson, C. Li, E. Trotter, A. Fletcher, G. Katrona, P. Bieniasz-Krzywiec, M. Howell, C. Writh, C. Miller, J. Brognard. \n\n",
    "Functional Screening of Cancer Datasets Identifies Mutational Hotspots and Novel Targets in the Tumour-Suppressing Kinome. \n\n", 
    "Script requires stringr and plyr packages, additionally set the working directory to the folder containing the SOURCE.csv files \n\n", 
    sep="")

# Load script dependencies
install.packages("stringr")
install.packages("plyr")
library("stringr", "plyr")


########################
# Merge data from TCGA and CCLE scripts
########################

# Merge the two SHORTEST truncating scores (using dplyr).
ccletruncSHORTEST <- read.csv (file = "TEMP_lengthcorrectedfreq_ccleALLSHORTEST.csv")
tcgatruncSHORTEST <- read.csv (file = "TEMP_lengthcorrectedfreq_tcgaALLSHORTEST.csv")
combinedTRUNCSHORTEST <- rbind(ccletruncSHORTEST, tcgatruncSHORTEST)
combinedTRUNCscoresSHORTEST <- ddply(combinedTRUNCSHORTEST,"Gene",numcolwise(sum))
write.csv (combinedTRUNCscoresSHORTEST, file ='outputs/mergedTRUNCscoresSHORTEST.csv')


# Merge the two LONGEST truncating scores (using dplyr).
ccletruncLONGEST <- read.csv (file = 'TEMP_lengthcorrectedfreq_ccleALLLONGEST.csv')
tcgatruncLONGEST <- read.csv (file = 'TEMP_lengthcorrectedfreq_tcgaALLLONGEST.csv')
combinedTRUNCLONGEST <- rbind(ccletruncLONGEST, tcgatruncLONGEST)
combinedTRUNCscoresLONGEST <- ddply(combinedTRUNCLONGEST,"Gene",numcolwise(sum))
write.csv (combinedTRUNCscoresLONGEST, file ='outputs/mergedTRUNCscoresLONGEST.csv')


# Calculate the mean of LONGEST and SHORTEST scores.
combinedTRUNCscoresMEAN <- combinedTRUNCscoresSHORTEST
combinedTRUNCscoresMEAN$longestSCORE <- as.numeric(combinedTRUNCscoresLONGEST$score)
colnames(combinedTRUNCscoresMEAN)[7] <- "shortestSCORE"
combinedTRUNCscoresMEAN$shortestSCORE <- as.numeric (combinedTRUNCscoresMEAN$shortestSCORE)
combinedTRUNCscoresMEAN$meanSCORE = (combinedTRUNCscoresMEAN$shortestSCORE + combinedTRUNCscoresMEAN$longestSCORE) / 2
write.csv (combinedTRUNCscoresMEAN, file ='outputs/mergedTRUNCscoresMEAN.csv')
SOURCE_Top30 <- combinedTRUNCscoresMEAN[with(combinedTRUNCscoresMEAN, order(-meanSCORE)),]
SOURCE_Top30 <- head(SOURCE_Top30, n=30L)
SOURCE_Top30 <- as.data.frame(SOURCE_Top30$Gene)
colnames(SOURCE_Top30)[1] <- "Gene"
write.csv (SOURCE_Top30, file="SOURCE_Top30.csv")


# Combine the data from the validate_CCLE and validate_TCGA files.
validateTCGA <- read.csv (file = "Validate_TCGA.csv")
validateCCLE <- read.csv (file = "Validate_CCLE.csv")
validateCCLE$study <- NULL
top30check <- rbind (validateTCGA, validateCCLE) 


# Check the top30 to make sure the all transcript analysis doesnt give any false positives.
top30 <- read.csv ('SOURCE_Top30.csv')
top30check <- merge (top30, top30check, by = 'Gene', nomatch = 0)
write.csv (top30check, file ='top30check.csv')
# NOTE: The above showed that the top 30 hits calculated were valiid.


########################
# Extract the catalytic fragment of the top 30 kinases from the truncation screen. 
########################

# NOTE: There are mutliple different transcripts for each kinase, the longest match was used.
# Rank by decreasing fragment size and remove duplicates, leaving only the longest transcript
FRAGEXTRACT <- genefreqLONGEST
FRAGEXTRACT$DUP= !duplicated(FRAGEXTRACT$Gene)
FRAGEXTRACT <- subset(FRAGEXTRACT, DUP == 'TRUE')
FRAGEXTRACTtop30 <- merge (top30, FRAGEXTRACT, by = 'Gene', nomatch = 0)
#AMHR2 was changed manually as it is incorrect
write.csv (FRAGEXTRACTtop30, file ='TEMP_check.csv')



