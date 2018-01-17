#################################################################################################################
# Script 2: Identifying truncation mutations within the CCLE database that occur before the end of the kinase   #
# domain.                                                                                                       #
#                                                                                                               #
# Automatic search function to detect CCLE trucation mutations occuring before the end of kinase domains.       #
#################################################################################################################

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
# Load Genbank data and determine motif positions
########################

# Loads and merges definitive gene names with GENBANK
genbank <- read.csv ( file = 'SOURCE_DEFINITIVE_GENBANK_TRIM.csv') #PROTEIN SEQ FOR ALL TRANSCRIPTS OF ALL GENES
names   <- read.csv (file = 'SOURCE_411kinases_genename.csv') # Loads genenames for the 411 kinases
genbankkinases = merge(genbank, names, by = 'Gene', nomatch = 0)


# Extract DFG motif information for the kinase sequence.
DFGmotifs  <- as.vector(genbankkinases$DFGmotif)


# Count the number of occurences of DFG motifs in each seq (using StringR package).
numberofDFG  <- str_count(genbankkinases$Protein_Seq,  DFGmotifs)
genbankkinases$numberofDFG  <- numberofDFG


# Define the correct motif to use.
# NOTE:  SOURCE_motif2use is a file with the correct motif to use when there are multiple. Correct motif 
motifToUse     <- read.csv ('SOURCE_motif2use.csv')
genbankkinases<- merge    (genbankkinases, motifToUse, by = c('f','Gene', 'Length', 'Protein_Seq', 'VAIKmotif', 'HRDmotif', 'DFGmotif'))


# Locate all motifs
DFGmotifs  <- as.vector(genbankkinases$DFGmotif)


########################
# Locate all DFG motifs in the peptide sequence
########################

# Locate the first DFG motif.
firstDFG <- str_locate(genbankkinases$Protein_Seq, DFGmotifs)
firstDFG <- as.data.frame(firstDFG)
genbankkinases$firstDFG <- firstDFG$start
firstDFG <- as.numeric(firstDFG$start)
endofstring = as.numeric(firstDFG +1000000)
beyondfirstDFG <- substr(genbankkinases$Protein_Seq, firstDFG + 1, endofstring) #prints the rest of the sequence after the motif


# Locate the second DFG motif.
secondDFG <- str_locate(beyondfirstDFG, DFGmotifs)
secondDFG <- as.data.frame(secondDFG)
genbankkinases$secondDFG <- secondDFG$start
secondDFG = as.numeric(genbankkinases$secondDFG + firstDFG)
genbankkinases$secondDFG <- secondDFG
beyondsecondDFG <- substr(genbankkinases$Protein_Seq, secondDFG + 1, endofstring) #prints the rest of the sequence after the motif


# Locate the third DFG motif.
thirdDFG <- str_locate(beyondsecondDFG, DFGmotifs)
thirdDFG <- as.data.frame(thirdDFG)
genbankkinases$thirdDFG <- thirdDFG$start
thirdDFG = as.numeric(genbankkinases$thirdDFG + secondDFG)
genbankkinases$thirdDFG <- thirdDFG


# Determine the correct DFG to use.
genbankkinases$actual_DFG_D <- ifelse (genbankkinases$DFGtouse == '1', genbankkinases$firstDFG, ifelse (genbankkinases$DFGtouse == '2', genbankkinases$secondDFG , genbankkinases$thirdDFG))
genbankkinases$actual_DFG_F <- genbankkinases$actual_DFG_D + 1
genbankkinases$actual_DFG_G <- genbankkinases$actual_DFG_D + 2


########################
# Locates APE and GxGxxG motif information
########################

# Determine an APE fragment based upon the DFG motif. 
# NOTE: Fragment houses the APE and GxGxxG motifs.
APEfragstart <- genbankkinases$actual_DFG_G + 10
genbankkinases$APEfragstart <- as.numeric(APEfragstart)
APEfragend <- genbankkinases$actual_DFG_G + 70
genbankkinases$APEfragend <- as.numeric(APEfragend)
APEfrag <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APEfragstart, genbankkinases$APEfragend )
genbankkinases$APEfrag <- APEfrag


# Locate the APE motif.
numAPE <- str_count(APEfrag,  'APE')
genbankkinases$numAPE <- as.numeric(numAPE)
APELoc <- str_locate(APEfrag, 'APE')
APELoc <- as.data.frame(APELoc)
genbankkinases$APELoc <- as.numeric(APELoc$start)
genbankkinases$locAPEprotein <- genbankkinases$APEfragstart + genbankkinases$APELoc - 1


# Locate the PE motif.
numPE <- str_count(APEfrag,  'PE')
genbankkinases$numPE <- as.numeric(numPE)
PELoc <- str_locate(APEfrag, 'PE')
PELoc <- as.data.frame(PELoc)
genbankkinases$PELoc <- as.numeric(PELoc$start)
genbankkinases$locPEprotein <- genbankkinases$APEfragstart + genbankkinases$PELoc - 2 # minus 2 as PE is 1 further back 


# Find GTxxxxxxE.
numGTx6E <- str_count(APEfrag,  'GT[A-Z]{6}E')
genbankkinases$numGTx6E  <- as.numeric(numGTx6E)
GTx6ELoc <- str_locate(APEfrag, 'GT[A-Z]{6}E')
GTx6ELoc <- as.data.frame(GTx6ELoc)
genbankkinases$GTx6ELoc <- as.numeric(GTx6ELoc$start)
genbankkinases$locGTx6Eprotein <- genbankkinases$APEfragstart + genbankkinases$GTx6ELoc + 5


# Find GTxxxxxNE.
numGTx5NE <- str_count(APEfrag,  'GT[A-Z]{5}NE')
genbankkinases$numGTx5NE  <- as.numeric(numGTx5NE)
GTx5NELoc <- str_locate(APEfrag, 'GT[A-Z]{5}NE')
GTx5NELoc <- as.data.frame(GTx5NELoc)
genbankkinases$GTx5NELoc <- as.numeric(GTx5NELoc$start)
genbankkinases$locGTx5NEprotein <- genbankkinases$APEfragstart + genbankkinases$GTx5NELoc + 5


# Find GTxxxxxxD.
numGTx6D <- str_count(APEfrag,  'GT[A-Z]{6}D')
genbankkinases$numGTx6D  <- as.numeric(numGTx6D)
GTx6DLoc <- str_locate(APEfrag, 'GT[A-Z]{6}D')
GTx6DLoc <- as.data.frame(GTx6DLoc)
genbankkinases$GTx6DLoc <- as.numeric(GTx6DLoc$start)
genbankkinases$locGTx6Dprotein <- genbankkinases$APEfragstart + genbankkinases$GTx6DLoc + 5


# Find AxE.
numAxE <- str_count(APEfrag,  'A[A-Z]E')
genbankkinases$numAxE <- as.numeric(numAxE)
AxELoc <- str_locate(APEfrag, 'A[A-Z]E')
AxELoc <- as.data.frame(AxELoc)
genbankkinases$AxELoc <- as.numeric(AxELoc$start)
genbankkinases$locAxEprotein <- genbankkinases$APEfragstart + genbankkinases$AxELoc - 1


# Find APD.
numAPD <- str_count(APEfrag,  'APD')
genbankkinases$numAPD  <- as.numeric(numAPD)
APDLoc <- str_locate(APEfrag, 'APD')
APDLoc <- as.data.frame(APDLoc)
genbankkinases$APDLoc <- as.numeric(APDLoc$start)
genbankkinases$locAPDprotein <- genbankkinases$APEfragstart + genbankkinases$APDLoc - 1


# Find PPD.
numPPD <- str_count(APEfrag,  'PPD')
genbankkinases$numPPD  <- as.numeric(numPPD)
PPDLoc <- str_locate(APEfrag, 'PPD')
PPDLoc <- as.data.frame(PPDLoc)
genbankkinases$PPDLoc <- as.numeric(PPDLoc$start)
genbankkinases$locPPDprotein <- genbankkinases$APEfragstart + genbankkinases$PPDLoc - 1


# Find GTxxY.
numGTxxY <- str_count(APEfrag,  'GT[A-Z]{2}Y')
genbankkinases$numGTxxY  <- as.numeric(numGTxxY)
GTxxYLoc <- str_locate(APEfrag, 'GT[A-Z]{2}Y')
GTxxYLoc <- as.data.frame(GTxxYLoc)
genbankkinases$GTxxYLoc <- as.numeric(GTxxYLoc$start)
genbankkinases$locGTxxYprotein <- genbankkinases$APEfragstart + genbankkinases$GTxxYLoc + 5


# Find PIR.
numPIR <- str_count(APEfrag,  'PIR')
genbankkinases$numPIR  <- as.numeric(numPIR)
PIRLoc <- str_locate(APEfrag, 'PIR')
PIRLoc <- as.data.frame(PIRLoc)
genbankkinases$PIRLoc <- as.numeric(PIRLoc$start)
genbankkinases$locPIRprotein <- genbankkinases$APEfragstart + genbankkinases$PIRLoc + 4


# Find Gx7E.
numGx7E <- str_count(APEfrag,  'G[A-Z]{7}E')
genbankkinases$numGx7E  <- as.numeric(numGx7E)
Gx7ELoc <- str_locate(APEfrag, 'G[A-Z]{7}E')
Gx7ELoc <- as.data.frame(Gx7ELoc)
genbankkinases$Gx7ELoc <- as.numeric(Gx7ELoc$start)
genbankkinases$locGx7Eprotein <- genbankkinases$APEfragstart + genbankkinases$Gx7ELoc + 5


# Find YxAP (captures MAPKAPK5).
numYxAP <- str_count(APEfrag,  'Y[A-Z]{1}AP')
genbankkinases$numYxAP  <- as.numeric(numYxAP)
YxAPLoc <- str_locate(APEfrag, 'Y[A-Z]{1}AP')
YxAPLoc <- as.data.frame(YxAPLoc)
genbankkinases$YxAPLoc <- as.numeric(YxAPLoc$start)
genbankkinases$YxAPLocprotein <- genbankkinases$APEfragstart + genbankkinases$YxAPLoc + 1


# Find WYxxPR (captures MAPK4 and MAPK6).
numWYxxPR <- str_count(APEfrag,  'WY[A-Z]{2}PR')
genbankkinases$numWYxxPR  <- as.numeric(numWYxxPR)
WYxxPRLoc <- str_locate(APEfrag, 'WY[A-Z]{2}PR')
WYxxPRLoc <- as.data.frame(WYxxPRLoc)
genbankkinases$WYxxPRLoc <- as.numeric(WYxxPRLoc$start)
genbankkinases$WYxxPRLocprotein <- genbankkinases$APEfragstart + genbankkinases$WYxxPRLoc + 2


########################
# Determines correct APE and location information
########################

# Choose correct APE.
genbankkinases$loc_actual_APE_A <- ifelse (genbankkinases$numAPE <1, genbankkinases$locGTx5NEprotein, genbankkinases$locAPEprotein)
genbankkinases$loc_actual_APE_B <- genbankkinases$loc_actual_APE_A
genbankkinases$loc_actual_APE_C   <- ifelse (genbankkinases$loc_actual_APE_A %in% NA, genbankkinases$locGTx6Eprotein, genbankkinases$loc_actual_APE_B)
genbankkinases$loc_actual_APE_D <- genbankkinases$loc_actual_APE_C 
genbankkinases$loc_actual_APE_E   <- ifelse (genbankkinases$loc_actual_APE_C %in% NA, genbankkinases$locPEprotein, genbankkinases$loc_actual_APE_D)
genbankkinases$loc_actual_APE_F <- genbankkinases$loc_actual_APE_E
genbankkinases$loc_actual_APE_G   <- ifelse (genbankkinases$loc_actual_APE_E %in% NA, genbankkinases$locGTx6Dprotein, genbankkinases$loc_actual_APE_F)
genbankkinases$loc_actual_APE_H <- genbankkinases$loc_actual_APE_G
genbankkinases$loc_actual_APE_I   <- ifelse (genbankkinases$loc_actual_APE_G %in% NA, genbankkinases$locGTxxYprotein, genbankkinases$loc_actual_APE_H)
genbankkinases$loc_actual_APE_J <- genbankkinases$loc_actual_APE_I
genbankkinases$loc_actual_APE_K   <- ifelse (genbankkinases$loc_actual_APE_I %in% NA, genbankkinases$locAPDprotein, genbankkinases$loc_actual_APE_J)
genbankkinases$loc_actual_APE_L  <- genbankkinases$loc_actual_APE_K
genbankkinases$loc_actual_APE_M   <- ifelse (genbankkinases$loc_actual_APE_K %in% NA, genbankkinases$locPPDprotein, genbankkinases$loc_actual_APE_L)
genbankkinases$loc_actual_APE_N  <- genbankkinases$loc_actual_APE_M
genbankkinases$loc_actual_APE_O   <- ifelse (genbankkinases$loc_actual_APE_M %in% NA, genbankkinases$WYxxPRLocprotein, genbankkinases$loc_actual_APE_N)
genbankkinases$loc_actual_APE_P <- genbankkinases$loc_actual_APE_O
genbankkinases$loc_actual_APE_Q   <- ifelse (genbankkinases$loc_actual_APE_O %in% NA, genbankkinases$locAxEprotein, genbankkinases$loc_actual_APE_P)
genbankkinases$loc_actual_APE_R <- genbankkinases$loc_actual_APE_Q
genbankkinases$loc_actual_APE_S   <- ifelse (genbankkinases$loc_actual_APE_Q %in% NA, genbankkinases$locPIRprotein, genbankkinases$loc_actual_APE_R)
genbankkinases$loc_actual_APE_T <- genbankkinases$loc_actual_APE_S
genbankkinases$loc_actual_APE_U   <- ifelse (genbankkinases$loc_actual_APE_S %in% NA, genbankkinases$locGx7Eprotein, genbankkinases$loc_actual_APE_T)
genbankkinases$loc_actual_APE_V <- genbankkinases$loc_actual_APE_U
genbankkinases$loc_actual_APE   <- ifelse (genbankkinases$loc_actual_APE_U %in% NA, genbankkinases$YxAPLocprotein, genbankkinases$loc_actual_APE_V)

genbankkinases$APE_E1loc   <- genbankkinases$loc_actual_APE + 2


# Define the end of the kinase domain (APE + defined end).
genbankkinases$endCAT <- as.numeric(genbankkinases$APE_E1loc)


########################
# Identifying Truncating Mutations
########################

mutDATA <- read.csv (file = 'SOURCE_ccletrimmed_namechanged_HRT18removed_HEC1Bremoved.csv')
colnames(mutDATA)[8] <- "case_id"
TRUNCATIONmutations = subset(mutDATA, mutation_type == 'Nonsense_Mutation') 

# Merge gene motif data with truncation mutations.
merged <- merge(TRUNCATIONmutations, genbankkinases, by = 'Gene', nomatch = 0)


# Extract Codon / Ref  from amino acid change.
# Note: Will match ref amino acid as well as codon (endCatalytic).
merged$Codon <- as.numeric(str_extract(merged$amino_acid_change, "[0-9]+"))
merged$RefAmino <- str_extract(merged$amino_acid_change, "[A-Z]")

# Extract predicted amino acid from each protein seq (calculated from Codon).
merged$codonConfirm <- substr(merged$Protein_Seq, merged$Codon, merged$Codon)

# Remove cases with more than one truncation mutation within the same sample (case_id).
mergedNODUPS <- merged[order(merged$Codon),]
mergedNODUPS$noMULTIPLEhits <- paste(mergedNODUPS$f, mergedNODUPS$case_id)
mergedNODUPS$transcriptDUP = !duplicated(mergedNODUPS$noMULTIPLEhits)
mergedNODUPS <- subset(mergedNODUPS, transcriptDUP == 'TRUE')

# Merge by Codon smaller than endCAT and Ref Amino acid is the same as locAmino.
truncSTRINGENT <- subset (mergedNODUPS, (Codon <= endCAT) & (RefAmino == codonConfirm))

# Rank matches by ascending endCAT result so that the shortest matching transcript is retained.
truncSTRINGENTSHORTEST <- truncSTRINGENT[order(truncSTRINGENT$endCAT),]

# Filter out duplicates of Gene, Codon and Sample (this leaves the shortest transcript).
truncSTRINGENTSHORTEST$duplicateREF <- paste(truncSTRINGENTSHORTEST$Gene, truncSTRINGENTSHORTEST$Codon, truncSTRINGENTSHORTEST$case_id, sep ='_')
truncSTRINGENTSHORTEST$mergedDUP= !duplicated(truncSTRINGENTSHORTEST$duplicateREF)
genefreqSHORTEST <- subset(truncSTRINGENTSHORTEST, mergedDUP == 'TRUE')
genefreqSHORTEST$transcriptFREQ <- paste(genefreqSHORTEST$f, genefreqSHORTEST$endCAT, sep = '_')
freqSTRINGENTSHORTEST <-table(genefreqSHORTEST$transcriptFREQ)
write.csv (freqSTRINGENTSHORTEST, file ='TEMP_freqTRANSCIPTSSHORTESTccleALL.csv')

# Length correct.
freqSHORTEST <- read.csv ('TEMP_freqTRANSCIPTSSHORTESTccleALL.csv')
colnames(freqSHORTEST)[2] <- "transcriptFREQ"
#extract gene and endCAT from each freq name
freqSHORTEST$Gene <- gsub    ( "[-].*$"     , "", freqSHORTEST$transcriptFREQ)
freqSHORTEST$endCAT <-  gsub ( ".*[_]"    , "", freqSHORTEST$transcriptFREQ)
freqSHORTEST$endCAT <- as.numeric(freqSHORTEST$endCAT)
#calculate length corrected freq for each transcript entry
freqSHORTEST$oneoverendCAT <- 1 / freqSHORTEST$endCAT
freqSHORTEST$score <- freqSHORTEST$Freq * freqSHORTEST$oneoverendCAT

# Add scores for each gene together to get total for shortest transcript (using dplyr).
freqADDEDSHORTEST <- ddply(freqSHORTEST,"Gene",numcolwise(sum))
write.csv (freqADDEDSHORTEST, file ='TEMP_lengthcorrectedfreq_ccleALLSHORTEST.csv')

# Add scores for each gene together to get total for longest transcript (using dplyr).
truncSTRINGENTLONGEST <- truncSTRINGENT[order(-truncSTRINGENT$endCAT),]
truncSTRINGENTLONGEST$duplicateREF <- paste(truncSTRINGENTLONGEST$Gene, truncSTRINGENTLONGEST$Codon, truncSTRINGENTLONGEST$case_id, sep ='_')
truncSTRINGENTLONGEST$mergedDUP= !duplicated(truncSTRINGENTLONGEST$duplicateREF)
genefreqLONGEST <- subset(truncSTRINGENTLONGEST, mergedDUP == 'TRUE')
genefreqLONGEST$transcriptFREQ <- paste(genefreqLONGEST$f, genefreqLONGEST$endCAT, sep = '_')
freqSTRINGENTLONGEST <-table(genefreqLONGEST$transcriptFREQ)
write.csv (freqSTRINGENTLONGEST, file ='TEMP_freqTRANSCIPTSLONGESTccleALL.csv')
freqLONGEST <- read.csv ('TEMP_freqTRANSCIPTSLONGESTccleALL.csv')
colnames(freqLONGEST)[2] <- "transcriptFREQ"
freqLONGEST$Gene <- gsub    ( "[-].*$"     , "", freqLONGEST$transcriptFREQ)
freqLONGEST$endCAT <-  gsub ( ".*[_]"    , "", freqLONGEST$transcriptFREQ)
freqLONGEST$endCAT <- as.numeric(freqLONGEST$endCAT)
freqLONGEST$oneoverendCAT <- 1 / freqLONGEST$endCAT
freqLONGEST$score <- freqLONGEST$Freq * freqLONGEST$oneoverendCAT
freqADDEDLONGEST <- ddply(freqLONGEST,"Gene",numcolwise(sum))
write.csv (freqADDEDLONGEST, file ='TEMP_lengthcorrectedfreq_ccleALLLONGEST.csv')

# Create a table to bind to CCLE genefreqSHORTEST to validate all mutations.
validateCCLE <- as.data.frame(genefreqSHORTEST$Gene)
colnames(validateCCLE)[1] <- "Gene"
validateCCLE$transcript <- genefreqSHORTEST$f
validateCCLE$mutation <- genefreqSHORTEST$amino_acid_change
validateCCLE$study <- genefreqSHORTEST$case_id
validateCCLE$case_id <- genefreqSHORTEST$case_id
validateCCLE$codon <- genefreqSHORTEST$Codon
validateCCLE$endCAT <- genefreqSHORTEST$endCAT
write.csv (validateCCLE, file = "Validate_CCLE.csv")


