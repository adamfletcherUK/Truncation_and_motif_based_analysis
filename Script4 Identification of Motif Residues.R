#################################################################################################################
# Script 4: Identification of Motif Residues                                                                    #
#                                                                                                               #
# Automatic search function to detect conserved functional motifs in the kinase domain of the kinome.           #
#################################################################################################################

# Intro text.
cat("A. Hudson, N. Stephenson, C. Li, E. Trotter, A. Fletcher, G. Katrona, P. Bieniasz-Krzywiec, M. Howell, C. Writh, C. Miller, J. Brognard. \n\n",
 "Functional Screening of Cancer Datasets Identifies Mutational Hotspots and Novel Targets in the Tumour-Suppressing Kinome. \n\n", 
 "Script requires stringr and plyr packages, additionally set the working directory to the folder containing the SOURCE.csv files \n\n", 
 sep="")


# Load script dependencies.
install.packages("stringr")
install.packages("plyr")
library("stringr", "plyr")


########################
# Load Genbank data and determine motifs positions
########################

# Load and merge definitive gene names with GENBANK.
genbank <- read.csv ( file = 'SOURCE_DEFINITIVE_GENBANK_TRIM.csv') #PROTEIN SEQ FOR ALL TRANSCRIPTS OF ALL GENES
names   <- read.csv (file = 'SOURCE_411kinases_genename.csv') # Loads genenames for the 411 kinases
genbankkinases = merge(genbank, names, by = 'Gene', nomatch = 0)


# Check how many unique kinases are present.
uniqkinase <- genbankkinases$Gene
unique(uniqkinase) 


# Extract known kinase motifs from the kinase sequence. 
# NOTE: Known kinase motifs are described in Manning et al. (2002), Science (https://www.ncbi.nlm.nih.gov/pubmed/12471243).
HRDmotifs  <- as.vector(genbankkinases$HRDmotif)
VAIKmotifs <- as.vector(genbankkinases$VAIKmotif)
DFGmotifs  <- as.vector(genbankkinases$DFGmotif)


# Count number of occurences of motif in each seq (using StringR package).
numberofDFG  <- str_count(genbankkinases$Protein_Seq,  DFGmotifs)
numberofVAIK <- str_count(genbankkinases$Protein_Seq, VAIKmotifs)
numberofHRD  <- str_count(genbankkinases$Protein_Seq,  HRDmotifs)
genbankkinases$numberofVAIK <- numberofVAIK
genbankkinases$numberofHRD  <- numberofHRD
genbankkinases$numberofDFG  <- numberofDFG


# Removes all entries with no VAIK, HRD or DFG motifs.
genbankkinases <- subset(genbankkinases, (numberofVAIK != 0))
genbankkinases <- subset(genbankkinases, (numberofHRD  != 0))
genbankkinases <- subset(genbankkinases, (numberofDFG  != 0))


# Define the correct motif to use.
# NOTE: SOURCE_motif2use is a file with the correct motif to use when there are multiple. Correct motif 
motifToUse     <- read.csv ('SOURCE_motif2use.csv')
genbankkinases<- merge    (genbankkinases, motifToUse, by = c('f','Gene', 'Length', 'Protein_Seq', 'VAIKmotif', 'HRDmotif', 'DFGmotif'))


# Locate all motifs.
HRDmotifs  <- as.vector(genbankkinases$HRDmotif)
VAIKmotifs <- as.vector(genbankkinases$VAIKmotif)
DFGmotifs  <- as.vector(genbankkinases$DFGmotif)


# Locate the position of the first VAIK motif .
# Note: One kinase is found to have two VAIK motifs through this screen - MAP4K1. Upon sequence inspection the first was found to be the catalytic VAIK.
firstVAIK <- str_locate(genbankkinases$Protein_Seq, VAIKmotifs)
firstVAIK <- as.data.frame(firstVAIK)
genbankkinases$firstVAIK <- as.numeric(firstVAIK$start)
genbankkinases$actual_lysine <- as.numeric(genbankkinases$firstVAIK + 3)


########################
# Locate all HRD motifs in the peptide sequence
########################

# Locate the first HRD Motif.
firstHRD <- str_locate(genbankkinases$Protein_Seq, HRDmotifs)
firstHRD <- as.data.frame(firstHRD)
genbankkinases$firstHRD <- firstHRD$start
firstHRD <- as.numeric(firstHRD$start)
endofstring = as.numeric(firstHRD +1000000)
#Locates first HRD motif and outputs the remaining protein sequence
beyondfirstHRD <- substr(genbankkinases$Protein_Seq, firstHRD + 1, endofstring) #prints the rest of the sequence after the motif


# Locates second HRD motif in protein sequence.
secondHRD <- str_locate(beyondfirstHRD, HRDmotifs)
secondHRD <- as.data.frame(secondHRD)
genbankkinases$secondHRD <- secondHRD$start
secondHRD = as.numeric(genbankkinases$secondHRD + firstHRD)
genbankkinases$secondHRD <- secondHRD


# Prints the rest of the sequence after the 2nd HRD motif.
beyondsecondHRD <- substr(genbankkinases$Protein_Seq, secondHRD + 1, endofstring) 


# Locates third HRD motif in protein sequence.
thirdHRD <- str_locate(beyondsecondHRD, HRDmotifs)
thirdHRD <- as.data.frame(thirdHRD)
genbankkinases$thirdHRD <- thirdHRD$start
thirdHRD = as.numeric(genbankkinases$thirdHRD + secondHRD)
genbankkinases$thirdHRD <- thirdHRD


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
# Chosing the correct motifs to use
########################

# Determine the correct HRD to use.
genbankkinases$actual_HRD_H <- ifelse (genbankkinases$HRDtouse == '1', genbankkinases$firstHRD, ifelse (genbankkinases$HRDtouse == '2', genbankkinases$secondHRD , genbankkinases$thirdHRD))
genbankkinases$actual_HRD_R <- genbankkinases$actual_HRD_H + 1
genbankkinases$actual_HRD_D <- genbankkinases$actual_HRD_H + 2


# Now finds the HRD+5 amino acid based upon the HRD information.
genbankkinases$actual_HRDp5 <- genbankkinases$actual_HRD_D + 5
genbankkinases$actual_HRDp5_letter     <- str_sub(genbankkinases$Protein_Seq, genbankkinases$actual_HRDp5, genbankkinases$actual_HRDp5)
write.csv (genbankkinases, file ='TEMP_HRDp5check.csv')

########################
# Find the HRD-6 motif
########################

# Find sequence location around the APE motif.
genbankkinases$HRD_m14loc   <- genbankkinases$actual_HRD_H - 14
genbankkinases$HRD_m13loc   <- genbankkinases$actual_HRD_H - 13
genbankkinases$HRD_m12loc   <- genbankkinases$actual_HRD_H - 12
genbankkinases$HRD_m11loc   <- genbankkinases$actual_HRD_H - 11
genbankkinases$HRD_m10loc   <- genbankkinases$actual_HRD_H - 10
genbankkinases$HRD_m9loc   <- genbankkinases$actual_HRD_H - 9
genbankkinases$HRD_m8loc   <- genbankkinases$actual_HRD_H - 8
genbankkinases$HRD_m7loc   <- genbankkinases$actual_HRD_H - 7
genbankkinases$HRD_m6loc   <- genbankkinases$actual_HRD_H - 6
genbankkinases$HRD_m5loc   <- genbankkinases$actual_HRD_H - 5
genbankkinases$HRD_m4loc   <- genbankkinases$actual_HRD_H - 4
genbankkinases$HRD_m3loc   <- genbankkinases$actual_HRD_H - 3
genbankkinases$HRD_m2loc   <- genbankkinases$actual_HRD_H - 2
genbankkinases$HRD_m1loc   <- genbankkinases$actual_HRD_H - 1
genbankkinases$HRD_H1loc   <- genbankkinases$actual_HRD_H
genbankkinases$HRD_R1loc   <- genbankkinases$actual_HRD_H +1
genbankkinases$HRD_D1loc   <- genbankkinases$actual_HRD_H +2


# Find amino acid identity for above sequence locations.
genbankkinases$HRD_m14       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m14loc , genbankkinases$HRD_m14loc)
genbankkinases$HRD_m13       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m13loc , genbankkinases$HRD_m13loc)
genbankkinases$HRD_m12       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m12loc , genbankkinases$HRD_m12loc)
genbankkinases$HRD_m11       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m11loc , genbankkinases$HRD_m11loc)
genbankkinases$HRD_m10       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m10loc , genbankkinases$HRD_m10loc)
genbankkinases$HRD_m9       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m9loc , genbankkinases$HRD_m9loc)
genbankkinases$HRD_m8       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m8loc , genbankkinases$HRD_m8loc)
genbankkinases$HRD_m7       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m7loc , genbankkinases$HRD_m7loc)
genbankkinases$HRD_m6       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m6loc , genbankkinases$HRD_m6loc)
genbankkinases$HRD_m5       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m5loc , genbankkinases$HRD_m5loc)
genbankkinases$HRD_m4       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m4loc , genbankkinases$HRD_m4loc)
genbankkinases$HRD_m3       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m3loc , genbankkinases$HRD_m3loc)
genbankkinases$HRD_m2       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m2loc , genbankkinases$HRD_m2loc)
genbankkinases$HRD_m1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_m1loc , genbankkinases$HRD_m1loc)
genbankkinases$HRD_H1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_H1loc , genbankkinases$HRD_H1loc)
genbankkinases$HRD_R1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_R1loc , genbankkinases$HRD_R1loc)
genbankkinases$HRD_D1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRD_D1loc , genbankkinases$HRD_D1loc)


# Look for a H in HRD-6, then -7, -8, -4, -9, -5 and -10 
# NOTE: If it cannot find one in these locations, it will output HRD-6

genbankkinases$HRDm6_loc <- ifelse (genbankkinases$HRD_m6 != 'H', genbankkinases$HRD_m7loc, genbankkinases$HRD_m6loc)
genbankkinases$HRDm6 <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRDm6_loc , genbankkinases$HRDm6_loc)

genbankkinases$HRDm6_loc <- ifelse (genbankkinases$HRDm6 != 'H', genbankkinases$HRD_m8loc, genbankkinases$HRDm6_loc)
genbankkinases$HRDm6 <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRDm6_loc , genbankkinases$HRDm6_loc)

genbankkinases$HRDm6_loc <- ifelse (genbankkinases$HRDm6 != 'H', genbankkinases$HRD_m4loc, genbankkinases$HRDm6_loc)
genbankkinases$HRDm6 <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRDm6_loc , genbankkinases$HRDm6_loc) 

genbankkinases$HRDm6_loc <- ifelse (genbankkinases$HRDm6 != 'H', genbankkinases$HRD_m9loc, genbankkinases$HRDm6_loc)
genbankkinases$HRDm6 <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRDm6_loc , genbankkinases$HRDm6_loc) 

genbankkinases$HRDm6_loc <- ifelse (genbankkinases$HRDm6 != 'H', genbankkinases$HRD_m5loc, genbankkinases$HRDm6_loc)
genbankkinases$HRDm6 <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRDm6_loc , genbankkinases$HRDm6_loc) 

genbankkinases$HRDm6_loc <- ifelse (genbankkinases$HRDm6 != 'H', genbankkinases$HRD_m10loc, genbankkinases$HRDm6_loc)
genbankkinases$HRDm6 <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRDm6_loc , genbankkinases$HRDm6_loc) 

genbankkinases$HRDm6_loc <- ifelse (genbankkinases$HRDm6 == 'H', genbankkinases$HRDm6_loc, genbankkinases$HRD_m6loc)
genbankkinases$HRDm6 <- str_sub(genbankkinases$Protein_Seq, genbankkinases$HRDm6_loc , genbankkinases$HRDm6_loc)


########################
# Locate APE and GxGxxG motif information
########################

# Determine an APE fragment based upon the DFG motif.
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


# Find sequence location around the APE motif.
genbankkinases$APE_m14loc   <- genbankkinases$loc_actual_APE - 14
genbankkinases$APE_m13loc   <- genbankkinases$loc_actual_APE - 13
genbankkinases$APE_m12loc   <- genbankkinases$loc_actual_APE - 12
genbankkinases$APE_m11loc   <- genbankkinases$loc_actual_APE - 11
genbankkinases$APE_m10loc   <- genbankkinases$loc_actual_APE - 10
genbankkinases$APE_m9loc   <- genbankkinases$loc_actual_APE - 9
genbankkinases$APE_m8loc   <- genbankkinases$loc_actual_APE - 8
genbankkinases$APE_m7loc   <- genbankkinases$loc_actual_APE - 7
genbankkinases$APE_m6loc   <- genbankkinases$loc_actual_APE - 6
genbankkinases$APE_m5loc   <- genbankkinases$loc_actual_APE - 5
genbankkinases$APE_m4loc   <- genbankkinases$loc_actual_APE - 4
genbankkinases$APE_m3loc   <- genbankkinases$loc_actual_APE - 3
genbankkinases$APE_m2loc   <- genbankkinases$loc_actual_APE - 2
genbankkinases$APE_m1loc   <- genbankkinases$loc_actual_APE - 1
genbankkinases$APE_A1loc   <- genbankkinases$loc_actual_APE
genbankkinases$APE_P1loc   <- genbankkinases$loc_actual_APE + 1
genbankkinases$APE_E1loc   <- genbankkinases$loc_actual_APE + 2


# Find amino acid identity for above sequence locations.
genbankkinases$APE_m14       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m14loc , genbankkinases$APE_m14loc)
genbankkinases$APE_m13       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m13loc , genbankkinases$APE_m13loc)
genbankkinases$APE_m12       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m12loc , genbankkinases$APE_m12loc)
genbankkinases$APE_m11       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m11loc , genbankkinases$APE_m11loc)
genbankkinases$APE_m10       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m10loc , genbankkinases$APE_m10loc)
genbankkinases$APE_m9       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m9loc , genbankkinases$APE_m9loc)
genbankkinases$APE_m8       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m8loc , genbankkinases$APE_m8loc)
genbankkinases$APE_m7       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m7loc , genbankkinases$APE_m7loc)
genbankkinases$APE_m6       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m6loc , genbankkinases$APE_m6loc)
genbankkinases$APE_m5       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m5loc , genbankkinases$APE_m5loc)
genbankkinases$APE_m4       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m4loc , genbankkinases$APE_m4loc)
genbankkinases$APE_m3       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m3loc , genbankkinases$APE_m3loc)
genbankkinases$APE_m2       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m2loc , genbankkinases$APE_m2loc)
genbankkinases$APE_m1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_m1loc , genbankkinases$APE_m1loc)
genbankkinases$APE_A1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_A1loc , genbankkinases$APE_A1loc)
genbankkinases$APE_P1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_P1loc , genbankkinases$APE_P1loc)
genbankkinases$APE_E1       <- str_sub(genbankkinases$Protein_Seq, genbankkinases$APE_E1loc , genbankkinases$APE_E1loc)


# Output the final APE motif sequence information.
genbankkinases$APE_MOTIF <- paste(genbankkinases$APE_A1, genbankkinases$APE_P1, genbankkinases$APE_E1)
genbankkinases$PE_MOTIF <- paste(genbankkinases$APE_P1, genbankkinases$APE_E1)


########################
# Finds the APE-6 motif
########################

# Locate a G in APE-6, then APE-9, APE-5, APE-8, APE-10 and APE-11. 
# NOTE: If it cannot find one in these locations, it will output APE-6

genbankkinases$correctedG_loc <- ifelse (genbankkinases$APE_m6 != 'G', genbankkinases$APE_m9loc, genbankkinases$APE_m6loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m5loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m8loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m10loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m11loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG != 'G', genbankkinases$APE_m7loc, genbankkinases$correctedG_loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

genbankkinases$correctedG_loc <- ifelse (genbankkinases$correctedG == 'G', genbankkinases$correctedG_loc, genbankkinases$APE_m6loc)
genbankkinases$correctedG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$correctedG_loc , genbankkinases$correctedG_loc)

#end of kinase domain = APE + defined end
genbankkinases$endCAT <- as.numeric(genbankkinases$APE_E1loc)


########################
# Finds the GxGxxG motif
########################

# Find GxGxxG.
GxGfragstart <- genbankkinases$actual_lysine - 45
genbankkinases$GxGfragstart <- as.numeric(GxGfragstart)
genbankkinases$GxGfragstart <- ifelse(genbankkinases$GxGfragstart < 1, 1,genbankkinases$GxGfragstart )

GxGfragend <- genbankkinases$actual_lysine - 4
genbankkinases$GxGfragend <- as.numeric(GxGfragend)
GxGfrag <- str_sub(genbankkinases$Protein_Seq, genbankkinases$GxGfragstart, genbankkinases$GxGfragend )
genbankkinases$GxGfrag <- GxGfrag


# Find GxGxxG.
numGxGxxG <- str_count(GxGfrag,  'G[A-Z]G[A-Z]{2}G')
genbankkinases$numGxGxxG <- as.numeric(numGxGxxG)
GxGxxGLoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]{2}G')
GxGxxGLoc <- as.data.frame(GxGxxGLoc)
genbankkinases$GxGxxGLoc <- as.numeric(GxGxxGLoc$start)
genbankkinases$locGxGxxGprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxxGLoc - 1


# Find GxGxF.
numGxGxF <- str_count(GxGfrag,  'G[A-Z]G[A-Z]F')
genbankkinases$numGxGxF <- as.numeric(numGxGxF)
GxGxFLoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]F')
GxGxFLoc <- as.data.frame(GxGxFLoc)
genbankkinases$GxGxFLoc <- as.numeric(GxGxFLoc$start)
genbankkinases$locGxGxFprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxFLoc - 1


# Find GxGxY.
numGxGxY <- str_count(GxGfrag,  'G[A-Z]G[A-Z]Y')
genbankkinases$numGxGxY <- as.numeric(numGxGxY)
GxGxYLoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]Y')
GxGxYLoc <- as.data.frame(GxGxYLoc)
genbankkinases$GxGxYLoc <- as.numeric(GxGxYLoc$start)
genbankkinases$locGxGxYprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxYLoc - 1


# Find GxFG.
numGxFG <- str_count(GxGfrag,  'G[A-Z]FG')
genbankkinases$numGxFG <- as.numeric(numGxFG)
GxFGLoc <- str_locate(GxGfrag, 'G[A-Z]FG')
GxFGLoc <- as.data.frame(GxFGLoc)
genbankkinases$GxFGLoc <- as.numeric(GxFGLoc$start)
genbankkinases$locGxFGprotein <- genbankkinases$GxGfragstart + genbankkinases$GxFGLoc - 3
#above is minus 3 rather than minus 1 to take account of the GxFG starting 2 on from GxGxxG


# Find GxGxxA.
numGxGxxA <- str_count(GxGfrag,  'G[A-Z]G[A-Z]{2}A')
genbankkinases$numGxGxxA <- as.numeric(numGxGxxA)
GxGxxALoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]{2}A')
GxGxxALoc <- as.data.frame(GxGxxALoc)
genbankkinases$GxGxxALoc <- as.numeric(GxGxxALoc$start)
genbankkinases$locGxGxxAprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxxALoc - 1


# Find GxGxxS.
numGxGxxS <- str_count(GxGfrag,  'G[A-Z]G[A-Z]{2}S')
genbankkinases$numGxGxxS <- as.numeric(numGxGxxS)
GxGxxSLoc <- str_locate(GxGfrag, 'G[A-Z]G[A-Z]{2}S')
GxGxxSLoc <- as.data.frame(GxGxxSLoc)
genbankkinases$GxGxxSLoc <- as.numeric(GxGxxSLoc$start)
genbankkinases$locGxGxxSprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGxxSLoc - 1


# Find SxGxxG.
numSxGxxG <- str_count(GxGfrag,  'S[A-Z]G[A-Z]{2}G')
genbankkinases$numSxGxxG <- as.numeric(numSxGxxG)
SxGxxGLoc <- str_locate(GxGfrag, 'S[A-Z]G[A-Z]{2}G')
SxGxxGLoc <- as.data.frame(SxGxxGLoc)
genbankkinases$SxGxxGLoc <- as.numeric(SxGxxGLoc$start)
genbankkinases$locSxGxxGprotein <- genbankkinases$GxGfragstart + genbankkinases$SxGxxGLoc - 1


# Find GxxG.
numGxxG <- str_count(GxGfrag,  'G[A-Z]{2}G')
genbankkinases$numGxxG <- as.numeric(numGxxG)
GxxGLoc <- str_locate(GxGfrag, 'G[A-Z]{2}G')
GxxGLoc <- as.data.frame(GxxGLoc)
genbankkinases$GxxGLoc <- as.numeric(GxxGLoc$start)
genbankkinases$locGxxGprotein <- genbankkinases$GxGfragstart + genbankkinases$GxxGLoc - 3


# Find G4xG.
numG4xG <- str_count(GxGfrag,  'G[A-Z]{4}G')
genbankkinases$numG4xG <- as.numeric(numG4xG)
G4xGLoc <- str_locate(GxGfrag, 'G[A-Z]{4}G')
G4xGLoc <- as.data.frame(G4xGLoc)
genbankkinases$G4xGLoc <- as.numeric(G4xGLoc$start)
genbankkinases$locG4xGprotein <- genbankkinases$GxGfragstart + genbankkinases$G4xGLoc - 1


# Find GxG.
numGxG <- str_count(GxGfrag,  'G[A-Z]G')
genbankkinases$numGxG <- as.numeric(numGxG)
GxGLoc <- str_locate(GxGfrag, 'G[A-Z]G')
GxGLoc <- as.data.frame(GxGLoc)
genbankkinases$GxGLoc <- as.numeric(GxGLoc$start)
genbankkinases$locGxGprotein <- genbankkinases$GxGfragstart + genbankkinases$GxGLoc - 1


# Find GxTxF (locates PIK3R4).
###################################### Did you check??? #########################################################
numGxTxF <- str_count(GxGfrag,  'G[A-Z]T[A-Z]F')
genbankkinases$numGxTxF <- as.numeric(numGxTxF)
GxTxFLoc <- str_locate(GxGfrag, 'G[A-Z]T[A-Z]F')
GxTxFLoc <- as.data.frame(GxTxFLoc)
genbankkinases$GxTxFLoc <- as.numeric(GxTxFLoc$start)
genbankkinases$locGxTxFprotein <- genbankkinases$GxGfragstart + genbankkinases$GxTxFLoc - 1
#################################################################################################################

# Find GG (only use as the last resort to capture AMHRH2 and ANKK1).
numGG <- str_count(GxGfrag, 'GG')
genbankkinases$numGG <- as.numeric(numGG)
GGLoc <- str_locate(GxGfrag, 'GG')
GGLoc <- as.data.frame(GGLoc)
genbankkinases$GGLoc <- as.numeric(GGLoc$start)
genbankkinases$locGGprotein <- genbankkinases$GxGfragstart + genbankkinases$GGLoc - 3


########################
# Determine the correct GxG
########################

genbankkinases$loc_actual_gxg_A <- ifelse (genbankkinases$numGxGxxG <1, genbankkinases$locGxGxFprotein, genbankkinases$locGxGxxGprotein)
genbankkinases$loc_actual_gxg_B <- genbankkinases$loc_actual_gxg_A
genbankkinases$loc_actual_gxg_C   <- ifelse (genbankkinases$loc_actual_gxg_A %in% NA, genbankkinases$locGxFGprotein, genbankkinases$loc_actual_gxg_B)
genbankkinases$loc_actual_gxg_D <- genbankkinases$loc_actual_gxg_C
genbankkinases$loc_actual_gxg_E   <- ifelse (genbankkinases$loc_actual_gxg_C %in% NA, genbankkinases$locGxGxYprotein, genbankkinases$loc_actual_gxg_D)
genbankkinases$loc_actual_gxg_F <- genbankkinases$loc_actual_gxg_E
genbankkinases$loc_actual_gxg_G   <- ifelse (genbankkinases$loc_actual_gxg_E %in% NA, genbankkinases$locGxGxxAprotein, genbankkinases$loc_actual_gxg_F)
genbankkinases$loc_actual_gxg_H <- genbankkinases$loc_actual_gxg_G
genbankkinases$loc_actual_gxg_I   <- ifelse (genbankkinases$loc_actual_gxg_G %in% NA, genbankkinases$locGxGxxSprotein, genbankkinases$loc_actual_gxg_H)
genbankkinases$loc_actual_gxg_J <- genbankkinases$loc_actual_gxg_I
genbankkinases$loc_actual_gxg_K   <- ifelse (genbankkinases$loc_actual_gxg_I %in% NA, genbankkinases$locSxGxxGprotein, genbankkinases$loc_actual_gxg_J)
genbankkinases$loc_actual_gxg_L <- genbankkinases$loc_actual_gxg_K
genbankkinases$loc_actual_gxg_M   <- ifelse (genbankkinases$loc_actual_gxg_K %in% NA, genbankkinases$locGxxGprotein, genbankkinases$loc_actual_gxg_L)
genbankkinases$loc_actual_gxg_N <- genbankkinases$loc_actual_gxg_M
genbankkinases$loc_actual_gxg_O   <- ifelse (genbankkinases$loc_actual_gxg_M %in% NA, genbankkinases$locG4xGprotein, genbankkinases$loc_actual_gxg_N)
genbankkinases$loc_actual_gxg_P <- genbankkinases$loc_actual_gxg_O
genbankkinases$loc_actual_gxg_Q   <- ifelse (genbankkinases$loc_actual_gxg_O %in% NA, genbankkinases$locGxGprotein, genbankkinases$loc_actual_gxg_P)
genbankkinases$loc_actual_gxg_R <- genbankkinases$loc_actual_gxg_Q
genbankkinases$loc_actual_gxg_S   <- ifelse (genbankkinases$loc_actual_gxg_Q %in% NA, genbankkinases$locGxTxFprotein, genbankkinases$loc_actual_gxg_R)
genbankkinases$loc_actual_gxg_T <- genbankkinases$loc_actual_gxg_S
genbankkinases$loc_actual_gxg   <- ifelse (genbankkinases$loc_actual_gxg_S %in% NA, genbankkinases$locGGprotein, genbankkinases$loc_actual_gxg_T)


# Extract GxGxxG motif and location/aa for each.
startGxG <- genbankkinases$loc_actual_gxg
endGxG <- genbankkinases$loc_actual_gxg + 5
genbankkinases$GxGxxGmotif <- str_sub(genbankkinases$Protein_Seq, startGxG, endGxG)


# Amino acid residues
genbankkinases$GxGxxG_G1     <- str_sub(genbankkinases$GxGxxGmotif, 1, 1)
genbankkinases$GxGxxG_X1     <- str_sub(genbankkinases$GxGxxGmotif, 2, 2)
genbankkinases$GxGxxG_G2     <- str_sub(genbankkinases$GxGxxGmotif, 3, 3)
genbankkinases$GxGxxG_X2     <- str_sub(genbankkinases$GxGxxGmotif, 4, 4)
genbankkinases$GxGxxG_X3     <- str_sub(genbankkinases$GxGxxGmotif, 5, 5)
genbankkinases$GxGxxG_G3     <- str_sub(genbankkinases$GxGxxGmotif, 6, 6)


# Amino acid positions
genbankkinases$GxGxxG_G1loc   <- genbankkinases$loc_actual_gxg
genbankkinases$GxGxxG_X1loc   <- genbankkinases$loc_actual_gxg + 1
genbankkinases$GxGxxG_G2loc   <- genbankkinases$loc_actual_gxg + 2
genbankkinases$GxGxxG_X2loc   <- genbankkinases$loc_actual_gxg + 3
genbankkinases$GxGxxG_X3loc   <- genbankkinases$loc_actual_gxg + 4
genbankkinases$GxGxxG_G3loc   <- genbankkinases$loc_actual_gxg + 5

genbankkinases$GXGXXG_MOTIF <- paste(genbankkinases$GxGxxG_G1, genbankkinases$GxGxxG_G2, genbankkinases$GxGxxG_G3)
genbankkinases$CATALYTIC_FRAG <- str_sub(genbankkinases$Protein_Seq, genbankkinases$GxGxxG_G1loc, genbankkinases$APE_E1loc )


########################
# Finds the Salt Bridge Information
########################

# Determines a Salt Bridge fragment based on Lysine.
# NOTE: This is used to determine the Salt Bridge E.
genbankkinases$saltbridgeEfragSTART <- genbankkinases$actual_lysine + 5
genbankkinases$saltbridgeEfragEND <- genbankkinases$actual_lysine + 30
genbankkinases$saltbridgeEfrag <- str_sub(genbankkinases$Protein_Seq, genbankkinases$saltbridgeEfragSTART, genbankkinases$saltbridgeEfragEND )
saltbridgeEfrag <- genbankkinases$saltbridgeEfrag


# Find ExxxL.
numExxxL <- str_count(saltbridgeEfrag,  'E[A-Z]{3}L')
genbankkinases$numExxxL <- as.numeric(numExxxL)
ExxxLLoc <- str_locate(saltbridgeEfrag,  'E[A-Z]{3}L')
ExxxLLoc <- as.data.frame(ExxxLLoc)
genbankkinases$ExxxLLoc <- as.numeric(ExxxLLoc$start)
genbankkinases$ExxxLLocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$ExxxLLoc - 1


#  Find ExxI.
numExxI <- str_count(saltbridgeEfrag,  'E[A-Z]{2}I')
genbankkinases$numExxI <- as.numeric(numExxI)
ExxILoc <- str_locate(saltbridgeEfrag,  'E[A-Z]{2}I')
ExxILoc <- as.data.frame(ExxILoc)
genbankkinases$ExxILoc <- as.numeric(ExxILoc$start)
genbankkinases$ExxILocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$ExxILoc - 1


# Find ExxxM.
numExxxM <- str_count(saltbridgeEfrag,  'E[A-Z]{3}M')
genbankkinases$numExxxM <- as.numeric(numExxxM)
ExxxMLoc <- str_locate(saltbridgeEfrag,  'E[A-Z]{3}M')
ExxxMLoc <- as.data.frame(ExxxMLoc)
genbankkinases$ExxxMLoc <- as.numeric(ExxxMLoc$start)
genbankkinases$ExxxMLocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$ExxxMLoc - 1


# Find ExxL.
numExxL <- str_count(saltbridgeEfrag,  'E[A-Z]{2}L')
genbankkinases$numExxL <- as.numeric(numExxL)
ExxLLoc <- str_locate(saltbridgeEfrag,  'E[A-Z]{2}L')
ExxLLoc <- as.data.frame(ExxLLoc)
genbankkinases$ExxLLoc <- as.numeric(ExxLLoc$start)
genbankkinases$ExxLLocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$ExxLLoc - 1


# Find QxxxE.
numQxxxE <- str_count(saltbridgeEfrag,  'Q[A-Z]{3}E')
genbankkinases$numQxxxE <- as.numeric(numQxxxE)
QxxxELoc <- str_locate(saltbridgeEfrag,  'Q[A-Z]{3}E')
QxxxELoc <- as.data.frame(QxxxELoc)
genbankkinases$QxxxELoc <- as.numeric(QxxxELoc$start)
genbankkinases$QxxxELocprotein <- genbankkinases$saltbridgeEfragSTART + genbankkinases$QxxxELoc + 3


# Choose correct saltbridge.
genbankkinases$loc_actual_saltbridge_A <- ifelse (genbankkinases$numExxxL <1, genbankkinases$ExxILocprotein, genbankkinases$ExxxLLocprotein)
genbankkinases$loc_actual_saltbridge_B <- genbankkinases$loc_actual_saltbridge_A
genbankkinases$loc_actual_saltbridge_C   <- ifelse (genbankkinases$loc_actual_saltbridge_A %in% NA, genbankkinases$ExxxMLocprotein, genbankkinases$loc_actual_saltbridge_B)
genbankkinases$loc_actual_saltbridge_D <- genbankkinases$loc_actual_saltbridge_C
genbankkinases$loc_actual_saltbridge_E   <- ifelse (genbankkinases$loc_actual_saltbridge_C %in% NA, genbankkinases$ExxLLocprotein, genbankkinases$loc_actual_saltbridge_D)
genbankkinases$loc_actual_saltbridge_F <- genbankkinases$loc_actual_saltbridge_E
genbankkinases$loc_actual_saltbridge   <- ifelse (genbankkinases$loc_actual_saltbridge_E %in% NA, genbankkinases$QxxxELocprotein, genbankkinases$loc_actual_saltbridge_F)
genbankkinases$checkE <- str_sub(genbankkinases$Protein_Seq, genbankkinases$loc_actual_saltbridge, genbankkinases$loc_actual_saltbridge)


# Outputs the final Genbank kinase information as a TEMP backup file.
write.csv (genbankkinases, file ='TEMP_genbankkinases.csv')


