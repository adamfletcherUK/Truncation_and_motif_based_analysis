#################################################################################################################
# Script 5: Identifies mutations within TCGA database that occur in hotspot kinase motifs specified in Script 3.#
#                                                                                                               #
# Automatic search function to detect mutations in the TCGA database that occur within hotspot kinase residues  #
# as defined by Script 3.                                                                                       #
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
# Merge Genbank Data with the TCGA Database and formatting
########################

# Merge with combined TCGA data.
comboDATA <- read.csv (file = 'SOURCE_allTCGA_jan2016.csv') 
kinaseLoc <- genbankkinases


# Removes all non missense mutations.
comboDATA = subset(comboDATA, mutation_type != '3UTR')
comboDATA = subset(comboDATA, mutation_type != '5UTR')
comboDATA = subset(comboDATA, mutation_type != 'Intron')
comboDATA = subset(comboDATA, mutation_type != 'Splice_Site_Ins')
comboDATA = subset(comboDATA, mutation_type != 'Splice_Site_Del')
comboDATA = subset(comboDATA, mutation_type != 'Splice_Site_SNP')
comboDATA = subset(comboDATA, mutation_type != 'Silent')
comboDATA = subset(comboDATA, mutation_type != 'Nonsense_Mutation')
comboDATA = subset(comboDATA, mutation_type != 'Frame_Shift_Del')
comboDATA = subset(comboDATA, mutation_type != 'Frame_Shift_Ins')
comboDATA = subset(comboDATA, mutation_type != 'In_Frame_Del')
comboDATA = subset(comboDATA, mutation_type != 'In_Frame_Ins')
comboDATA = subset(comboDATA, mutation_type != 'outofframe')
comboDATA = subset(comboDATA, mutation_type != 'Flank')
comboDATA = subset(comboDATA, mutation_type != 'Splice_Site')
comboDATA = subset(comboDATA, mutation_type != 'De_novo_Start_InFrame')
comboDATA = subset(comboDATA, mutation_type != 'Stop_Codon_DNP')
comboDATA = subset(comboDATA, mutation_type != 'Stop_Codon_Ins')


# Extract the mutation position data from comboDATA.
CodonCombo <- str_extract(comboDATA$amino_acid_change, "[0-9]+")
comboDATA$Codon <- CodonCombo
RefAmino <- str_extract(comboDATA$amino_acid_change, "[A-Z]")
comboDATA$AA <- RefAmino
VarAmino <- str_sub(comboDATA$amino_acid_change, -1, -1)
comboDATA$VarAA <- VarAmino


########################
# Calculates the number of mutations occuring within hotspot residues from the TCGA database. 
########################

# GXGXXG_G1.
kinaseLoc$Codon <- kinaseLoc$GxGxxG_G1loc
kinaseLoc$AA <- kinaseLoc$GxGxxG_G1
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_G1 <- merged
write.csv (merged_G1, file='outputs/tcga/merged_GxGxxG_G1_all.csv')


# GXGXXG_G2.
kinaseLoc$Codon <- kinaseLoc$GxGxxG_G2loc
kinaseLoc$AA <- kinaseLoc$GxGxxG_G2
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_G2 <- merged
write.csv (merged_G2, file='outputs/tcga/merged_GxGxxG_G2_all.csv')


# GXGXXG_G3.
kinaseLoc$Codon <- kinaseLoc$GxGxxG_G3loc
kinaseLoc$AA <- kinaseLoc$GxGxxG_G3
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_G3 <- merged
write.csv (merged_G3, file='outputs/tcga/merged_GxGxxG_G3_all.csv')


# Critical lysine.
kinaseLoc$Codon <- kinaseLoc$actual_lysine
kinaseLoc$AA <- 'K'
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_K <- merged
write.csv (merged_K, file='outputs/tcga/critical_lysine_all.csv')


# DFG_D.
kinaseLoc$Codon <- kinaseLoc$actual_DFG_D
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_DFG_D, kinaseLoc$actual_DFG_D)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_DFG_D <- merged
write.csv (merged_DFG_D, file='outputs/tcga/merged_DFG_D_all.csv')


# DFG_F.
kinaseLoc$Codon <- kinaseLoc$actual_DFG_F
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_DFG_F, kinaseLoc$actual_DFG_F)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_DFG_F <- merged
write.csv (merged_DFG_F, file='outputs/tcga/merged_DFG_F_all.csv')


# DFG_G.
kinaseLoc$Codon <- kinaseLoc$actual_DFG_G
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_DFG_G, kinaseLoc$actual_DFG_G)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_DFG_G <- merged
write.csv (merged_DFG_G, file='outputs/tcga/merged_DFG_G_all.csv')


# HRD_H.
kinaseLoc$Codon <- kinaseLoc$actual_HRD_H
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_HRD_H, kinaseLoc$actual_HRD_H)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRD_H <- merged
write.csv (merged_HRD_H, file='outputs/tcga/merged_HRD_H_all.csv')


# HRD_R.
kinaseLoc$Codon <- kinaseLoc$actual_HRD_R
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_HRD_R, kinaseLoc$actual_HRD_R)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRD_R <- merged
write.csv (merged_HRD_R, file='outputs/tcga/merged_HRD_R_all.csv')


# HRD_D.
kinaseLoc$Codon <- kinaseLoc$actual_HRD_D
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_HRD_D, kinaseLoc$actual_HRD_D)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRD_D <- merged
write.csv (merged_HRD_D, file='outputs/tcga/merged_HRD_D_all.csv')


# HRD+5.
kinaseLoc$Codon <- kinaseLoc$actual_HRD_D + 5
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$actual_HRD_D, kinaseLoc$actual_HRD_D)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRD_Dp5 <- merged
write.csv (merged_HRD_Dp5, file='outputs/tcga/merged_HRD_Dp5_all.csv')


# HRD-6.
kinaseLoc$Codon <- kinaseLoc$HRDm6_loc
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$HRDm6_loc, kinaseLoc$HRDm6_loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_HRDm6 <- merged
write.csv (merged_HRDm6, file='outputs/tcga/merged_HRDm6_all.csv')


# APE_A.
kinaseLoc$Codon <- kinaseLoc$APE_A1loc 
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$APE_A1loc, kinaseLoc$APE_A1loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_APE_A <- merged
write.csv (merged_APE_A, file='outputs/tcga/merged_APE_A_all.csv')


# APE_P.
kinaseLoc$Codon <- kinaseLoc$APE_P1loc 
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$APE_P1loc, kinaseLoc$APE_P1loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_APE_P <- merged
write.csv (merged_APE_P, file='outputs/tcga/merged_APE_P_all.csv')


# APE_E.
kinaseLoc$Codon <- kinaseLoc$APE_E1loc 
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$APE_E1loc, kinaseLoc$APE_E1loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_APE_E <- merged
write.csv (merged_APE_E, file='outputs/tcga/merged_APE_E_all.csv')


# Corrected_G.
kinaseLoc$Codon <- kinaseLoc$correctedG_loc
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$correctedG_loc, kinaseLoc$correctedG_loc)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_correctedG <- merged
write.csv (merged_correctedG, file='outputs/tcga/merged_correctedG_all.csv')


# Saltbridge_E
kinaseLoc$Codon <- kinaseLoc$loc_actual_saltbridge
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$loc_actual_saltbridge, kinaseLoc$loc_actual_saltbridge)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_saltbridge <- merged
write.csv (merged_saltbridge, file='outputs/tcga/merged_saltbridge_all.csv')


# Critical lysine minus 2 = V[A]IK.
kinaseLoc$Codon <- kinaseLoc$actual_lysine - 2
kinaseLoc$AA <- str_sub(kinaseLoc$Protein_Seq, kinaseLoc$Codon, kinaseLoc$Codon)
merged <- merge(comboDATA, kinaseLoc, by = c('Gene', 'Codon', 'AA') )
merged$pastedref <- paste(merged$Gene, merged$Codon, merged$case_id, merged$AA, sep = "_" )
merged$DUP= !duplicated(merged$pastedref)
merged <- subset(merged, DUP == 'TRUE')
merged_Kminus2 <- merged
write.csv (merged_Kminus2, file='outputs/tcga/critical_lysine_minus2.csv')


# Bind tables together.
total <- rbind (merged_G2, merged_G3, merged_Kminus2, merged_HRD_H, merged_HRD_R, merged_HRD_D, merged_HRD_Dp5, merged_HRDm6, merged_DFG_D, merged_DFG_G, merged_APE_A, merged_APE_P, merged_APE_E, merged_saltbridge, merged_correctedG)
write.csv (total, file ='outputs/tcga/total_ALLTRANSCRIPTS.csv')

freqscreen <- table(total$Gene)
write.csv (freqscreen, file ='SOURCE_tcgafreqscreen_ALLTRANSCRIPTS.csv')         


