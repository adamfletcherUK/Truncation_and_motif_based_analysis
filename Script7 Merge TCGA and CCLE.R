#################################################################################################################
# Script 7: Merges the hotspot mutational data from CCLE and TCGA databases outputted in Scripts 5 and 6.       #
#                                                                                                               #
# Automatic merge function to combine the CCLE and TCGA hotspot mutational data from Scripts 5 and 6.           #
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
# Merge CCLE and TCGA outputs and orders based on total frequency of mutations per kinase.
########################

# Read CCLE and TCGA frequency files.
tcgafreq <- read.csv ( file = 'SOURCE_tcgafreqscreen_ALLTRANSCRIPTS.csv')
cclefreq <- read.csv ( file = 'SOURCE_cclefreqscreen_ALLTRANSCRIPTS.csv') #Error (File not found)


# Merge CCLE and TCGA data.
totalfreq = merge(tcgafreq, cclefreq, by = 'Var1', nomatch = 0)


# Remove unnecessary columns.
totalfreq$X.x = NULL
totalfreq$X.y = NULL


# Find the total frequency of mutations within each kinase and rank in descending order.
totalfreq$TotalFreq = with(totalfreq, Freq.x+Freq.y)
totalfreq <- totalfreq[with(totalfreq, order(-TotalFreq)),]


# Write file
write.csv (totalfreq, file ='outputs/totalfreq.csv')

