#load packages
library(readxl)
library(dplyr)
library(tidyverse)
library(matrixStats)

#import data file with total spacer counts for each sample (2 conditions, 2 technical replicates, 2 biological replicates = 8 samples total)
combined <- read_excel("C:/Users/samli/Downloads/CRISPRi_Library_Screen_Analysis/CRISPRi_Library_Screen_Combined_Python_Output_Files.xlsx", 
                       sheet = "combined")

#ensure spacer count columns for each sample are read as integers
combined[, c(2:9)] <- sapply(combined[, c(2:9)], as.integer)

#calculate and display the total spacer count for each sample
mapply(sum,combined[,c(2:9)])

#calculate and display the ratio of the spacer count sums for each untreated/Sal4-treated pair
sum(combined[,c(2)])/ sum(combined[,c(3)]); sum(combined[,c(4)])/ sum(combined[,c(5)]);sum(combined[,c(6)])/ sum(combined[,c(7)]);sum(combined[,c(8)])/ sum(combined[,c(9)])

#divide the untreated column reads by the ratio of the untreated/Sal4-treated sums and insert those values in a new column before the Sal4-treated reads for each sample in a new table
normalized <-add_column(combined, NU = combined[,c(2)]/(sum(combined[,c(2)])/ sum(combined[,c(3)])), .before = 3)
normalized <-add_column(normalized, NU2 = normalized[,c(5)]/(sum(normalized[,c(5)])/ sum(normalized[,c(6)])), .before = 6)
normalized <-add_column(normalized, NU3 = normalized[,c(8)]/(sum(normalized[,c(8)])/ sum(normalized[,c(9)])), .before = 9)
normalized <-add_column(normalized, NU4 = normalized[,c(11)]/(sum(normalized[,c(11)])/ sum(normalized[,c(12)])), .before = 12)

#make a new data table with only the normalized untreated and Sal4-treated reads
normalized2 <- normalized[,c(1,3,4,6,7,9,10,12,13)]

#convert data frames within the data frame into columns
normalized2 <- do.call(data.frame, normalized2)

#rename columns
colnames(normalized2)[2] <- "E1_P2_NU_R1"
colnames(normalized2)[3] <- "E1_P2_Sal4_R1"
colnames(normalized2)[4] <- "E1_P2_NU_R2"
colnames(normalized2)[5] <- "E1_P2_Sal4_R2"
colnames(normalized2)[6] <- "E2_P2_NU_R1"
colnames(normalized2)[7] <- "E2_P2_Sal4_R1"
colnames(normalized2)[8] <- "E2_P2_NU_R2"
colnames(normalized2)[9] <- "E2_P2_Sal4_R2"

#round all values in numeric columns to the nearest integer
rounded <- normalized2 %>% mutate_if(is.numeric, round, 0)

#make a new table with the untreated and Sal4-treated columns for each paired sample
sp1 <- rounded[,c(1,2,3)]
sp2 <- rounded[,c(1,4,5)]
sp3 <- rounded[,c(1,6,7)]
sp4 <- rounded[,c(1,8,9)]

#exclude all spacers with less than 100 total read counts
sp1f100 <- filter(sp1, sp1[,c(2)]>=100 & sp1[,c(3)]>=100)
sp2f100 <- filter(sp2, sp2[,c(2)]>=100 & sp2[,c(3)]>=100)
sp3f100 <- filter(sp3, sp3[,c(2)]>=100 & sp3[,c(3)]>=100)
sp4f100 <- filter(sp4, sp4[,c(2)]>=100 & sp4[,c(3)]>=100)

#calculate the ratio of Sal4-treated/untreated reads for each spacer and add this to a new column
sp1f100r <-add_column(sp1f100, Ratio1 = sp1f100[,c(3)]/sp1f100[,c(2)], .before = 4)
sp2f100r <-add_column(sp2f100, Ratio2 = sp2f100[,c(3)]/sp2f100[,c(2)], .before = 4)
sp3f100r <-add_column(sp3f100, Ratio3 = sp3f100[,c(3)]/sp3f100[,c(2)], .before = 4)
sp4f100r <-add_column(sp4f100, Ratio4 = sp4f100[,c(3)]/sp4f100[,c(2)], .before = 4)

### ENRICHED GENES ANALYSIS ###
#exclude all spacers with a ratio of less than two
sp1r2 <- filter(sp1f100r, sp1f100r[,c(4)]>=2)
sp2r2 <- filter(sp2f100r, sp2f100r[,c(4)]>=2)
sp3r2 <- filter(sp3f100r, sp3f100r[,c(4)]>=2)
sp4r2 <- filter(sp4f100r, sp4f100r[,c(4)]>=2)

#calculate the log2 value of each spacer ratio and add this to a new column
sp1r2log2 <-add_column(sp1r2, Log2Ratio1 = log(sp1r2$Ratio1,2), .before = 5)
sp2r2log2 <-add_column(sp2r2, Log2Ratio2 = log(sp2r2$Ratio2,2), .before = 5)
sp3r2log2 <-add_column(sp3r2, Log2Ratio3 = log(sp3r2$Ratio3,2), .before = 5)
sp4r2log2 <-add_column(sp4r2, Log2Ratio4 = log(sp4r2$Ratio4,2), .before = 5)

#intersect filtered sample read count lists by matching First 28 nt (unique) values
intersect12e <- inner_join(sp1r2log2,sp2r2log2, by = "First_28_nt_.unique.")
intersect34e <- inner_join(sp3r2log2,sp4r2log2, by = "First_28_nt_.unique.")
intersecte <- inner_join(intersect12e,intersect34e, by = "First_28_nt_.unique.")

###

### REDUCED GENES ANALYSIS ###
#exclude all spacers with a ratio of greater than 0.5
sp1r0.5 <- filter(sp1f100r, sp1f100r[,c(4)]<=0.5)
sp2r0.5 <- filter(sp2f100r, sp2f100r[,c(4)]<=0.5)
sp3r0.5 <- filter(sp3f100r, sp3f100r[,c(4)]<=0.5)
sp4r0.5 <- filter(sp4f100r, sp4f100r[,c(4)]<=0.5)

#calculate the log2 value of each spacer ratio and add this to a new column
sp1r0.5log2 <-add_column(sp1r0.5, Log2_Ratio1 = log(sp1r0.5$Ratio1,2), .before = 5)
sp2r0.5log2 <-add_column(sp2r0.5, Log2_Ratio2 = log(sp2r0.5$Ratio2,2), .before = 5)
sp3r0.5log2 <-add_column(sp3r0.5, Log2_Ratio3 = log(sp3r0.5$Ratio3,2), .before = 5)
sp4r0.5log2 <-add_column(sp4r0.5, Log2_Ratio4 = log(sp4r0.5$Ratio4,2), .before = 5)

#intersect filtered sample read count lists by matching First 28 nt (unique) values
intersect12r <- inner_join(sp1r0.5log2,sp2r0.5log2, by = "First_28_nt_.unique.")
intersect34r <- inner_join(sp3r0.5log2,sp4r0.5log2, by = "First_28_nt_.unique.")
intersectr <- inner_join(intersect12r,intersect34r, by = "First_28_nt_.unique.")

###

#rename column 1 on combined enriched/reduced lists
colnames(intersecte)[1] <- "First_28_nt_(unique)"
colnames(intersectr)[1] <- "First_28_nt_(unique)"

#calculate the average Log2 Ratio value for each spacer across the four samples and add values to a new column
intersecte2 <-add_column(intersecte, Average_Log2_Ratio = rowMeans(intersecte[ , c(5, 9, 13, 17)]), .before = 18)
intersectr2 <-add_column(intersectr, Average_Log2_Ratio = rowMeans(intersectr[ , c(5, 9, 13, 17)]), .before = 18)

#write standard deviation calculation function
rowVars <- function(x, na.rm=F) {rowSums((x - rowMeans(x, na.rm=na.rm))^2, na.rm=na.rm) / (ncol(x) - 1)}

#calculate the standard deviation of the Log2 Ratio values for each spacer across the four samples and add to a new column
intersecte2 <-add_column(intersecte2, Standard_Deviation = sqrt(rowVars(intersecte2[ , c(5, 9, 13, 17)])), .before = 19)
intersectr2 <-add_column(intersectr2, Standard_Deviation = sqrt(rowVars(intersectr2[ , c(5, 9, 13, 17)])), .before = 19)

#import spacer information and assignment list excel file
spacers <- read_excel("C:/Users/samli/Downloads/CRISPRi_Library_Screen_Analysis/CRISPRi_Library_Spacer_Assignment_List.xlsx")

#intersect spacer info and combined enriched/reduced read count lists by matching First 28 nt (unique) values
enriched_spacers <- inner_join(spacers,intersecte2, by = "First_28_nt_(unique)")
reduced_spacers <- inner_join(spacers,intersectr2, by = "First_28_nt_(unique)")

#import STm 14028s genome annotation excel file
genome <- read_excel("C:/Users/samli/Downloads/CRISPRi_Library_Screen_Analysis/STm_14028s_Full_Annotated_Genome_Information.xlsx")

#intersect annotated enriched/reduced spacer and genome info lists by Old_Locus_Tag values
enriched_list <- right_join(genome,enriched_spacers, by = "Old_Locus_Tag")
reduced_list <- right_join(genome,reduced_spacers, by = "Old_Locus_Tag")

#remove unneeded/duplicated columns
final_enriched <- enriched_list[ -c(1, 3, 8, 10, 14, 15) ]
final_reduced <- reduced_list[ -c(1, 3, 8, 10, 14, 15) ]

#reorder columns
final_enriched <- final_enriched[, c(13,14,18,19,16,17,15,1,6,20,37,38,7,9,3,4,5,8,2,10,11,12,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]
final_reduced <- final_reduced[, c(13,14,18,19,16,17,15,1,6,20,37,38,7,9,3,4,5,8,2,10,11,12,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)]

#rename columns
colnames(final_enriched)[4] <- "Spacer_ID"
colnames(final_enriched)[5] <- "Spacer_Coordinate"
colnames(final_enriched)[6] <- "Location"
colnames(final_enriched)[7] <- "Strand"
colnames(final_enriched)[8] <- "New_Locus_Tag"
colnames(final_enriched)[10] <- "Gene_Name"
colnames(final_enriched)[13] <- "Predicted_Function"

colnames(final_reduced)[4] <- "Spacer_ID"
colnames(final_reduced)[5] <- "Spacer_Coordinate"
colnames(final_reduced)[6] <- "Location"
colnames(final_reduced)[7] <- "Strand"
colnames(final_reduced)[8] <- "New_Locus_Tag"
colnames(final_reduced)[10] <- "Gene_Name"
colnames(final_reduced)[13] <- "Predicted_Function"

#make a new data frame with rows that have multiple Old_Locus_Tag entries in the enriched/reduced lists
enriched_mult <- final_enriched %>% 
  group_by(Old_Locus_Tag) %>% 
  filter(n()>1)

reduced_mult <- final_reduced %>% 
  group_by(Old_Locus_Tag) %>% 
  filter(n()>1)

#exclude all instances of Promoter from the Old_Locus_Tag column in the duplicated enriched/reduced lists
enriched_mult <- enriched_mult %>%
  filter(Old_Locus_Tag!='Promoter')

reduced_mult <- reduced_mult %>%
  filter(Old_Locus_Tag!='Promoter')

#export select data frames into Excel file
openxlsx::write.xlsx(list(final_enriched = final_enriched,
                enriched_mult = enriched_mult,
                final_reduced = final_reduced,
                reduced_mult = reduced_mult),
           'C:/Users/samli/Downloads/CRISPRi_Library_Screen_Analysis/CRISPRi_Library_Screen_Analysis.xlsx')
