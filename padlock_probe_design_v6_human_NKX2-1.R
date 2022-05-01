# load the necessary libraries
library("plyr") # if not installed, run first the: install.packages("plyr")
library("dplyr") # if not installed, run first the: install.packages("dplyr")
library("spgs") # if not installed, run first the: install.packages("spgs")
library("rmelting")
# if not installed, run first the: 
#if (!require("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")

#BiocManager::install("rmelting")

# set the working directory
setwd("C:/Users/alex/Desktop/")
#====================================================================================================================
# The present script uses specific Taqman probe sequences, which have been designed with IDT PrimerQuest Tool, to create padlock probes.
# During the analysis, we should set some parameters, according to the results of the melting temperature (Tm) calculations.
# The parameters to be set are fount in lines: 39, 41, 43, 48-51, 54, 57, 60, 112, 115, 179. The result of the pipeline is 
# a *.csv file with all the information about the padlock probe and its detection oligo. 
# Considering that, we have to set some parameters, during the analysis, further automation would be problematic.
#====================================================================================================================
{
#Buffer sequence
A <- "TCCTCTATGATTACTGAC"
# Backbone-1 (B1)
B1 <- "TGCGTCTATTTAGTGGAGCC"
# Backbone-2(B2) for highly abundant genes
B2 <- "GTATATCTTGCGTGACGCTG"
anchors <- data.frame(matrix(data = 0, nrow = 2, ncol=2))
anchors[1,1] <- "B1"
anchors[1,2] <- B1
anchors[2,1] <- "B2"
anchors[2,2] <- B2

# Second buffer sequence
D <- "ATGCCTATCTTCTTT"
}
#=========================================
# Set the information about the targeted gene
target <- "NKX2_1"
# Set the NCBI number of the targeted mRNA (not gene)
NCBI_no <- "NM_003317.4"
# Set the species of the targeted gene (h: human, m: mouse, r: rat)
species <- "h"


# Provide the Taqman probe sequences for probe design. It is important to use the proper orientation of the sequence. BLAST to transcriptome will 
# help checking the orientation and obtain the proper sequence.
Taq1 <- "CACACTCTGCCAGCAAAGAGGACTCGCTTGTAAATACCAG" #1719-1758
Taq2 <- "CTGACATGGCTCTGGACTCTAAAGACCAAACTTCACTCTG" #1677-1716
Taq3 <- "CATCCAATCTCAAGGAATCTTTAAGCAGAGAAGGGCATAA" #1555-1594
Taq4 <- "GAAAGAGTCTTCAACATAACCCACTTGTCACTGACACAAAG" #1802-1842

# Set the sequence for probe design
Probe <- 1

# Set the sequence for the anchor primer (1: low and medium abundance, or 2: high abundance)
B_n <- 1

# set the fluorophore of the detection oligo, between the [FAM], [Cy3], [TxRed], [Cy5], [AF750]
fluor <- "[Cy3]"
#-------------------------------
# backbone construction
{B <- anchors[B_n,2] 
backbone <- paste(A, B, D)
# Check that the backbone length is 55nt (53 plus 2 "spaces")
nchar(backbone)


#==================================================
probe_list <- list(Taq1, Taq2, Taq3, Taq4)
names(probe_list) <- c("Probe1", "Probe2", "Probe3", "Probe4")
Taq<- probe_list[[Probe]]

# Design of padlock arms
arm_list <- data.frame(matrix(nrow = 20, ncol =26))
arm_list<- setNames(arm_list, c("Seq", "Tm"))
for (i in 0:9){
        a <- substr(Taq, 1, (nchar(Taq)-i))
        t1 <- melting(sequence = a, nucleic.acid.conc = 5*10^(-8), 
                      hybridisation.type = "dnadna", Na.conc = 1.25*10^(-3), Mg.conc=10*10^(-3), formamide.conc=20, Tris.conc=20*10^(-3), K.conc=75*10^(-3),correction.formamide="lincorr")
        arm_list[(2*i+1),1] <- t1[["Environment"]][["Sequence"]]
        arm_list[(2*i+1),2] <-t1[["Results"]][["Melting temperature (C)"]]
        b <- substr(Taq, (i+1), nchar(Taq))
        t2 <- melting(sequence = b, nucleic.acid.conc = 5*10^(-8), 
                      hybridisation.type = "dnadna", Na.conc = 1.25*10^(-3), Mg.conc=10*10^(-3), formamide.conc=20, Tris.conc=20*10^(-3), K.conc=75*10^(-3),correction.formamide="lincorr")
        arm_list[(2*i+2),1] <- t2[["Environment"]][["Sequence"]]
        arm_list[(2*i+2),2] <-t2[["Results"]][["Melting temperature (C)"]]
        rm(a, b, t1, t2)
}

# The script calculates the Tm values of the arms, for 6 nucleotides around the middle of the targeted sequence. 
for (i in 1:12)
for (j in 1:20){
m <- as.integer(nchar(arm_list[j,1])/2)
a5 <- substr(arm_list[j,1], 1, m+i-6)

t_a5 <- melting(sequence = a5, nucleic.acid.conc = 5*10^(-8), 
                hybridisation.type = "dnadna", Na.conc = 1.25*10^(-3), Mg.conc=10*10^(-3), formamide.conc=20, Tris.conc=20*10^(-3), K.conc=75*10^(-3),correction.formamide="lincorr")
        arm_list[j,(2*i)+1] <-t_a5[["Results"]][["Melting temperature (C)"]]
        names(arm_list)[(2*i)+1] <- paste0("5arm-nt",i)
        
a3 <- substr(arm_list[j,1], m+i-5, nchar(arm_list[j,1]))
t_a3 <- melting(sequence = a3, nucleic.acid.conc = 5*10^(-8), 
                hybridisation.type = "dnadna", Na.conc = 1.25*10^(-3), Mg.conc=10*10^(-3), formamide.conc=20, Tris.conc=20*10^(-3), K.conc=75*10^(-3),correction.formamide="lincorr")
        arm_list[j,(2*i)+2] <-t_a3[["Results"]][["Melting temperature (C)"]]
        names(arm_list)[(2*i)+2] <- paste0("3arm-nt",i)
}
rm(m, a5, a3, t_a5, t_a3, i, j)
}
# The Tm of the targeted sequence should be 60+/-3oC. 
# Set "M" between 1 and 20 based on the Tm results of the "arm_list". The average Tm of 12 randomly selected is 60oC but we have also working probes with Tm ~55oC
M <- 7

#set the optimal cutting position in order to have equal Tm values for the two branches, close to 45oC.
k <- 4

{m <- as.integer(nchar(arm_list[M,1])/2)
a5 <- substr(arm_list[M,1], 1, m+k-6)
a3 <- substr(arm_list[M,1], m+k-5, nchar(arm_list[M,1]))
c_a3 <-reverseComplement(a3, case="upper")
t_c_a3 <- melting(sequence = c_a3, nucleic.acid.conc = 5*10^(-8), 
                  hybridisation.type = "dnadna", Na.conc = 1.25*10^(-3), Mg.conc=10*10^(-3), formamide.conc=20, Tris.conc=20*10^(-3), K.conc=75*10^(-3),correction.formamide="lincorr")

c_a5 <-reverseComplement(a5, case="upper")
t_c_a5 <- melting(sequence = c_a5, nucleic.acid.conc = 5*10^(-8), 
                  hybridisation.type = "dnadna", Na.conc = 1.25*10^(-3), Mg.conc=10*10^(-3), formamide.conc=20, Tris.conc=20*10^(-3), K.conc=75*10^(-3),correction.formamide="lincorr")

#===========================================
#===========================================
padlock <- paste(c_a5, backbone, c_a3)

padlock1 <- paste0(c_a5, A, B, D, c_a3)
padlock_info <- melting(sequence = padlock1, nucleic.acid.conc = 5*10^(-8), 
                        hybridisation.type = "dnadna", Na.conc = 1.25*10^(-3), Mg.conc=10*10^(-3), formamide.conc=20, Tris.conc=20*10^(-3), K.conc=75*10^(-3),correction.formamide="lincorr")


#===========================================================

# Set the padlock specific barcode. It is the same sequence as the recognized by the arms
C<- c(arm_list[M,1])

# measure the Tm of the candidate barcode and choose the proper length for Tm ~36.2oC 
# The length of the default barcode is 23nt. 
barcode_list <- data.frame(matrix(nrow = 42, ncol =2))
barcode_list<- setNames(barcode_list, c("Seq", "Tm"))

a <- C
b <- C
t1 <- melting(sequence = a, nucleic.acid.conc = 5*10^(-8), 
              hybridisation.type = "dnadna", Na.conc = 39*10^(-3), formamide.conc=20,correction.formamide="lincorr")
barcode_list[1,1] <- t1[["Environment"]][["Sequence"]]
barcode_list[1,2] <-t1[["Results"]][["Melting temperature (C)"]]
barcode_list[2,1] <- t1[["Environment"]][["Sequence"]]
barcode_list[2,2] <-t1[["Results"]][["Melting temperature (C)"]]


for (i in 1:20){
        a <- substr(b, 1, (nchar(b)-1))
        t1 <- melting(sequence = a, nucleic.acid.conc = 5*10^(-8), 
                      hybridisation.type = "dnadna", Na.conc = 39*10^(-3), formamide.conc=20,correction.formamide="lincorr")
        barcode_list[(2*i+1),1] <- t1[["Environment"]][["Sequence"]]
        barcode_list[(2*i+1),2] <-t1[["Results"]][["Melting temperature (C)"]]
        
        b <- substr(a, 2, (nchar(a)))
        t2 <- melting(sequence = b, nucleic.acid.conc = 5*10^(-8), 
                      hybridisation.type = "dnadna", Na.conc = 39*10^(-3), formamide.conc=20,correction.formamide="lincorr")
        barcode_list[(2*i+2),1] <- t2[["Environment"]][["Sequence"]]
        barcode_list[(2*i+2),2] <-t2[["Results"]][["Melting temperature (C)"]]
        rm(t1, t2)
}

# select the optimal barcode length according to: 
opt <- barcode_list[barcode_list$Tm >=35.7 & barcode_list$Tm <=36.7,]
print(opt)
}
# if no value is returned, check the "backbone_Tm_list" and
# select the rowname with Tm close to the 36.2oC  
L=22
{
barcode <- barcode_list[L,]

#=====================================================
# prepare a direct detection oligo, that recognises the gene specific backbone

{direct_oligo <- paste(reverseComplement(barcode$Seq, case="upper"), fluor)

#====================================================
        
        
#====================================================
# Create a summary file
summary <- data.frame(matrix(nrow = 14, ncol =2))

summary[[1,1]] <-"mRNA target"
summary[[1,2]] <- target

summary[[2,1]] <-"Species"
summary[[2,2]] <- species

summary[[3,1]] <-"NCBI ID"
summary[[3,2]] <- NCBI_no


summary[[4,1]] <- "Anchor Number"
summary[[4,2]] <- paste0("A",B_n)

summary[[5,1]] <- "targeted sequence"
summary[[5,2]] <- c(arm_list[M,1])

summary[[6,1]] <- "padlock Sequence"
summary[[6,2]] <- padlock

summary[[7,1]] <- "padlock Tm (oC)"
summary[[7,2]] <- padlock_info[["Results"]][["Melting temperature (C)"]]

summary[[8,1]] <- "5'-Arm Tm"
summary[[8,2]] <- t_c_a5[["Results"]][["Melting temperature (C)"]]

summary[[9,1]] <- "3'-Arm Tm"
summary[[9,2]] <- t_c_a3[["Results"]][["Melting temperature (C)"]]

summary[[10,1]] <- "targeted sequence Tm (oC)"
summary[[10,2]] <- arm_list[M,2]

summary[[11,1]] <- "Padlock Name"
summary[[11,2]] <- paste0(target, "_", species, "_pr", Probe, "_A", B_n, "_SplR_D")

summary[[12,1]] <- "Detection oligo Name"
summary[[12,2]] <- paste0(target, "_", species, "_oligo", Probe, "_", fluor)

summary[[13,1]] <- "Detection oligo Sequence"
summary[[13,2]] <- direct_oligo

summary[[14,1]] <- "Detection oligo Tm (oC)"
summary[[14,2]] <- barcode$Tm



# export the summary (by default in the "documents" directory)
name = paste0(target,"_", NCBI_no,"_padlock_", Probe, "_summary.csv")
write.csv(summary, file=name)
}
}
