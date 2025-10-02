library(bio3d)

args = commandArgs(trailingOnly=TRUE)

if(length(args) !=3){
  print("Error: correct usage is")
  print("Rscript Identification_interfaces.R file1.pdb file2.pdb Interaction_thresold")
  stop()
}

#read the input argument
name_prot_1 <- args[1]
name_prot_2 <- args[2]
THR_INT <- as.numeric(args[3])

# read the pdb files
prot_1 <- read.pdb(name_prot_1)
prot_1 <- prot_1$atom
prot_1[is.na(prot_1$insert),]$insert <- ""

prot_2 <- read.pdb(name_prot_2)
prot_2 <- prot_2$atom
prot_2[is.na(prot_2$insert),]$insert <- ""

# compute the distance matrix
mat_dist <- dist.xyz(prot_1[,9:11], prot_2[,9:11])

# define the proteins binding sites
prot_1_BS <- prot_1[apply(mat_dist,1,min) < THR_INT,]
prot_2_BS <- prot_2[apply(mat_dist,2,min) < THR_INT,]

prot_1_BS <- unique(prot_1_BS[,5:8])
prot_2_BS <- unique(prot_2_BS[,5:8])

# write the output
write.csv(prot_1_BS,"Files/interface_prot_1.csv", quote = F, row.names = F)
write.csv(prot_2_BS,"Files/interface_prot_2.csv", quote = F, row.names = F)

print("The interface files have been written")
