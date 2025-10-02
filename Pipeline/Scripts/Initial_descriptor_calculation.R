current_dir = getwd()
setwd("Scripts")

library(bio3d)
source("my_functions.R")


args = commandArgs(trailingOnly=TRUE)

if(length(args) !=4){
  print("Error: correct usage is")
  print("Rscript Initial_descriptor_calculation.R file1.pdb(fixed) file2.pdb(mutable) inter_file_1.csv(fixed) inter_file2.csv(mutable)")
  stop()
}

#read the input argument
name_prot_1 <- paste0('../',args[1])
name_prot_2 <- paste0('../',args[2])
name_int_1 <- paste0('../',args[3])
name_int_2 <- paste0('../',args[4])


####################################
## Surfaces calculation ############
####################################

print("###### Surfaces calculation #######")

tmp <- strsplit(name_prot_1,"/")[[1]]
data_path_1 <- paste(paste(tmp[1:(length(tmp)-1)],collapse = "/"),"/",sep="")
if(substring(name_prot_1,1,1) == "/"){
  data_path_1 <- paste("/",data_path_1,sep="")
}
data_name_1 <- tmp[length(tmp)]
output_path <- "../Files/Surf_tmp/"
output_name_1 <- substring(data_name_1,1,nchar(data_name_1)-4)

get.surf(data_path_1,data_name_1,1.4,output_path, output_name_1 )
surf_1 <- read.csv(paste(output_path,output_name_1,".csv",sep=""),stringsAsFactors = F)


tmp <- strsplit(name_prot_2,"/")[[1]]
data_path_2 <- paste(paste(tmp[1:(length(tmp)-1)],collapse = "/"),"/",sep="")
if(substring(name_prot_2,1,1) == "/"){
  data_path_2 <- paste("/",data_path_2,sep="")
}
data_name_2 <- tmp[length(tmp)]
output_path <- "../Files/Surf_tmp/"
output_name_2 <- substring(data_name_2,1,nchar(data_name_2)-4)

get.surf(data_path_2,data_name_2,1.4,output_path, output_name_2 )
surf_2 <- read.csv(paste(output_path,output_name_2,".csv",sep=""),stringsAsFactors = F)

###############################################################
###############################################################
## identification Interface Surface ###########################
###############################################################

Int_1 <- read.csv(name_int_1)
Int_1[is.na(Int_1$insert),]$insert <- ""
Int_2 <- read.csv(name_int_2)
Int_2[is.na(Int_2$insert),]$insert <- ""

res_id <- substring(surf_1[,1],1,3)
res_no <- substring(surf_1[,1],5,nchar(surf_1[,1])-1)
res_chain <- substring(surf_1[,1],nchar(surf_1[,1]),nchar(surf_1[,1]))

surf_1_BS <- surf_1[paste(res_id,res_no,res_chain) %in%  paste(Int_1$resid, Int_1$resno, Int_1$chain),]

res_id <- substring(surf_2[,1],1,3)
res_no <- substring(surf_2[,1],5,nchar(surf_2[,1])-1)
res_chain <- substring(surf_2[,1],nchar(surf_2[,1]),nchar(surf_2[,1]))

surf_2_BS <- surf_2[paste(res_id,res_no,res_chain) %in%  paste(Int_2$resid, Int_2$resno, Int_2$chain),]


write.csv(surf_1_BS,paste("../Files/Surf_BS_1/",paste(substring(data_name_1,1,nchar(data_name_1)-4),"_bs.csv",sep=""),sep=""),quote = F,row.names = F)
write.csv(surf_2_BS,paste("../Files/Surf_BS_2/",paste(substring(data_name_2,1,nchar(data_name_2)-4),"_bs.csv",sep=""),sep=""),quote = F,row.names = F)

print("The interface surfaces have been written")

################################################################################################
####### Zernike Calculations ###################################################################
################################################################################################
print("####### Zernike Calculations ###########")

system("ls ../Files/Surf_BS_1/ > auxiliary_files/names_1.txt")
system("ls ../Files/Surf_BS_2/ > auxiliary_files/names_2.txt")

system("python3 Z2D_verso_1.py auxiliary_files/names_1.txt")
system("python3 Z2D_verso_-1.py auxiliary_files/names_2.txt")


################################################################################################
####### calculation of Zernike distances #######################################################
################################################################################################

name_file_zernike_1 <- list.files("../Files/ZernDescr/", pattern = substring(data_name_1,1,nchar(data_name_1)-4))
name_file_zernike_2 <- list.files("../Files/ZernDescr/", pattern = substring(data_name_2,1,nchar(data_name_2)-4))


ZD_1 <- read.csv(paste("../Files/ZernDescr/",name_file_zernike_1,sep=""),header=F)
ZD_2 <- read.csv(paste("../Files/ZernDescr/",name_file_zernike_2,sep=""),header=F)

Zern_dist <- as.numeric(dist(rbind(ZD_1[,1],ZD_2[,1])))

##############################################################################################
##### Electrostatics Calculation #############################################################
##############################################################################################

print("##### Electrostatics Calculation #####")

system(paste("pdb2pqr", name_prot_1,
             paste("../Files/PQR/",substring(data_name_1,1,nchar(data_name_1)-4),".pqr",sep=""),
             "--ff=CHARMM","--keep-chain", sep = " "))



system(paste("pdb2pqr", name_prot_2,
             paste("../Files/PQR/",substring(data_name_2,1,nchar(data_name_2)-4),".pqr",sep=""),
             "--ff=CHARMM","--keep-chain", sep = " "))


system("rm ../Files/PQR/*.log")


name_file_elec_1 <- list.files("../Files/PQR/", pattern = substring(data_name_1,1,nchar(data_name_1)-4))
name_file_elec_2 <- list.files("../Files/PQR/", pattern = substring(data_name_2,1,nchar(data_name_2)-4))

PQR_1 <- read.pdb(paste("../Files/PQR/",name_file_elec_1,sep=""))
PQR_1 <- PQR_1$atom
PQR_2 <- read.pdb(paste("../Files/PQR/",name_file_elec_2,sep=""))
PQR_2 <- PQR_2$atom

residui_PQR_1<- unique(paste(PQR_1$chain,PQR_1$resno,PQR_1$resid)) 

PQR_1_CG <- c()
for(k in 1:length(residui_PQR_1)){
  aus <- strsplit(residui_PQR_1[k]," ")[[1]]
  single_res <- PQR_1[PQR_1$chain %in% aus[1] & PQR_1$resno %in% aus[2] & PQR_1$resid %in% aus[3],]
  
  single_res_backbone <- single_res[single_res$elety %in% c("N","CA","C","O"),]
  riga_bb <- single_res_backbone[2,]
  riga_bb[,9:11] <- apply(single_res_backbone[,9:11],2,mean)
  riga_bb[,12] <- sum(single_res_backbone[,12])
  
  single_res_sidechain <- single_res[! single_res$elety %in% c("N","CA","C","O"),]
  riga_sc <- single_res_sidechain[1,]
  riga_sc$elety <- "B"
  riga_sc[,9:11] <-apply(single_res_sidechain[substring(single_res_sidechain$elety,1,1) != "H",9:11], 2,mean) 
  if(aus[3] =="GLY"){riga_sc[,9:11] <- apply(single_res_sidechain[,9:11], 2,mean)}
  riga_sc[,12] <- sum(single_res_sidechain[,12])
  
  
  PQR_1_CG <- rbind(PQR_1_CG,riga_bb,riga_sc)
  
}

residui_PQR_2<- unique(paste(PQR_2$chain,PQR_2$resno,PQR_2$resid)) 

PQR_2_CG <- c()
for(k in 1:length(residui_PQR_2)){
  aus <- strsplit(residui_PQR_2[k]," ")[[1]]
  single_res <- PQR_2[PQR_2$chain %in% aus[1] & PQR_2$resno %in% aus[2] & PQR_2$resid %in% aus[3],]
  
  single_res_backbone <- single_res[single_res$elety %in% c("N","CA","C","O"),]
  riga_bb <- single_res_backbone[2,]
  riga_bb[,9:11] <- apply(single_res_backbone[,9:11],2,mean)
  riga_bb[,12] <- sum(single_res_backbone[,12])
  
  single_res_sidechain <- single_res[! single_res$elety %in% c("N","CA","C","O"),]
  riga_sc <- single_res_sidechain[1,]
  riga_sc$elety <- "B"
  riga_sc[,9:11] <-apply(single_res_sidechain[substring(single_res_sidechain$elety,1,1) != "H",9:11], 2,mean) 
  if(aus[3] =="GLY"){riga_sc[,9:11] <- apply(single_res_sidechain[,9:11], 2,mean)}
  riga_sc[,12] <- sum(single_res_sidechain[,12])
  
  
  PQR_2_CG <- rbind(PQR_2_CG,riga_bb,riga_sc)
  
}


Mat_dist <- dist.xyz(PQR_1_CG[,9:11],PQR_2_CG[,9:11])

PQR_1_CG_BS <- PQR_1_CG[apply(Mat_dist,1,min) < 12,]

PQR_2_CG_BS <- PQR_2_CG[apply(Mat_dist,2,min) < 12,]

Mat_dist_BS <- dist.xyz(PQR_1_CG_BS[,9:11], PQR_2_CG_BS[,9:11])

En_tot <- c()
const = 1./(6.24151**2 * 8.85418782*4*pi)*0.000239006*6.022*10**(10 + 23 + 12 - 36)
for(i in 1:nrow(PQR_1_CG_BS)){
  for(j in 1:nrow(PQR_2_CG_BS)){
  
    En <-  const* ((PQR_1_CG_BS[i,]$o * PQR_2_CG_BS[j,]$o) / (Mat_dist_BS[i,j])) 
    En_tot <- rbind(En_tot,c(i,j,En))
    
  }
}

En_INT <- sum(En_tot[,3])

#####################################################################
####### hydrophobicity calculation ##################################
#####################################################################
print("####### hydrophobicity calculation ########")

hyd_scale <- c(0.47, 1.60, 2.03, 3.20, 0.59, 1.26, 3.32,
               0.80, 2.29, 0.00, 0.79, 2.76, 0.38, 0.24,
               0.55, 1.46, 0.66, 0.65, 0.83, 0.46)
names(hyd_scale) <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                      "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                      "PRO", "SER", "THR", "TRP", "TYR", "VAL")

inv_res <- table(surf_1_BS$Res)
inv_res <- inv_res[inv_res != 0]
names(inv_res) <- substring(names(inv_res),1,3)

inv_res_hyd <- hyd_scale[names(inv_res)]
hydro_patch_1 <- sum(inv_res_hyd * inv_res)
hydro_patch_1 <- round(hydro_patch_1 / nrow(surf_1_BS),3)

inv_res <- table(surf_2_BS$Res)
inv_res <- inv_res[inv_res != 0]
names(inv_res) <- substring(names(inv_res),1,3)

inv_res_hyd <- hyd_scale[names(inv_res)]
hydro_patch_2 <- sum(inv_res_hyd * inv_res)
hydro_patch_2 <- round(hydro_patch_2 / nrow(surf_2_BS),3)


hydro_similarity <- abs(hydro_patch_1 - hydro_patch_2)

#####################################################################
#####################################################################

diff_descriptors <- rbind(round(Zern_dist,3),round(En_INT,3), round(hydro_similarity,3))
write.csv(diff_descriptors, "../Files/Initial_descriptors.csv", quote=F, row.names = F)

setwd(current_dir)
