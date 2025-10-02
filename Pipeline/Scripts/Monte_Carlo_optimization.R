#####################################################################
#### manage the user input ##########################################
#####################################################################

args = commandArgs(trailingOnly=TRUE)

if(length(args) !=6){
  print("Error: correct usage is")
  print("Rscript MonteCarlo_Optimization.R  mutable.pdb inter_file_mutable.csv surface_BS_fixed.csv Zern_file_fixed.csv fixed.pqr initial_descriptors.csv")
  stop()
}


#read the input argument
name_mut <- paste0('../',args[1])
name_int_mut <- paste0('../',args[2])
name_surf_BS_fixed <- paste0('../',args[3])
name_Zern_fixed <- paste0('../',args[4])
name_PQR_fixed <- paste0('../',args[5])
name_descr_ini <- paste0('../',args[6])



#########################################
#########################################
#### libraries and parameters ###########
#########################################
current_dir = getwd()
setwd("Scripts")

library(bio3d)
library(seqinr)
source("my_functions.R")


##################################
# weight of the cost function ####
A <- 10
B <- 0.04
C <- 32.1
D <- 0.5

#################################
# parameter of the simulation ###
#################################
num_step <- 1500
Beta_range <- seq(0.5,2,0.5)


################################
current_dir <- getwd()

#######################################################################################
### initialization of all the variables necessary to the simulation ###################
#######################################################################################

possible_Res <- c("GLY", "ALA", "VAL", "LEU", "ILE", "MET", "SER", "PRO", "THR", "CYS",
                   "ASN", "GLN", "PHE", "TYR", "TRP", "LYS", "HIS", "ARG", "ASP", "GLU")

hyd_scale <- c(0.47, 1.60, 2.03, 3.20, 0.59, 1.26, 3.32,
               0.80, 2.29, 0.00, 0.79, 2.76, 0.38, 0.24,
               0.55, 1.46, 0.66, 0.65, 0.83, 0.46)
names(hyd_scale) <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                      "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                      "PRO", "SER", "THR", "TRP", "TYR", "VAL")


# initial sequence
struct_ini <- read.pdb(name_mut)
struct_ini <- struct_ini$atom
seq_ini <- my.get.seq(struct_ini)

# initial descriptors
descriptors_ini <- t(read.csv(name_descr_ini))
# interface where to find the mutable residues
Interface_ini <- read.csv(name_int_mut)

# variables to be varied during the cycle
seq_cycle <- seq_ini
dist_zern_cycle <- descriptors_ini[1,1]
EN_tot_cycle <- descriptors_ini[1,2]
dist_hydro_cycle <- descriptors_ini[1,3]
Interface_cycle <- Interface_ini


df_tot <- c("0","-","-","-","-",dist_zern_cycle,"-", EN_tot_cycle,"-",dist_hydro_cycle,"-","-","-")
names(df_tot) <- c("NoStep", "Beta", "Substituted Residue", "ResNo", "Inserted Residue",
                      "Selected Zernike", "Proposed Zernike", "Selected Elec", "Proposed Elec", 
                      "Selected Hydro", "Proposed Hydro", "Mutations number", "Acceptance")
for(Beta in 1:length(Beta_range)){
  for(step in 1:num_step){
    
    actual_step <- step + num_step*Beta - num_step
    print(actual_step)
    ####################################
    ## Mutation modeling ###############
    ####################################
    
    pos_res_mut <- sample(nrow(Interface_cycle),1)
    poss_Res_random <- possible_Res[! possible_Res %in% Interface_cycle[pos_res_mut,1]]
    res_to_ins <- sample(poss_Res_random,1)
    tmp <- seq_cycle[names(seq_cycle) %in% Interface_cycle[pos_res_mut,3]]
    
    seq_valutanda <- seq_cycle
    seq_valutanda[names(seq_cycle) %in% Interface_cycle[pos_res_mut,3]] <- res_to_ins
    div_seq <- unname(which(seq_valutanda != seq_ini))
    
    Cterm = c()
    Nterm = c()
    inside = c()
    
    Cterm = div_seq[which(div_seq == length(seq_cycle))]
    Nterm = div_seq[which(div_seq == 1)]
    inside = div_seq[which((div_seq != 1) & (div_seq != length(seq_cycle)))]
    
    if(length(Nterm != 0)){
      Nterm = c(Nterm, Nterm + 1)
    }
    if(length(Cterm != 0)){
      Cterm = c(Cterm, Cterm - 1)
    }
    if(length(inside != 0)){
      inside = c(inside, inside + 1, inside - 1)
    }
    
    div_seq = c(Nterm,inside,Cterm)
    
    seq_valutanda <- paste(substring(seq_valutanda, 1,1), tolower(substring(seq_valutanda,2,3)), sep="")
    seq_valutanda <- tolower(a(seq_valutanda))
    seq_valutanda[div_seq] <- toupper(seq_valutanda[div_seq])
    seq_valutanda <- paste(seq_valutanda, collapse = "")
    
    write(seq_valutanda, paste("../Files/sequences/seq_step_", actual_step,".txt", sep=""))
    
    old_res <- Interface_cycle[pos_res_mut,1]
    new_res <- res_to_ins
    num_res <- Interface_cycle[pos_res_mut,3]
    
    #### structure prediction via Scwrl4 ####
    print('+++ structure prediction +++')
    setwd("path_to_Scwrl_dir") ###############remember to check the input and output files path (Scwrl location)
    
    input_file <- paste(current_dir,name_mut,sep="/")
    
    mutated_seq <- list.files(paste(current_dir,"../Files/sequences/",sep = "/"), pattern = paste("seq_step_",actual_step,".txt",sep=""))
    
    system(paste("/path_to_dir_/scwrl4/Scwrl4 -i ", input_file, " -o ",
                 paste(paste(current_dir,"../Files/structures/",sep="/"),"/struct_step_",actual_step,".pdb",sep=""),
                 " -s ",paste(current_dir,"../Files/sequences/",sep = "/"),mutated_seq," -h",sep=""))
    
    setwd(current_dir)
    
    ##########################################################################################################
    ##########################################################################################################
    
    ####################################
    ## New Descriptor Calculation ######
    ####################################
    ####################################
    
    system("rm ../Files/Surf_tmp/*")
    system("rm ../Files/Surf_BS_2/*")
    
    # new surface calculator#
    
    data_path_mut <- "../Files/structures/"
    data_name_mut <- paste("struct_step_",actual_step,".pdb",sep="")
    output_path <- "../Files/Surf_tmp/"
    output_name_mut <- substring(data_name_mut,1,nchar(data_name_mut)-4)
    
    get.surf(data_path_mut,data_name_mut,1.4,output_path, output_name_mut )
    surf_mut <- read.csv(paste(output_path,output_name_mut,".csv",sep=""),stringsAsFactors = F)
    
    # New BS identification #
    
    res_id <- substring(surf_mut[,1],1,3)
    res_no <- substring(surf_mut[,1],5,nchar(surf_mut[,1])-1)
    res_chain <- substring(surf_mut[,1],nchar(surf_mut[,1]),nchar(surf_mut[,1]))
    
    surf_mut_BS <- surf_mut[paste(res_no,res_chain) %in%  paste( Interface_cycle$resno, Interface_cycle$chain),]
    write.csv(surf_mut_BS,paste("../Files/Surf_BS_2/",paste(substring(data_name_mut,1,nchar(data_name_mut)-4),"_bs.csv",sep=""),sep=""),quote = F,row.names = F)
    
    # Zernike Calculations #
    print("++++ Zernike calculations ++++")
    system("ls ../Files/Surf_BS_2/ > auxiliary_files/names_2.txt")
    
    system("python3 Z2D_verso_-1.py auxiliary_files/names_2.txt")
    
    name_file_zernike_mut <- list.files("../Files/ZernDescr/", pattern = paste(substring(data_name_mut,1,nchar(data_name_mut)-4),"_",sep=""))
    
    
    ZD_fixed <- read.csv(name_Zern_fixed, header=F)
    ZD_mut <- read.csv(paste("../Files/ZernDescr/", name_file_zernike_mut,sep=""), header=F)
    
    Zern_dist_mut <- as.numeric(dist(rbind(ZD_fixed[,1],ZD_mut[,1])))
    
    #### Electrostatics Calculation #############################################################
    ##############################################################################################
    print('++++ Electrostatics calculation ++++')
    system(paste("pdb2pqr", paste("../Files/structures/", data_name_mut, sep=""),
                 paste("../Files/PQR/",substring(data_name_mut,1,nchar(data_name_mut)-3),"pqr",sep=""),
                 "--ff=CHARMM","--keep-chain", sep = " "))
    system("rm ../Files/PQR/*.log")
    
    name_PQR_mut <- list.files("../Files/PQR/", pattern = paste(substring(data_name_mut,1,nchar(data_name_mut)-4),".pqr",sep=""))
    
    PQR_fixed <- read.pdb(name_PQR_fixed)
    PQR_1 <- PQR_fixed$atom
    PQR_mut <- read.pdb(paste("../Files/PQR/",name_PQR_mut,sep=""))
    PQR_2 <- PQR_mut$atom
    
    ########################################
    # Coarse-Grained construction
    
    residues_PQR_1<- unique(paste(PQR_1$chain,PQR_1$resno,PQR_1$resid)) 
    
    PQR_1_CG <- c()
    for(k in 1:length(residues_PQR_1)){
      aus <- strsplit(residues_PQR_1[k]," ")[[1]]
      single_res <- PQR_1[PQR_1$chain %in% aus[1] & PQR_1$resno %in% aus[2] & PQR_1$resid %in% aus[3],]
      
      single_res_backbone <- single_res[single_res$elety %in% c("N","CA","C","O"),]
      row_bb <- single_res_backbone[2,]
      row_bb[,9:11] <- apply(single_res_backbone[,9:11],2,mean)
      row_bb[,12] <- sum(single_res_backbone[,12])
      
      single_res_sidechain <- single_res[! single_res$elety %in% c("N","CA","C","O"),]
      row_sc <- single_res_sidechain[1,]
      row_sc$elety <- "B"
      row_sc[,9:11] <-apply(single_res_sidechain[substring(single_res_sidechain$elety,1,1) != "H",9:11], 2,mean) 
      if(aus[3] =="GLY"){row_sc[,9:11] <- apply(single_res_sidechain[,9:11], 2,mean)}
      row_sc[,12] <- sum(single_res_sidechain[,12])
      
      
      PQR_1_CG <- rbind(PQR_1_CG,row_bb,row_sc)
      
    }
    
    residues_PQR_2<- unique(paste(PQR_2$chain,PQR_2$resno,PQR_2$resid)) 
    
    PQR_2_CG <- c()
    for(k in 1:length(residues_PQR_2)){
      aus <- strsplit(residues_PQR_2[k]," ")[[1]]
      single_res <- PQR_2[PQR_2$chain %in% aus[1] & PQR_2$resno %in% aus[2] & PQR_2$resid %in% aus[3],]
      
      single_res_backbone <- single_res[single_res$elety %in% c("N","CA","C","O"),]
      row_bb <- single_res_backbone[2,]
      row_bb[,9:11] <- apply(single_res_backbone[,9:11],2,mean)
      row_bb[,12] <- sum(single_res_backbone[,12])
      
      single_res_sidechain <- single_res[! single_res$elety %in% c("N","CA","C","O"),]
      row_sc <- single_res_sidechain[1,]
      row_sc$elety <- "B"
      row_sc[,9:11] <-apply(single_res_sidechain[substring(single_res_sidechain$elety,1,1) != "H",9:11], 2,mean) 
      if(aus[3] =="GLY"){row_sc[,9:11] <- apply(single_res_sidechain[,9:11], 2,mean)}
      row_sc[,12] <- sum(single_res_sidechain[,12])
      
      
      PQR_2_CG <- rbind(PQR_2_CG,row_bb,row_sc)
      
    }
    
    # we defined electrostatics interface using 12Å threshold
    
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
    
    En_INT_mut <- sum(En_tot[,3])
    
    
    #####################################################################
    ####### hydrophobicity calculation ##################################
    #####################################################################
    print("++++ Hydropathy calculations ++++")
    surf_fixed_BS <- read.csv(name_surf_BS_fixed)
    
    
    
    inv_res <- table(surf_fixed_BS$Res)
    inv_res <- inv_res[inv_res != 0]
    names(inv_res) <- substring(names(inv_res),1,3)
    inv_res_hyd <- hyd_scale[names(inv_res)]
    hydro_patch_fixed <- sum(inv_res_hyd * inv_res)
    hydro_patch_fixed <- round(hydro_patch_fixed / nrow(surf_fixed_BS),3)
    
    inv_res <- table(surf_mut_BS$Res)
    inv_res <- inv_res[inv_res != 0]
    names(inv_res) <- substring(names(inv_res),1,3)
    inv_res_hyd <- hyd_scale[names(inv_res)]
    hydro_patch_mut <- sum(inv_res_hyd * inv_res)
    hydro_patch_mut <- round(hydro_patch_mut / nrow(surf_mut_BS),3)
    
    
    hydro_similarity <- abs(hydro_patch_fixed - hydro_patch_mut)
    
    
    ###########################################################
    ###########################################################
    ### mutation evaluation ###################################
    ###########################################################
    
    
    new_descriptord <- c(round(Zern_dist_mut,3),round(En_INT_mut,3), round(hydro_similarity,3))
    
    Zern_bal <- Zern_dist_mut - dist_zern_cycle
    EN_bal <- En_INT_mut - EN_tot_cycle
    hydro_bal <- hydro_similarity - dist_hydro_cycle
    
    seq_new <- strsplit(seq_valutanda,"")
    seq_new <- toupper(aaa(toupper(seq_new[[1]])))
    names(seq_new) <- names(seq_cycle)
    
    DELTA_MUT <- 0
    if( (length(which(! seq_new == seq_ini)) - length(which(! seq_cycle == seq_ini))) == 1){
      DELTA_MUT <- (length(which(! seq_new == seq_ini)))
    }
    DELTA_MUT_tot <- length(which(! seq_new == seq_ini))
    
    ### Cost function definition
    Delta_Ene_tot <- A*Zern_bal + B*EN_bal + C*hydro_bal + D*DELTA_MUT^(2)
    
    if(Delta_Ene_tot <= 0){
      
      seq_cycle <- seq_new
      dist_zern_cycle <- Zern_dist_mut
      EN_tot_cycle <- En_INT_mut 
      dist_hydro_cycle <- hydro_similarity
      Interface_cycle[pos_res_mut,1] <- res_to_ins
      acceptance <- 1
      
    }
    if(Delta_Ene_tot > 0){
      probability <- exp(- Beta_range[Beta]*Delta_Ene_tot)
      probCasual <- runif(1,0,1)
      
      if(probCasual < probability){
        seq_cycle <- seq_new
        dist_zern_cycle <- Zern_dist_mut
        EN_tot_cycle <- En_INT_mut 
        dist_hydro_cycle <- hydro_similarity
        Interface_cycle[pos_res_mut,1] <- res_to_ins
        acceptance <- 1
      }
      if(probCasual > probability){
        seq_cycle <- seq_cycle
        dist_zern_cycle <- dist_zern_cycle
        EN_tot_cycle <- EN_tot_cycle
        dist_hydro_cycle <- dist_hydro_cycle
        acceptance <- 0
      }
      
    }
    
    row_to_save <- c(actual_step, Beta_range[Beta], as.character(old_res), num_res, new_res,
                      round(as.numeric(dist_zern_cycle),3), round(as.numeric(Zern_dist_mut),3),
                      round(as.numeric(EN_tot_cycle),3), round(as.numeric(En_INT_mut),3),
                      round(as.numeric(dist_hydro_cycle),3), round(as.numeric(hydro_similarity),3),
                      DELTA_MUT_tot,acceptance)
    
    df_tot <- rbind(df_tot, row_to_save)
    colnames(df_tot) <- c("NoStep", "Beta", "Substituted Residue", "ResNo", "Inserted Residue",
                             "Selected Zernike", "Proposed Zernike", "Selected Elec", "Proposed Elec", 
                             "Selected Hydro", "Proposed Hydro", "Mutations number", "Acceptance")
    
    write.csv(df_tot,"../Files/MonteCarlo_record.csv", quote = F, row.names = F)
  }
}
