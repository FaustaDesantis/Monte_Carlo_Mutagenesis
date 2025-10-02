#library(tidyverse)

get_mol_weig <- function(x){
  
  nomi_atomi <- c("H","C", "O", "N", "I", "S", "P", "F", "B", "E")
  peso_atomi <- c(1.00,12.01, 15.99, 14.00, 126.90, 32.06, 30.97, 18.99, 10.81, 12.01)
  
  peso_mol <- 0
  for (j in 1:nrow(x)) {
    if(substring(x$elety[j], 1, 1) %in% nomi_atomi){
      peso_mol <- peso_mol + peso_atomi[nomi_atomi %in% substring(x$elety[j], 1, 1)]
    }
    if(!substring(x$elety[j], 1, 1) %in% nomi_atomi){
      peso_mol <- peso_mol + 12.01
    }
    
  }
  
  return(peso_mol)
}


my.get.seq <- function(data_pdb){
  CA <- data_pdb[data_pdb$elety == "CA",]
  seq <- CA$resid
  num_seq <- CA$resno
  names(seq) <- num_seq
  
  return(seq)
}


pdb_format<-function(x, name) {
  x <- x[,1:13]
  matri<-x[,9:11]
  
  coord<-as.vector(apply(matri, 1, function(x) list=c(x[1], x[2], x[3])))
  
  if(length(x[is.na(x$o),]$o) > 0){
    x[is.na(x$o),]$o <- 1
  }
  if(length(x[is.na(x$b),]$b) > 0){
    x[is.na(x$b),]$b <- 1
  }
  write.pdb(xyz = coord, type=as.character(x$type), resno=as.character(x$resno), resid=as.character(x$resid), eleno=as.character(x$eleno), elety=as.character(x$elety), chain=as.character(x$chain), insert=as.character(x$insert),o = as.character(x$o), b = as.character(x$b), file=name)
  
}

clear_dms_file <- function(surf) {
  
  surf_final <- surf %>% 
    mutate(B = gsub("[[:punct:]]", "", as.character(B))) %>%    # keep only alphanumeric characters
    filter(str_sub(G, 1, 1) == "S") %>%                         # keep only surface points
    unite(Res, c(1,2)) %>%                                      # join res name and number
    select(c(1, 3, 4, 5, 8, 9, 10)) %>%                         # select useful columns
    set_names(c("Res","x","y","z","Nx","Ny","Nz"))              # rename columns
  
  return(surf_final)
}

clear_dms_file <- function(surf) {
  
  surf <- as.data.frame(surf)
  
  # keep only alphanumeric characters
  surf$B <- gsub("[[:punct:]]", "", as.character(surf$B))
  
  # keep only surface points
  surf <- surf[substring(surf$G, 1, 1) == "S", ]
  
  # select and rename columns
  surf_final <- cbind(paste(surf[, 1], surf[, 2], sep="_"), surf[, 4:6], surf[, 9:11])
  colnames(surf_final) <- c("Res","x","y","z","Nx","Ny","Nz")
  
  return(surf_final)
}
get.surf <- function(data_path, data_name, probe_radius, out_path, out_name){
  
  system(paste("dms",paste(data_path,data_name,sep=""), "-w", probe_radius ,"-n -a -o", paste(out_path,substring(data_name,1,nchar(data_name)-3), "dms", sep=""), sep = " "))
  
  
  surf_name <- list.files(out_path, pattern = paste(substring(data_name,1,nchar(data_name)-3),"dms",sep=""))
  
  verifica <- read.csv(paste(out_path,surf_name,sep=""), sep="", col.names = LETTERS[1:11], header = F, colClasses = "character")
  
  verifica <- as.matrix(verifica)
  verifica[is.na(verifica[,3]),3] <- "X"
  verifica <- as.data.frame(verifica)
  
  surf <- verifica
  
  da_verificare <- which(nchar(as.character(verifica$C))> 3)
  if(length(da_verificare)>0 & length(da_verificare) < nrow(verifica)){
    da_tenere <- verifica[-da_verificare,]
    da_cambiare <- verifica[da_verificare,]
    colIV <- substring(da_cambiare[,3],nchar(as.character(da_cambiare[,3]))-7,nchar(as.character(da_cambiare[,3])))
    colIII <- substring(da_cambiare[,3],1,nchar(as.character(da_cambiare[,3]))-8)
    da_inserire <- cbind(da_cambiare[,1:2], colIII,colIV,da_cambiare[,4:10])
    colnames(da_inserire) <- LETTERS[1:11]
    surf <- rbind(da_tenere,da_inserire)
  }
  if(length(da_verificare)>0 & length(da_verificare) == nrow(verifica)){
    da_cambiare <- verifica[da_verificare,]
    colIV <- substring(da_cambiare[,3],nchar(as.character(da_cambiare[,3]))-7,nchar(as.character(da_cambiare[,3])))
    colIII <- substring(da_cambiare[,3],1,nchar(as.character(da_cambiare[,3]))-8)
    da_inserire <- cbind(da_cambiare[,1:2], colIII,colIV,da_cambiare[,4:10])
    colnames(da_inserire) <- LETTERS[1:11]
    surf <- rbind(da_inserire)
  }
  
  surf_final <- clear_dms_file(surf)
  
  system(paste("rm ",out_path,surf_name,sep=""))
  write.csv(surf_final, paste(out_path,out_name, ".csv",sep=""), row.names = F, quote=F)
  
}
