############################################################################################
# align sepcies 
## functions 
align_specy_table <- function(path, replace_pattern2='', ref_name, prefix='')
{
  read_as_list <- function(path, seq='', prefix=''){
    ipak <- function(pkg){
      new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
      if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
      sapply(pkg, require, character.only = TRUE)
    }
    packages <- c('dplyr', 'purrr', 'fs')
    ipak(packages)
    mini_pre.bed <- list()
    bed_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
    bed_path <- as_fs_path(bed_path)
    mini_pre.bed <- bed_path %>% map(.f = function(path){read.delim(path, header = FALSE, stringsAsFactors = F, sep = seq)})
    names(mini_pre.bed) <- bed_path
    pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
    for(i in 1:2)
    {names(mini_pre.bed) <- sub(names(mini_pre.bed), pattern = pattern[i], replacement = '')}
    return(mini_pre.bed)
    sapply(packages, require, character.only = TRUE)
    lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
  }
  A <- read_as_list(path = path, prefix = prefix)
  class <- as.data.frame(matrix(nrow = length(A), ncol = 3))
  for(i in seq_along(A)) {
    B <- as.data.frame(A[[i]])
    class[i,'libs'] <- names(A)[i]
    unmap <- length(B[,1][B[,1] %in% '*'])
    class[i,1] <- unmap
    virus <- length(B[,1][B[,1] %in% ref_name])
    class[i,2] <- virus
    class[i,3] <- nrow(B)-unmap-virus
  }
  colnames(class)[1:3] <- c('Unmap', 'Virus', 'Host')
  return(class)
}



class <- align_specy_table(path = 'Analysis/MHV/bed/host/genome/'
                           , ref_name = 'NC_001846.1')

write.csv(x = class, file = '~/Analysis/MHV/result/align_species.csv')



############################################################################################
#viral reads fragmens distribution
## functions
read_as_list <- function(path, seq='', prefix='')
{
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('dplyr', 'purrr', 'fs')
  ipak(packages)
  mini_pre.bed <- list()
  bed_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
  bed_path <- as_fs_path(bed_path)
  mini_pre.bed <- bed_path %>% map(.f = function(path){read.delim(path, header = FALSE, stringsAsFactors = F, sep = seq)})
  names(mini_pre.bed) <- bed_path
  pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
  for(i in 1:2)
  {names(mini_pre.bed) <- sub(names(mini_pre.bed), pattern = pattern[i], replacement = '')}
  return(mini_pre.bed)
  sapply(packages, require, character.only = TRUE)
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}

fragment_table <- function(list, max_frag=5, re_organ=T, proportion=T) {
  frag <- as.data.frame(matrix('', nrow = length(bed),ncol = max_frag+2))
  colnames(frag) <- c(seq(1,max_frag), paste('>=', (max_frag+1), sep = ''), 'libs')
  index <- c(seq(1,max_frag), paste('>=', (max_frag+1), sep = ''))
  for(k in seq_along(list)) {
    BED <- list[[k]]
    BED[,2] <- BED[,2]+1
    BED <- BED[order(BED$V4), ]
    BED$inedx <- seq(1, nrow(BED))
    name <- as.data.frame(table(BED$V4))
    A <- as.data.frame(table(name$Freq), stringsAsFactors = F)
    A$Var1 <- as.numeric(A$Var1)
    row <- nrow(A[A$Var1 <= max_frag,])
    A[(row+1), 'Freq'] <- colSums(as.matrix(A[(row+1):nrow(A), 'Freq']))
    A[(row+1),"Var1"] <- paste('>=', (max_frag+1))
    A <- A[-c((row+2):nrow(A)),]
    for(j in seq(1,max_frag+1)){
      value <- A$Freq[A$Var1 %in% index[j]]
      frag[k,j] <- ifelse(identical(value, numeric(0)), yes = 0, no = value)
      frag[k,'libs'] <- names(list)[k]
    }
    if(proportion==T) {
      frag[k,1:(max_frag+1)] <- round(as.numeric(frag[k,1:(max_frag+1)])/colSums(as.matrix(as.numeric(frag[k,1:(max_frag+1)]))), digits = 2)
    }
    
  }
  rownames(frag) <- frag[,'libs']
  if(re_organ==T) {
    A <- frag %>% 
      dplyr::mutate(ID = 1:n()) %>%  
      tidyr::gather(key, value, c(seq(1,(max_frag+1)))) 
    A$key <- factor(A$key, levels =  c(seq(1,max_frag), paste('>=', (max_frag+1), sep = '')))
    A$value <- as.numeric(A$value)
    return(A)
  }
  else{
    return(frag)
  }
}

## ploting

bed <- read_as_list(path='~/Analysis/MHV/bed/cell_liver/v2/') 

A <- fragment_table(list = bed, max_frag = 40 , proportion = F, re_organ = F)
write.csv(A, file = '~/Analysis/MHV/result/MHV_viral_fragment.csv')




############################################################################################
# sgRNA classification table
## functions
read_as_list <- function(path, seq='', prefix='')
{
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('dplyr', 'purrr', 'fs')
  ipak(packages)
  mini_pre.bed <- list()
  bed_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
  bed_path <- as_fs_path(bed_path)
  mini_pre.bed <- bed_path %>% map(.f = function(path){read.delim(path, header = FALSE, stringsAsFactors = F, sep = seq)})
  names(mini_pre.bed) <- bed_path
  pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
  for(i in 1:2)
  {names(mini_pre.bed) <- sub(names(mini_pre.bed), pattern = pattern[i], replacement = '')}
  return(mini_pre.bed)
  sapply(packages, require, character.only = TRUE)
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}


# function 5: format organization depending on split frag
specific_frag_format_organ <- function(list, frag)
{
  specific_fragment <- list()
  for(j in seq_along(list))
  {
    BED <- list[[j]]
    BED[,2] <- BED[,2]+1
    name <- as.data.frame(table(BED[,4]), stringsAsFactors =F)
    split <- as.data.frame(table(name$Freq), stringsAsFactors = F)
    name <- name$Var1[name$Freq %in% frag]
    if(length(name) >0 )
    {
      BED <- BED[BED[,4] %in% name,]
      BED <- BED[order(BED[,4], decreasing = F),]
      fragment <- list()
      if(frag==1) {
        B <- BED
      }
      else if(frag!=1) {
        for(i in 1:frag) {
          if(i ==1)
          {A <- BED[seq(i,nrow(BED), frag),][,1:3]}
          if(i ==frag)
          {A <- BED[seq(i,nrow(BED), frag),][,2:6]}
          if(i>1 & i<frag)
          {A <- BED[seq(i,nrow(BED), frag),][,2:3]}
          fragment[[i]] <- A
        }
        B <- as.data.frame(matrix(data = 0, nrow = nrow(fragment[[1]]), ncol = 0))
        for(i in rev(seq_along(fragment))) {
          A <- fragment[[i]]
          B <- cbind(A,B)
        }
      }
      
      colnames(B)[seq(2,(ncol(B)-3), 2)] <- paste('start', seq(1,frag,1), sep = '')
      colnames(B)[seq(3,(ncol(B)-2), 2)] <- paste('end', seq(1,frag,1), sep = '')
      specific_fragment[[j]] <- B
      names(specific_fragment)[j] <- names(list)[j]
    }
    else(next)
  }
  return(specific_fragment)
}



# Function3: identify the sgRNA (for files from Function5)
## gff_file have to be the required format
## sgmRNA_BA <- sgmRNA_boundary_allow
specific_reads_sgRNA_identify <- function(list, gap=50, gff_path, sgmRNA_BA=50
                                          , leader, UTR5, UTR3, frag, out_prefix=''
                                          , whole_genome_tolerance=20, length
                                          , isrecom=F, write=F, write_path=getwd())
{ 
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('plyr', 'dplyr')
  ipak(packages)
  
  gff <- read.csv(file = gff_path, header = F)
  
  if(frag==2)
  {
    output_list <- list()
    for(i in seq_along(list)) {
      A <- list[[i]]
      A$gap <- A[,'start2']-A[,'end1']
      B <- A
      B$small_gap <- (B[,'gap'] < gap)
      small_gap <- B[B$small_gap ==T,-c(ncol(B))]
      if(nrow(small_gap)>0) {
        small_gap$class <- c('small_gap')
        small_gap$sub_class <- c('')
      }
      if(nrow(small_gap)==0) { 
        small_gap <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      
      B <- B[B$small_gap ==F,-c(ncol(B))]
      B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
      UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
      if(nrow(UTR3_F)>0) {
        C <- UTR3_F
        C$UTR5 <- C[,'start1'] <= UTR5
        D <- C[C$UTR5 ==T,-c(ncol(C))]
        if(nrow(D)>0) {
          D$class <- c('Truncated Di')
          D$sub_class <- c('delete 3\'UTR Di')
        }
        if(nrow(D)==0) {
          D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
        E <- C[C$UTR5 ==F,-c(ncol(C))]
        if(nrow(E)>0)  {
          E$class <- c('Truncated Di')
          E$sub_class <- c('delete 3\' 5\'UTR Di')
        }
        if(nrow(E)==0) {
          E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
        delete_UTR_Di <- rbind(D,E)
      }
      
      B <- B[B$UTR3 ==T,-c(ncol(B))]
      if(nrow(B)>0) {
        C <- B
        C$UTR5 <- C[,'start1'] < UTR5
        delete_5UTR_Di <- C[C$UTR5 ==F,-c(ncol(C))]
        if(nrow(delete_5UTR_Di)>0) {
          delete_5UTR_Di$class <- c('Truncated Di')
          delete_5UTR_Di$sub_class <- c('delete 5\'UTR Di')
        }
        if(nrow(delete_5UTR_Di)==0) {
          delete_5UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
      }
      
      C <- C[C$UTR5 ==T, -c(ncol(C))]
      C$leader <- C$end1 <= leader
      sgmRNA <- C[C$leader == T, -c(ncol(C))]
      if(nrow(sgmRNA)>0) {
        sgmRNA$class <- c('sgmRNA')
        sgmRNA$sub_class <- c('')
      }
      if(nrow(sgmRNA)==0) {sgmRNA <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      
      DiRNA <- C[C$leader ==F,-c(ncol(C))]
      if(nrow(DiRNA)>0) {
        DiRNA$class <- c('DiRNA')
        DiRNA$sub_class <- c('')
      }
      if(nrow(DiRNA)==0) {
        DiRNA <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      class <- rbind(small_gap, delete_UTR_Di, delete_5UTR_Di, DiRNA, sgmRNA)
      
      B <- class
      rownames(B) <- seq(1,nrow(B))
      if(nrow(B[B$class %in% 'sgmRNA', ])>0) {
        A <- B[!(B$class %in% 'sgmRNA'),]
        sgmRNA <- list()
        D <- B[B$class %in% c('sgmRNA'), ]
        if(nrow(D)>0) {
          for(j in seq_len(nrow(gff))) {
            D$sub_class_bulin <- (gff[j,4]-sgmRNA_BA)<D[,paste('start', 2, sep = '')] & D[,paste('start', 2, sep = '')]<(gff[j,4]+sgmRNA_BA)
            E <- D[D$sub_class_bulin==T,]
            E <- E[,-c(ncol(E))]
            if(nrow(E)>0) {
              E$sub_class <- gff[j,9]
            }
            sgmRNA[[j]] <- E
            D <- D[D$sub_class_bulin ==F,]
          }
          D <- D[D$sub_class_bulin ==F,]
          D <- D[,-c(ncol(D))]
          if(nrow(D)>0) {
            D$sub_class <- c('uncanonical sgRNAs')
          }
          sgmRNA[[(length(sgmRNA)+1)]] <- D
          D <- do.call(rbind, sgmRNA)
          B <- rbind(A,D)
          if(isTRUE(isrecom)) {
            B <- B[order(B$freq, decreasing = T),]
          }
        }
      }
      output_list[[i]] <- B
      names(output_list)[[i]] <- names(list)[[i]]
    }
  }
  
  if(frag>2) {
    output_list <- list()
    for(i in seq_along(list)) {
      frag_other <- list[[i]]
      if(nrow(frag_other)>0) {
        B <- frag_other
        dfgap <- as.data.frame(matrix(data = '', nrow = nrow(B), ncol = 1))
        for(j in 2:frag) {
          dfgap[,(j-1)] <- frag_other[,paste('end', j, sep = '')]-frag_other[,paste('start', (j-1), sep = '')]
        }
        dfgap$min <- apply(dfgap, MARGIN = 1,FUN = min)
        B$gap <- dfgap$min
        B$small_gap <- B[,'gap'] <gap
        small_gap <- B[B$small_gap ==T,-c(ncol(B))]
        if(nrow(small_gap)>0) {
          small_gap$class <- c('small_gap')
          small_gap$sub_class <- c('')
        }
        if(nrow(small_gap)==0) {
          small_gap <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        
        B <- B[B$small_gap ==F,-c(ncol(B))]
        B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
        UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
        if(nrow(UTR3_F)>0) {
          C <- UTR3_F
          C$UTR5 <- C[,'start1'] <= UTR5
          D <- C[C$UTR5 ==T,-c(ncol(C))]
          if(nrow(D)>0) {
            D$class <- c('Truncated Di')
            D$sub_class <- c('delete 3\'UTR Di')
          }
          if(nrow(D)==0) {
            D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          E <- C[C$UTR5 ==F,-c(ncol(C))]
          if(nrow(E)>0) {
            E$class <- c('Truncated Di')
            E$sub_class <- c('delete 3\' 5\'UTR Di')
          }
          if(nrow(E)==0) {
            E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          delete_UTR_Di <- rbind(D,E)
        }
        if(nrow(UTR3_F)==0) {
          delete_UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        
        UTR3_T <- B[B$UTR3 ==T,-c(ncol(B))]
        if(nrow(UTR3_T)>0) {
          C <- UTR3_T
          C$UTR5 <- C[,'start1'] > UTR5
          D <- C[C$UTR5 ==T,-c(ncol(C))]
          if(nrow(D)>0) {
            D$class <- c('Truncated Di')
            D$sub_class <- c('delete 5\'UTR Di')
          }
          if(nrow(D)==0) { 
            D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          E <- C[C$UTR5 ==F,-c(ncol(C))]
          if(nrow(E)>0) {
            E$class <- c('DiRNA')
            E$sub_class <- c('')
          }
          if(nrow(E)==0) {
            E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          delete_UTR5_Di <- rbind(D,E)
        }
        if(nrow(UTR3_T)==0) {
          delete_UTR5_Di <-  as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        class <- rbind(delete_UTR5_Di, delete_UTR_Di, small_gap)
        if(isTRUE(isrecom)) {class <- class[order(class$freq, decreasing = T),]
        }
      }
      output_list[[i]] <- class
      names(output_list)[[i]] <- names(list)[[i]]
    }
  }
  
  if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False')) {
    return(output_list)
  }
  if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true')) {
    if(isFALSE(grepl(write_path, pattern = '&/'))){
      write_path <- paste(write_path,'/', sep = '')
    }
    lapply(seq_along(output_list)
           , function(i) 
             write.csv(output_list[[i]]
                       , file = paste0(write_path
                                       ,names(output_list[i]),'_', frag, 'fragment_', out_prefix, ".csv"),
                       row.names = FALSE, quote = F))
  }
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}



# function6 : construct recombination point sort and new seq name
specific_frag_recom_table <- function(list, frag, out_prefix='', write=F, write_path=getwd())
{
  recom <- list()
  for(i in seq_along(list)) {
    B <- list[[i]]
    B[B$sub_class %in% '','sub_class'] <- '.'
    colNames <-  colnames(B)[c(2:(frag*2+1),(frag*2+6), (frag*2+7))] # could be any number of column names here
    B$recom <- do.call(paste, c(B[colNames], sep="-"))
    recom_table <- as.data.frame(table(B$recom), stringsAsFactors = F)
    recom_table <- recom_table[order(recom_table$Freq, decreasing = T),]
    freq <- recom_table$Freq
    recom_table <- as.data.frame(do.call(rbind, strsplit(recom_table$Var1, split = '-')))
    recom_table$Freq <- freq
    recom[[i]] <- recom_table
  }
  names(recom) <- names(list)
  
  if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False')) {
    return(recom)
  }
  if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true')) {
    if(isFALSE(grepl(write_path, pattern = '&/'))){
      write_path <- paste(write_path,'/', sep = '')
    }
    lapply(seq_along(recom)
           , function(i) 
             write.csv(recom[[i]]
                       , file = paste0(write_path
                                       ,names(recom[i]),'_frag', frag, '_recombination',out_prefix, ".csv"),
                       row.names = FALSE, quote = F))
  }
}

# function : one fragment frequency table
one_frag_recom_table <- function(list, frag, out_prefix='', write=F, write_path=getwd())
{
  recom <- list()
  for(i in seq_along(list)) {
    B <- list[[i]]
    B[B$sub_class %in% '','sub_class'] <- '.'
    colNames <-  colnames(B)[c(2,3,8,9)] # could be any number of column names here
    B$recom <- do.call(paste, c(B[colNames], sep="-"))
    recom_table <- as.data.frame(table(B$recom), stringsAsFactors = F)
    recom_table <- recom_table[order(recom_table$Freq, decreasing = T),]
    freq <- recom_table$Freq
    recom_table <- as.data.frame(do.call(rbind, strsplit(recom_table$Var1, split = '-')))
    recom_table$Freq <- freq
    recom[[i]] <- recom_table
  }
  names(recom) <- names(list)
  
  if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False')) {
    return(recom)
  }
  if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true')) {
    if(isFALSE(grepl(write_path, pattern = '&/'))){
      write_path <- paste(write_path,'/', sep = '')
    }
    lapply(seq_along(recom)
           , function(i) 
             write.csv(recom[[i]]
                       , file = paste0(write_path
                                       ,names(recom[i]),'_frag', frag, '_recombination',out_prefix, ".csv"),
                       row.names = FALSE, quote = F))
  }
}


## frag2-frag5
library(dplyr)
bed <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/')
for(j in seq(2,5,1)) {
  bed%>%
    specific_frag_format_organ(frag = j) %>% 
    specific_reads_sgRNA_identify(gap = 50, length = 31335
                                  , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                  , leader = 72, UTR5 = 210, UTR3 = 31035, frag = j) %>%
    specific_frag_recom_table(frag = j, write = T, out_prefix = 'MHV'
                              , write_path = '~/Analysis/MHV/result/recombination/')
}


## frag1
lib_recom_matrix <- function(lib_list, frag, libs) {
  recom_table_list <- list()
  recom_list <- list()
  inter_count <- ''
  x=0
  for(i in seq_along(lib_list)) {
    lib1 <- lib_list[[i]]
    lib1 <- lib1[!(lib1$class %in% 'small_gap'), ]
    lib1[rownames(lib1[lib1$sub_class %in% '',]),'sub_class'] <- ' '
    if(frag !=1) {
      colNames <-  colnames(lib1)[c(2:(frag*2+1),(frag*2+6), (frag*2+7))]
    }
    if(frag ==1) {
      colNames <-  colnames(lib1)[c(2,3,8,9)]
    }
    lib1$recom <-  do.call(paste, c(lib1[colNames], sep="-"))
    lib1 <- as.data.frame(table(lib1$recom), stringsAsFactors = F)
    recom_list[[i]] <- lib1$Var1
    recom_table_list[[i]] <- lib1
  }
  names(recom_list) <- names(lib_list)
  names(recom_table_list) <- names(lib_list)
  
  all_recom <- Reduce(union, recom_list)
  lib_matrix <- as.data.frame(matrix(data = 0, nrow = length(all_recom), ncol = libs))
  lib_matrix$index <- seq(1, nrow(lib_matrix))
  rownames(lib_matrix) <- all_recom; colnames(lib_matrix) <- names(lib_list)
  
  for(i in seq_along(recom_table_list)) {
    B <- recom_table_list[[i]]
    B <- B[order(B$Var1, decreasing = T), ]
    index <- lib_matrix[rownames(lib_matrix) %in% B$Var1, (ncol(lib_matrix))]
    C <- lib_matrix[lib_matrix[,ncol(lib_matrix)] %in% index,]
    C <- C[order(rownames(C), decreasing = T), ]
    lib_matrix <- lib_matrix[!(lib_matrix[,ncol(lib_matrix)] %in% index), ]
    print(table(rownames(C)== B$Var1))
    C[,i] <- B$Freq
    lib_matrix <- rbind(lib_matrix,C)
  }
  lib_matrix <- lib_matrix[, -c(ncol(lib_matrix))]
}


lin_frag1_identification <- function(list, recom_matrix, gff_path, libs
                                     , UTR5, UTR3, length, sgmRNA_BA=50
                                     , whole_genome_tolerance=20) {
  unc <- list()
  for(j in seq(1,libs,2)) {
    B <- recom_matrix[,j:(j+1)]
    C <- as.data.frame(ifelse(test = B>=2, yes = 1, no = 0))
    C$con <- rowSums(as.matrix(C[,1:2]))
    C <- C[C$con > 0, ]
    C <- C[grep(rownames(C), pattern = 'uncanonical'), ]
    C <- C[C$con >1, ]
    C <- as.data.frame(do.call(rbind, strsplit(rownames(C), split = '-')))
    unc[[j]] <- rep(as.numeric(unique(C$V3)))
    unc[[j+1]] <- rep(as.numeric(unique(C$V3)))
  }
  
  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('plyr', 'dplyr')
  ipak(packages)
  
  gff <- read.csv(file = gff_path, header = F)
  
  output_list <- list()
  for(i in seq_along(list))  {
    f1 <- list[[i]]
    f1$length <- f1[,3]-f1[,2]
    f1$complete <- f1$length > (length-whole_genome_tolerance)
    complete <- f1[f1$complete %in% TRUE,-c(ncol(f1))]
    if(nrow(complete)>0){
      complete$class <- c('complete genome')
    }
    if(nrow(complete)==0){
      complete <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- f1[f1$complete %in% FALSE,-c(ncol(f1))]
    f1$utr3 <- f1$end> UTR3
    delete_utr3 <- f1[f1$utr3 %in% FALSE,-c(ncol(f1))]
    if(nrow(delete_utr3)>0){
      delete_utr3$utr5 <- delete_utr3$start < UTR5
      delete_utr3_utr5 <- delete_utr3[delete_utr3$utr5 %in% FALSE, -c(ncol(delete_utr3))]
      if(nrow(delete_utr3_utr5)>0){
        delete_utr3_utr5$class <- c('delete 3\'UTR and 5\' UTR genome')
      }
      if(nrow(delete_utr3_utr5)==0)
      {delete_utr3_utr5 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
      
      delete_utr3 <- delete_utr3[delete_utr3$utr5 %in% TRUE, -c(ncol(delete_utr3))]
      if(nrow(delete_utr3)>0){
        delete_utr3$class <- c('delete 3\'UTR genome')
      }
      if(nrow(delete_utr3)==0)
      {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    }
    if(nrow(delete_utr3)==0)
    {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- f1[f1$utr3 %in% TRUE,-c(ncol(f1))]
    f1$utr5 <- f1$start < UTR5
    delete_utr5 <- f1[f1$utr5 %in% FALSE,-c(ncol(f1))]
    if(nrow(delete_utr5)>0)
    {
      delete_utr5$class <- c('delete 5\' terminal')
    }
    if(nrow(delete_utr5)==0)
    {delete_utr5 <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    near_complete <- f1[f1$utr5 %in% T,-c(ncol(f1))]
    if(nrow(near_complete)>0)
    {
      near_complete$class <- c('near complete genome')
    }
    if(nrow(near_complete)==0)
    {near_complete <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- rbind(delete_utr3, delete_utr3_utr5, delete_utr5, near_complete)
    output_list[[i]] <- f1
    names(output_list)[[i]] <- names(list)[[i]]
  }
  
  for(k in 1:length(output_list)) {
    B <- output_list[[k]]
    rownames(B) <- seq(1,nrow(B))
    if(nrow(B[B$class %in% 'delete 5\' terminal', ])>0) {
      A <- B[!(B$class %in% 'delete 5\' terminal'),]
      sgmRNA <- list()
      A$sub_class <- c('')
      D <- B[B$class %in% c('delete 5\' terminal'), ]
      if(nrow(D)>0) {
        for(j in seq_len(nrow(gff))){
          D$sub_class_bulin <- (gff[j,4]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(gff[j,4]+sgmRNA_BA)
          E <- D[D$sub_class_bulin==T,]
          E <- E[,-c(ncol(E))]
          if(nrow(E)>0) {
            E$sub_class <- paste('delete leader_', gff[j,9], sep = '')
          }
          sgmRNA[[j]] <- E
          D <- D[D$sub_class_bulin ==F,]
        }
        D <- D[D$sub_class_bulin ==F,]
        D <- D[,-c(ncol(D))]
        if(nrow(D)>0) {
          unc_list <- list()
          unc_cor <- unc[[k]]
          if(length(unc_cor) != 0) {
            for(j in seq_len(length(unc_cor))) {
              D$sub_class_bulin <- (unc_cor[j]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(unc_cor[j]+sgmRNA_BA)
              E <- D[D$sub_class_bulin==T,-c(ncol(D))]
              if(nrow(E)>0) {
                E$sub_class <- paste('delete leader uncanonincal_', unc_cor[j], sep = '')
              }
              unc_list[[j]] <- E
              D <- D[D$sub_class_bulin ==F,]
            }
            D <- D[,-c(ncol(D))]
            if(nrow(D)>0) {
              D$sub_class <- c('delete 5\' UTR genome')
            }
            D <- rbind(do.call(rbind, sgmRNA), do.call(rbind, unc_list), D, A)
            output_list[[k]] <- D 
          }
          else if(length(unc_cor) == 0){
            if(nrow(D)>0) {
              D$sub_class <- c('delete 5\' UTR genome')
            }
            D <- rbind(do.call(rbind, sgmRNA),D, A)
            output_list[[k]] <- D 
          }
        }
      }
    }
    if(nrow(B[B$class %in% 'delete 5\' terminal', ])==0) {
      B$sub_class <- c('')
      B <- B[order(B$freq, decreasing = T), ]
      output_list[[k]] <- B
    }
  }
  return(output_list)
}


frag2 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 2) %>% 
  specific_reads_sgRNA_identify(gap = 50, length = 31335
                                , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                , leader = 72, UTR5 = 210, UTR3 = 31035, frag = 2)

lib_matrix <- lib_recom_matrix(lib_list = frag2, frag = 2, libs=4)

frag1 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 1)  %>%
  lin_frag1_identification(recom_matrix = lib_matrix, libs = 4
                           , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                           , UTR5 = 210, UTR3 = 31035, length = 31335)%>%
  specific_frag_recom_table(frag = 1, write = T, out_prefix = 'MHV'
                            , write_path = '~/Analysis/MHV/result/')



############################################################################################
# flanking seuqence
read_as_list <- function(path, seq='', prefix='')
{
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('dplyr', 'purrr', 'fs')
  ipak(packages)
  mini_pre.bed <- list()
  bed_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
  bed_path <- as_fs_path(bed_path)
  mini_pre.bed <- bed_path %>% map(.f = function(path){read.delim(path, header = FALSE, stringsAsFactors = F, sep = seq)})
  names(mini_pre.bed) <- bed_path
  pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
  for(i in 1:2)
  {names(mini_pre.bed) <- sub(names(mini_pre.bed), pattern = pattern[i], replacement = '')}
  return(mini_pre.bed)
  sapply(packages, require, character.only = TRUE)
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}


# function 5: format organization depending on split frag
specific_frag_format_organ <- function(list, frag)
{
  specific_fragment <- list()
  for(j in seq_along(list))
  {
    BED <- list[[j]]
    BED[,2] <- BED[,2]+1
    name <- as.data.frame(table(BED[,4]), stringsAsFactors =F)
    split <- as.data.frame(table(name$Freq), stringsAsFactors = F)
    name <- name$Var1[name$Freq %in% frag]
    if(length(name) >0 )
    {
      BED <- BED[BED[,4] %in% name,]
      BED <- BED[order(BED[,4], decreasing = F),]
      fragment <- list()
      if(frag==1) {
        B <- BED
      }
      else if(frag!=1) {
        for(i in 1:frag) {
          if(i ==1)
          {A <- BED[seq(i,nrow(BED), frag),][,1:3]}
          if(i ==frag)
          {A <- BED[seq(i,nrow(BED), frag),][,2:6]}
          if(i>1 & i<frag)
          {A <- BED[seq(i,nrow(BED), frag),][,2:3]}
          fragment[[i]] <- A
        }
        B <- as.data.frame(matrix(data = 0, nrow = nrow(fragment[[1]]), ncol = 0))
        for(i in rev(seq_along(fragment))) {
          A <- fragment[[i]]
          B <- cbind(A,B)
        }
      }
      
      colnames(B)[seq(2,(ncol(B)-3), 2)] <- paste('start', seq(1,frag,1), sep = '')
      colnames(B)[seq(3,(ncol(B)-2), 2)] <- paste('end', seq(1,frag,1), sep = '')
      specific_fragment[[j]] <- B
      names(specific_fragment)[j] <- names(list)[j]
    }
    else(next)
  }
  return(specific_fragment)
}



# Function3: identify the sgRNA (for files from Function5)
## gff_file have to be the required format
## sgmRNA_BA <- sgmRNA_boundary_allow
specific_reads_sgRNA_identify <- function(list, gap=50, gff_path, sgmRNA_BA=50
                                          , leader, UTR5, UTR3, frag, out_prefix=''
                                          , whole_genome_tolerance=20, length
                                          , isrecom=F, write=F, write_path=getwd())
{ 
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('plyr', 'dplyr')
  ipak(packages)
  
  gff <- read.csv(file = gff_path, header = F)
  
  if(frag==2)
  {
    output_list <- list()
    for(i in seq_along(list)) {
      A <- list[[i]]
      A$gap <- A[,'start2']-A[,'end1']
      B <- A
      B$small_gap <- (B[,'gap'] < gap)
      small_gap <- B[B$small_gap ==T,-c(ncol(B))]
      if(nrow(small_gap)>0) {
        small_gap$class <- c('small_gap')
        small_gap$sub_class <- c('')
      }
      if(nrow(small_gap)==0) { 
        small_gap <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      
      B <- B[B$small_gap ==F,-c(ncol(B))]
      B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
      UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
      if(nrow(UTR3_F)>0) {
        C <- UTR3_F
        C$UTR5 <- C[,'start1'] <= UTR5
        D <- C[C$UTR5 ==T,-c(ncol(C))]
        if(nrow(D)>0) {
          D$class <- c('Truncated Di')
          D$sub_class <- c('delete 3\'UTR Di')
        }
        if(nrow(D)==0) {
          D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
        E <- C[C$UTR5 ==F,-c(ncol(C))]
        if(nrow(E)>0)  {
          E$class <- c('Truncated Di')
          E$sub_class <- c('delete 3\' 5\'UTR Di')
        }
        if(nrow(E)==0) {
          E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
        delete_UTR_Di <- rbind(D,E)
      }
      
      B <- B[B$UTR3 ==T,-c(ncol(B))]
      if(nrow(B)>0) {
        C <- B
        C$UTR5 <- C[,'start1'] < UTR5
        delete_5UTR_Di <- C[C$UTR5 ==F,-c(ncol(C))]
        if(nrow(delete_5UTR_Di)>0) {
          delete_5UTR_Di$class <- c('Truncated Di')
          delete_5UTR_Di$sub_class <- c('delete 5\'UTR Di')
        }
        if(nrow(delete_5UTR_Di)==0) {
          delete_5UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
      }
      
      C <- C[C$UTR5 ==T, -c(ncol(C))]
      C$leader <- C$end1 <= leader
      sgmRNA <- C[C$leader == T, -c(ncol(C))]
      if(nrow(sgmRNA)>0) {
        sgmRNA$class <- c('sgmRNA')
        sgmRNA$sub_class <- c('')
      }
      if(nrow(sgmRNA)==0) {sgmRNA <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      
      DiRNA <- C[C$leader ==F,-c(ncol(C))]
      if(nrow(DiRNA)>0) {
        DiRNA$class <- c('DiRNA')
        DiRNA$sub_class <- c('')
      }
      if(nrow(DiRNA)==0) {
        DiRNA <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      class <- rbind(small_gap, delete_UTR_Di, delete_5UTR_Di, DiRNA, sgmRNA)
      
      B <- class
      rownames(B) <- seq(1,nrow(B))
      if(nrow(B[B$class %in% 'sgmRNA', ])>0) {
        A <- B[!(B$class %in% 'sgmRNA'),]
        sgmRNA <- list()
        D <- B[B$class %in% c('sgmRNA'), ]
        if(nrow(D)>0) {
          for(j in seq_len(nrow(gff))) {
            D$sub_class_bulin <- (gff[j,4]-sgmRNA_BA)<D[,paste('start', 2, sep = '')] & D[,paste('start', 2, sep = '')]<(gff[j,4]+sgmRNA_BA)
            E <- D[D$sub_class_bulin==T,]
            E <- E[,-c(ncol(E))]
            if(nrow(E)>0) {
              E$sub_class <- gff[j,9]
            }
            sgmRNA[[j]] <- E
            D <- D[D$sub_class_bulin ==F,]
          }
          D <- D[D$sub_class_bulin ==F,]
          D <- D[,-c(ncol(D))]
          if(nrow(D)>0) {
            D$sub_class <- c('uncanonical sgRNAs')
          }
          sgmRNA[[(length(sgmRNA)+1)]] <- D
          D <- do.call(rbind, sgmRNA)
          B <- rbind(A,D)
          if(isTRUE(isrecom)) {
            B <- B[order(B$freq, decreasing = T),]
          }
        }
      }
      output_list[[i]] <- B
      names(output_list)[[i]] <- names(list)[[i]]
    }
  }
  
  if(frag>2) {
    output_list <- list()
    for(i in seq_along(list)) {
      frag_other <- list[[i]]
      if(nrow(frag_other)>0) {
        B <- frag_other
        dfgap <- as.data.frame(matrix(data = '', nrow = nrow(B), ncol = 1))
        for(j in 2:frag) {
          dfgap[,(j-1)] <- frag_other[,paste('end', j, sep = '')]-frag_other[,paste('start', (j-1), sep = '')]
        }
        dfgap$min <- apply(dfgap, MARGIN = 1,FUN = min)
        B$gap <- dfgap$min
        B$small_gap <- B[,'gap'] <gap
        small_gap <- B[B$small_gap ==T,-c(ncol(B))]
        if(nrow(small_gap)>0) {
          small_gap$class <- c('small_gap')
          small_gap$sub_class <- c('')
        }
        if(nrow(small_gap)==0) {
          small_gap <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        
        B <- B[B$small_gap ==F,-c(ncol(B))]
        B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
        UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
        if(nrow(UTR3_F)>0) {
          C <- UTR3_F
          C$UTR5 <- C[,'start1'] <= UTR5
          D <- C[C$UTR5 ==T,-c(ncol(C))]
          if(nrow(D)>0) {
            D$class <- c('Truncated Di')
            D$sub_class <- c('delete 3\'UTR Di')
          }
          if(nrow(D)==0) {
            D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          E <- C[C$UTR5 ==F,-c(ncol(C))]
          if(nrow(E)>0) {
            E$class <- c('Truncated Di')
            E$sub_class <- c('delete 3\' 5\'UTR Di')
          }
          if(nrow(E)==0) {
            E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          delete_UTR_Di <- rbind(D,E)
        }
        if(nrow(UTR3_F)==0) {
          delete_UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        
        UTR3_T <- B[B$UTR3 ==T,-c(ncol(B))]
        if(nrow(UTR3_T)>0) {
          C <- UTR3_T
          C$UTR5 <- C[,'start1'] > UTR5
          D <- C[C$UTR5 ==T,-c(ncol(C))]
          if(nrow(D)>0) {
            D$class <- c('Truncated Di')
            D$sub_class <- c('delete 5\'UTR Di')
          }
          if(nrow(D)==0) { 
            D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          E <- C[C$UTR5 ==F,-c(ncol(C))]
          if(nrow(E)>0) {
            E$class <- c('DiRNA')
            E$sub_class <- c('')
          }
          if(nrow(E)==0) {
            E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          delete_UTR5_Di <- rbind(D,E)
        }
        if(nrow(UTR3_T)==0) {
          delete_UTR5_Di <-  as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        class <- rbind(delete_UTR5_Di, delete_UTR_Di, small_gap)
        if(isTRUE(isrecom)) {class <- class[order(class$freq, decreasing = T),]
        }
      }
      output_list[[i]] <- class
      names(output_list)[[i]] <- names(list)[[i]]
    }
  }
  
  if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False')) {
    return(output_list)
  }
  if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true')) {
    if(isFALSE(grepl(write_path, pattern = '&/'))){
      write_path <- paste(write_path,'/', sep = '')
    }
    lapply(seq_along(output_list)
           , function(i) 
             write.csv(output_list[[i]]
                       , file = paste0(write_path
                                       ,names(output_list[i]),'_', frag, 'fragment_', out_prefix, ".csv"),
                       row.names = FALSE, quote = F))
  }
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}




# function6 : construct recombination point sort and new seq name
specific_frag_recom_freq <- function(list, frag) {
  recom <- list()
  for(i in seq_along(list))
  {
    B <- list[[i]]
    B[B$sub_class %in% '','sub_class'] <- '.'
    colNames <-  colnames(B)[c(2:(frag*2+1),(frag*2+6), (frag*2+7))] # could be any number of column names here
    B$recom <- do.call(paste, c(B[colNames], sep="-"))
    recom_table <- as.data.frame(table(B$recom), stringsAsFactors = F)
    recom_table <- recom_table[order(recom_table$Freq, decreasing = T),]
    frequency <- recom_table$Freq
    out <- as.data.frame(do.call(rbind, strsplit(recom_table$Var1, split = '-')))
    out$name <- recom_table$Var1
    out$freq <- recom_table$Freq
    
    recom_table <- as.data.frame(matrix('U00735.2', nrow = nrow(out), ncol = 1))
    recom_table[,2:(1+(frag*4))] <- out[,1:ncol(out)]
    colnames(recom_table)[seq(2, (1+(frag*2)), 2)] <- paste('start', seq(1,frag), sep = '')
    colnames(recom_table)[seq(3, (1+(frag*2)), 2)] <- paste('end', seq(1,frag), sep ='')
    recom_table$score <- 0
    recom_table$strand <- c('+')
    for(j in seq(2, (frag*2+1),1))
    {recom_table[,j] <- as.numeric(recom_table[,j])}
    D <- 0
    for(k in seq(3, (frag*2+1),2))
    {D <- recom_table[,k]-recom_table[,(k-1)]+(frag-1)+D}
    recom_table$length <- D
    recom_table$gap <- recom_table[,frag*2]-recom_table[,3]
    recom_table$frag <- frag
    recom[[i]] <- recom_table
    names(recom)[i] <- names(list)[i]
  }
  return(recom)
}



# function6 : construct recombination point sort and new seq name
one_frag_recom_freq <- function(list) {
  recom <- list()
  for(i in seq_along(list))
  {
    B <- list[[i]]
    B[B$sub_class %in% '','sub_class'] <- '.'
    colNames <-  colnames(B)[c(2:(1*2+1),(1*2+6), (1*2+7))] # could be any number of column names here
    B$recom <- do.call(paste, c(B[colNames], sep="-"))
    recom_table <- as.data.frame(table(B$recom), stringsAsFactors = F)
    recom_table <- recom_table[order(recom_table$Freq, decreasing = T),]
    frequency <- recom_table$Freq
    out <- as.data.frame(do.call(rbind, strsplit(recom_table$Var1, split = '-')))
    out$name <- recom_table$Var1
    out$freq <- recom_table$Freq
    
    recom_table <- as.data.frame(matrix('U00735.2', nrow = nrow(out), ncol = 1))
    recom_table[,2:(1+6)] <- out[,1:ncol(out)]
    colnames(recom_table)[seq(2, (1+(1*2)), 2)] <- paste('start', 1, sep = '')
    colnames(recom_table)[seq(3, (1+(1*2)), 2)] <- paste('end', 1, sep ='')
    recom_table$score <- 0
    recom_table$strand <- c('+')
    for(j in seq(2, (1*2+1),1))
    {recom_table[,j] <- as.numeric(recom_table[,j])}
    D <- 0
    recom_table$length <- recom_table$end1- recom_table$start1+1
    recom_table$frag <- 1
    recom_table <- recom_table[order(recom_table$freq, decreasing = T), ]
    recom[[i]] <- recom_table
    names(recom)[i] <- names(list)[i]
  }
  return(recom)
}




# Function 7: get flanking sequence based on input frags and format from function 5
get_flanking <- function(list, frag, flanking=20, genome_length
                         , writeto_bed=FALSE, write_path, chr) {
  list.getfasta<- list()
  for(i in seq_along(list)) {
    A <- list[[i]]
    for(k in seq(2,(1+frag*2))) {
      A[,k] <- as.numeric(A[,k])
    }
    for(j in seq(2,(frag*2+1))) {
      B <- as.data.frame(A[,c(j, (frag*2+4))])
      B$class <- B[,1] <= flanking
      C <- as.data.frame(B[B$class ==T,][,c(1,2)])
      if(nrow(C)>0) {
        D <- as.data.frame(matrix(data = 0, nrow = nrow(C), ncol = 1))
        D$UTR3 <- C[,1]+flanking
        D$name <- C[,2]
        C <- D
      }
      B <- as.data.frame(B[B$class ==F, ][,1:2])
      B$class <- B[,1] >= genome_length-flanking
      D <- as.data.frame(B[B$class ==T,][,1:2])
      if(nrow(D)>0)  {
        E <- as.data.frame(matrix(data = 0, nrow = nrow(D), ncol = 1))
        E[,1] <- D[,1]-flanking
        E$UTR3 <- genome_length
        E$name <- D[,2]
        D <- E
      }
      B <- as.data.frame(B[B$class ==F,][,1:2])
      E <- as.data.frame(matrix(data = 0, nrow = nrow(B), ncol = 1))
      E[,1] <- B[,1]-flanking
      E$UTR3 <- B[,1]+flanking
      E$name <- B[,2]
      B <- rbind(C,D,E)
      G <- as.data.frame(matrix(data = chr, nrow = nrow(B), ncol = 1))
      B <- cbind(G,B)
      B$score <- 0
      B$strand <- c('+')
      list.getfasta[[(2*frag*i)-(2*frag-(j-1))]] <- B
      names(list.getfasta)[(2*frag*i)-(2*frag-(j-1))] <- paste(names(list)[i], colnames(A)[j], sep = '_')
    }
  }
  if(writeto_bed== FALSE| writeto_bed==c('false')| writeto_bed==c('F')| writeto_bed==c('f')| writeto_bed==c('False'))
  {return(list.getfasta)}
  else if(writeto_bed==TRUE | writeto_bed==c('True')| writeto_bed==c('T')| writeto_bed==c('t')| writeto_bed==c('true'))
  {lapply(seq_along(list.getfasta)
          , function(i) write.table(list.getfasta[[i]], append = F
                                    , file = paste0(write_path,names(list.getfasta[i])
                                                    ,'_','frag', frag,".bed"),
                                    row.names = FALSE, sep = '\t', col.names = F, quote = F))}
  else(print('Anstwer TRUE(T) and FALSE(F) only'))
}



# function 8: flanking sequence fasta file organization
org_flanking_fasta <- function(path, seq='\n', prefix='') {
  read_as_list <- function(path, seq='', prefix=''){
    ipak <- function(pkg){
      new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
      if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
      sapply(pkg, require, character.only = TRUE)
    }
    packages <- c('dplyr', 'purrr', 'fs')
    ipak(packages)
    mini_pre.bed <- list()
    bed_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
    bed_path <- as_fs_path(bed_path)
    mini_pre.bed <- bed_path %>% map(.f = function(path){read.delim(path, header = FALSE, stringsAsFactors = F, sep = seq)})
    names(mini_pre.bed) <- bed_path
    pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
    for(i in 1:2)
    {names(mini_pre.bed) <- sub(names(mini_pre.bed), pattern = pattern[i], replacement = '')}
    return(mini_pre.bed)
    sapply(packages, require, character.only = TRUE)
    lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
  }
  
  flanking_fasta <- read_as_list(path = path, seq=seq, prefix=prefix)
  for(i in seq_along(flanking_fasta)) {
    fasta <- flanking_fasta[[i]]
    fasta1 <- as.data.frame(fasta[seq(1,nrow(fasta),2),])
    fasta1$seq <- fasta[seq(2,nrow(fasta),2),]
    flanking_fasta[[i]] <- fasta1
  }
  return(flanking_fasta)
}


# function 9: final output file (recombination point with flanking region)
final_output_generation <- function(list_recom, list_flanking_fa, frag, write=FALSE, write_path)
{
  out_list <- list()
  for(i in seq_along(list_recom))
  {
    recom <- list_recom[[i]]
    recom <- recom[order(recom$name),]
    out <- as.data.frame(matrix(data = 0, ncol = (frag*2*2+4), nrow(recom)))
    out[,1] <- frag
    out[,seq(2, frag*2*2, 2)] <- recom[,c(seq(2,(frag*2+1)))]
    colnames(out)[seq(2, frag*2*2, 4)] <- paste('start', seq(1,frag), sep = '')
    colnames(out)[seq(4, frag*2*2, 4)] <- paste('end', seq(1,frag), sep = '')
    subset_flanking <- grep(names(list_flanking_fa), pattern = names(list_recom[i]), value = T)
    subset_flanking <- list_flanking_fa[subset_flanking]
    # start flanking
    for(j in seq(3, (frag*2*2), 4)) {
      flanking <- subset_flanking[[(j+1)/4+frag]]
      print((j+1)/4+frag)
      flanking <- flanking[order(flanking[,1]),]
      out[,j] <- flanking[,2]
    }
    colnames(out)[seq(3, frag*2*2, 4)] <- paste('start', seq(1,frag), '_flanking', sep = '')
    # end flanking
    for(k in seq(5, (frag*2*2+1), 4))  {
      flanking <- subset_flanking[[(k-1)/4]]
      print((k-1)/4)
      flanking <- flanking[order(flanking[,1]),]
      out[,k] <- flanking[,2]
    }
    colnames(out)[seq(5, (frag*2*2+1), 4)] <- paste('end', seq(1,frag), '_flanking', sep = '')
    out[,(frag*2*2+2)] <- recom$freq
    out[,(frag*2*2+3)] <- recom$length
    out[,(frag*2*2+4)] <- recom[,(frag*2+2)]
    out[,(frag*2*2+5)] <- recom[,(frag*2+3)]
    colnames(out)[c(1,(frag*2*2+2),(frag*2*2+3),(frag*2*2+4),(frag*2*2+5))] <- 
      c('frag', 'freq', 'length', 'class', 'sub_class')
    out <- out[order(out[,(frag*2*2+2)], decreasing = T),]
    out_list[[i]] <- out
    names(out_list)[i] <- names(list_recom)[i]
  }
  if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False'))
  {return(out_list)}
  else if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true'))
  {lapply(seq_along(out_list), function(i) write.csv(out_list[[i]],
                                                     file = paste0(write_path,
                                                                   names(out_list[i]),'_', 'frag', frag, ".csv"),
                                                     row.names = FALSE, quote = F))}
  
}


## frag2-frag5
library(dplyr)
#### write.flanking fa
bed <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/')
for(j in seq(2,5,1)) {
  write_path <- paste('~/Analysis/MHV/flanking/MHV_frag', j, '/',sep = '')
  dir.create(write_path)
  bed%>%
    specific_frag_format_organ(frag = j) %>% 
    specific_reads_sgRNA_identify(gap = 50, length = 31335
                                  , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                  , leader = 72, UTR5 = 210, UTR3 = 31035, frag = j) %>%
    specific_frag_recom_freq(frag = j) %>%
    get_flanking(frag = j, flanking = 20, 
                 genome_length = 31335, writeto_bed = T
                 , write_path = write_path
                 , chr = 'NC_001846.1')
}


#### generate the final output flanking sequence
bed <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/')
for(j in seq(2,5,1)) {
  frag <- bed%>%
    specific_frag_format_organ(frag = j) %>% 
    specific_reads_sgRNA_identify(gap = 50, length = 31335
                                  , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                  , leader = 72, UTR5 = 210, UTR3 = 31035, frag = j) %>%
    specific_frag_recom_freq(frag = j)
  
  
  frag_flanking <- org_flanking_fasta(path = paste('~/Analysis/MHV/flanking/MHV_frag', j, '/', sep='')
                                      , prefix = '.fa')
  
  final_output_generation(list_recom = frag, frag = j
                          , list_flanking_fa = frag_flanking
                          , write = T, write_path = '~/Analysis/MHV/result/flanking_sequence/')
}



## frag1 flanking seuqence
lib_recom_matrix <- function(lib_list, frag, libs) {
  recom_table_list <- list()
  recom_list <- list()
  inter_count <- ''
  x=0
  for(i in seq_along(lib_list)) {
    lib1 <- lib_list[[i]]
    lib1 <- lib1[!(lib1$class %in% 'small_gap'), ]
    lib1[rownames(lib1[lib1$sub_class %in% '',]),'sub_class'] <- ' '
    if(frag !=1) {
      colNames <-  colnames(lib1)[c(2:(frag*2+1),(frag*2+6), (frag*2+7))]
    }
    if(frag ==1) {
      colNames <-  colnames(lib1)[c(2,3,8,9)]
    }
    lib1$recom <-  do.call(paste, c(lib1[colNames], sep="-"))
    lib1 <- as.data.frame(table(lib1$recom), stringsAsFactors = F)
    recom_list[[i]] <- lib1$Var1
    recom_table_list[[i]] <- lib1
  }
  names(recom_list) <- names(lib_list)
  names(recom_table_list) <- names(lib_list)
  
  all_recom <- Reduce(union, recom_list)
  lib_matrix <- as.data.frame(matrix(data = 0, nrow = length(all_recom), ncol = libs))
  lib_matrix$index <- seq(1, nrow(lib_matrix))
  rownames(lib_matrix) <- all_recom; colnames(lib_matrix) <- names(lib_list)
  
  for(i in seq_along(recom_table_list)) {
    B <- recom_table_list[[i]]
    B <- B[order(B$Var1, decreasing = T), ]
    index <- lib_matrix[rownames(lib_matrix) %in% B$Var1, (ncol(lib_matrix))]
    C <- lib_matrix[lib_matrix[,ncol(lib_matrix)] %in% index,]
    C <- C[order(rownames(C), decreasing = T), ]
    lib_matrix <- lib_matrix[!(lib_matrix[,ncol(lib_matrix)] %in% index), ]
    print(table(rownames(C)== B$Var1))
    C[,i] <- B$Freq
    lib_matrix <- rbind(lib_matrix,C)
  }
  lib_matrix <- lib_matrix[, -c(ncol(lib_matrix))]
}


lin_frag1_identification <- function(list, recom_matrix, gff_path, libs
                                     , UTR5, UTR3, length, sgmRNA_BA=50
                                     , whole_genome_tolerance=20) {
  unc <- list()
  for(j in seq(1,libs,2)) {
    B <- recom_matrix[,j:(j+1)]
    C <- as.data.frame(ifelse(test = B>=2, yes = 1, no = 0))
    C$con <- rowSums(as.matrix(C[,1:2]))
    C <- C[C$con > 0, ]
    C <- C[grep(rownames(C), pattern = 'uncanonical'), ]
    C <- C[C$con >1, ]
    C <- as.data.frame(do.call(rbind, strsplit(rownames(C), split = '-')))
    unc[[j]] <- rep(as.numeric(unique(C$V3)))
    unc[[j+1]] <- rep(as.numeric(unique(C$V3)))
  }
  
  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('plyr', 'dplyr')
  ipak(packages)
  
  gff <- read.csv(file = gff_path, header = F)
  
  output_list <- list()
  for(i in seq_along(list))  {
    f1 <- list[[i]]
    f1$length <- f1[,3]-f1[,2]
    f1$complete <- f1$length > (length-whole_genome_tolerance)
    complete <- f1[f1$complete %in% TRUE,-c(ncol(f1))]
    if(nrow(complete)>0){
      complete$class <- c('complete genome')
    }
    if(nrow(complete)==0){
      complete <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- f1[f1$complete %in% FALSE,-c(ncol(f1))]
    f1$utr3 <- f1$end> UTR3
    delete_utr3 <- f1[f1$utr3 %in% FALSE,-c(ncol(f1))]
    if(nrow(delete_utr3)>0){
      delete_utr3$utr5 <- delete_utr3$start < UTR5
      delete_utr3_utr5 <- delete_utr3[delete_utr3$utr5 %in% FALSE, -c(ncol(delete_utr3))]
      if(nrow(delete_utr3_utr5)>0){
        delete_utr3_utr5$class <- c('delete 3\'UTR and 5\' UTR genome')
      }
      if(nrow(delete_utr3_utr5)==0)
      {delete_utr3_utr5 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
      
      delete_utr3 <- delete_utr3[delete_utr3$utr5 %in% TRUE, -c(ncol(delete_utr3))]
      if(nrow(delete_utr3)>0){
        delete_utr3$class <- c('delete 3\'UTR genome')
      }
      if(nrow(delete_utr3)==0)
      {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    }
    if(nrow(delete_utr3)==0)
    {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- f1[f1$utr3 %in% TRUE,-c(ncol(f1))]
    f1$utr5 <- f1$start < UTR5
    delete_utr5 <- f1[f1$utr5 %in% FALSE,-c(ncol(f1))]
    if(nrow(delete_utr5)>0)
    {
      delete_utr5$class <- c('delete 5\' terminal')
    }
    if(nrow(delete_utr5)==0)
    {delete_utr5 <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    near_complete <- f1[f1$utr5 %in% T,-c(ncol(f1))]
    if(nrow(near_complete)>0)
    {
      near_complete$class <- c('near complete genome')
    }
    if(nrow(near_complete)==0)
    {near_complete <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- rbind(delete_utr3, delete_utr3_utr5, delete_utr5, near_complete)
    output_list[[i]] <- f1
    names(output_list)[[i]] <- names(list)[[i]]
  }
  
  for(k in 1:length(output_list)) {
    B <- output_list[[k]]
    rownames(B) <- seq(1,nrow(B))
    if(nrow(B[B$class %in% 'delete 5\' terminal', ])>0) {
      A <- B[!(B$class %in% 'delete 5\' terminal'),]
      sgmRNA <- list()
      A$sub_class <- c('')
      D <- B[B$class %in% c('delete 5\' terminal'), ]
      if(nrow(D)>0) {
        for(j in seq_len(nrow(gff))){
          D$sub_class_bulin <- (gff[j,4]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(gff[j,4]+sgmRNA_BA)
          E <- D[D$sub_class_bulin==T,]
          E <- E[,-c(ncol(E))]
          if(nrow(E)>0) {
            E$sub_class <- paste('delete leader_', gff[j,9], sep = '')
          }
          sgmRNA[[j]] <- E
          D <- D[D$sub_class_bulin ==F,]
        }
        D <- D[D$sub_class_bulin ==F,]
        D <- D[,-c(ncol(D))]
        if(nrow(D)>0) {
          unc_list <- list()
          unc_cor <- unc[[k]]
          if(length(unc_cor) != 0) {
            for(j in seq_len(length(unc_cor))) {
              D$sub_class_bulin <- (unc_cor[j]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(unc_cor[j]+sgmRNA_BA)
              E <- D[D$sub_class_bulin==T,-c(ncol(D))]
              if(nrow(E)>0) {
                E$sub_class <- paste('delete leader uncanonincal_', unc_cor[j], sep = '')
              }
              unc_list[[j]] <- E
              D <- D[D$sub_class_bulin ==F,]
            }
            D <- D[,-c(ncol(D))]
            if(nrow(D)>0) {
              D$sub_class <- c('delete 5\' UTR genome')
            }
            D <- rbind(do.call(rbind, sgmRNA), do.call(rbind, unc_list), D, A)
            output_list[[k]] <- D 
          }
          else if(length(unc_cor) == 0){
            if(nrow(D)>0) {
              D$sub_class <- c('delete 5\' UTR genome')
            }
            D <- rbind(do.call(rbind, sgmRNA),D, A)
            output_list[[k]] <- D 
          }
        }
      }
    }
    if(nrow(B[B$class %in% 'delete 5\' terminal', ])==0) {
      B$sub_class <- c('')
      B <- B[order(B$freq, decreasing = T), ]
      output_list[[k]] <- B
    }
  }
  return(output_list)
}

frag2 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 2) %>% 
  specific_reads_sgRNA_identify(gap = 50, length = 31335
                                , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                , leader = 72, UTR5 = 210, UTR3 = 31035, frag = 2)

lib_matrix <- lib_recom_matrix(lib_list = frag2, frag = 2, libs=4)

frag1 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 1)  %>%
  lin_frag1_identification(recom_matrix = lib_matrix, libs = 4
                           , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                           , UTR5 = 210, UTR3 = 31035, length = 31335) %>%
  one_frag_recom_freq()

write_path <- paste('~/Analysis/MHV/flanking/MHV_frag', 1, '/',sep = '')
dir.create(write_path)

get_flanking(list = frag1, frag = 1, flanking = 20, 
             genome_length = 31335, writeto_bed = T
             , write_path = write_path
             , chr = 'NC_001846.1')


# function 9: final output file (recombination point with flanking region)
oneF_final_output_generation <- function(list_recom, list_flanking_fa, write=FALSE, write_path)
{
  out_list <- list()
  for(i in seq_along(list_recom))
  {
    recom <- list_recom[[i]]
    recom <- recom[order(recom$name),]
    out <- as.data.frame(matrix(data = 0, ncol = 8, nrow(recom)))
    out[,1] <- 1
    out[,c(2, 4)] <- recom[,c(2,3)]
    colnames(out)[2] <- paste('start',1 , sep = '')
    colnames(out)[4] <- paste('end', 1, sep = '')
    # start flanking
    flanking <- list_flanking_fa[[2*i]]
    flanking <- flanking[order(flanking[,1]),]
    out[,3] <- flanking[,2]
    colnames(out)[3] <- paste('start', 1, '_flanking', sep = '')
    flanking <- list_flanking_fa[[2*i-1]]
    flanking <- flanking[order(flanking[,1]),]
    out[,5] <- flanking[,2]
    colnames(out)[5] <- paste('end', 1, '_flanking', sep = '')
    out[,c(6,7,8,9)] <- recom[,c(7,10,4,5)]
    colnames(out)[c(6,7,8,9)] <- c('freq', 'length', 'class', 'sub_class')
    out <- out[order(out[,'freq'], decreasing = T),]
    out_list[[i]] <- out
    names(out_list)[i] <- names(list_recom)[i]
  }
  if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False'))
  {return(out_list)}
  else if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true'))
  {lapply(seq_along(out_list), function(i) write.csv(out_list[[i]],
                                                     file = paste0(write_path,
                                                                   names(out_list[i]),'_', 'frag1', ".csv"),
                                                     row.names = FALSE, quote = F))}
  
}


frag_flanking <- org_flanking_fasta(path = '~/Analysis/MHV/flanking/MHV_frag1/', prefix = '.fa')

out_list <- oneF_final_output_generation(list_recom = frag1
                                         , list_flanking_fa = frag_flanking
                                         , write = T
                                         , write_path = '~/Analysis/MHV/result/f1_0515/')


############################################################################################
## reproduction
read_as_list <- function(path, seq='', prefix='')
{
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('dplyr', 'purrr', 'fs')
  ipak(packages)
  mini_pre.bed <- list()
  bed_path <- paste(path ,list.files(path = path, pattern = prefix), sep = '')
  bed_path <- as_fs_path(bed_path)
  mini_pre.bed <- bed_path %>% map(.f = function(path){read.delim(path, header = FALSE, stringsAsFactors = F, sep = seq)})
  names(mini_pre.bed) <- bed_path
  pattern <- c('.*\\/', '\\.[a-z\\.]+', '\\.[a-z]+')
  for(i in 1:2)
  {names(mini_pre.bed) <- sub(names(mini_pre.bed), pattern = pattern[i], replacement = '')}
  return(mini_pre.bed)
  sapply(packages, require, character.only = TRUE)
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}


# function 5: format organization depending on split frag
specific_frag_format_organ <- function(list, frag)
{
  specific_fragment <- list()
  for(j in seq_along(list))
  {
    BED <- list[[j]]
    BED[,2] <- BED[,2]+1
    name <- as.data.frame(table(BED[,4]), stringsAsFactors =F)
    split <- as.data.frame(table(name$Freq), stringsAsFactors = F)
    name <- name$Var1[name$Freq %in% frag]
    if(length(name) >0 )
    {
      BED <- BED[BED[,4] %in% name,]
      BED <- BED[order(BED[,4], decreasing = F),]
      fragment <- list()
      if(frag==1) {
        B <- BED
      }
      else if(frag!=1) {
        for(i in 1:frag) {
          if(i ==1)
          {A <- BED[seq(i,nrow(BED), frag),][,1:3]}
          if(i ==frag)
          {A <- BED[seq(i,nrow(BED), frag),][,2:6]}
          if(i>1 & i<frag)
          {A <- BED[seq(i,nrow(BED), frag),][,2:3]}
          fragment[[i]] <- A
        }
        B <- as.data.frame(matrix(data = 0, nrow = nrow(fragment[[1]]), ncol = 0))
        for(i in rev(seq_along(fragment))) {
          A <- fragment[[i]]
          B <- cbind(A,B)
        }
      }
      
      colnames(B)[seq(2,(ncol(B)-3), 2)] <- paste('start', seq(1,frag,1), sep = '')
      colnames(B)[seq(3,(ncol(B)-2), 2)] <- paste('end', seq(1,frag,1), sep = '')
      specific_fragment[[j]] <- B
      names(specific_fragment)[j] <- names(list)[j]
    }
    else(next)
  }
  return(specific_fragment)
}



# Function3: identify the sgRNA (for files from Function5)
## gff_file have to be the required format
## sgmRNA_BA <- sgmRNA_boundary_allow
specific_reads_sgRNA_identify <- function(list, gap=50, gff_path, sgmRNA_BA=50
                                          , leader, UTR5, UTR3, frag, out_prefix=''
                                          , whole_genome_tolerance=20, length
                                          , isrecom=F, write=F, write_path=getwd())
{ 
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('plyr', 'dplyr')
  ipak(packages)
  
  gff <- read.csv(file = gff_path, header = F)
  
  if(frag==2)
  {
    output_list <- list()
    for(i in seq_along(list)) {
      A <- list[[i]]
      A$gap <- A[,'start2']-A[,'end1']
      B <- A
      B$small_gap <- (B[,'gap'] < gap)
      small_gap <- B[B$small_gap ==T,-c(ncol(B))]
      if(nrow(small_gap)>0) {
        small_gap$class <- c('small_gap')
        small_gap$sub_class <- c('')
      }
      if(nrow(small_gap)==0) { 
        small_gap <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      
      B <- B[B$small_gap ==F,-c(ncol(B))]
      B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
      UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
      if(nrow(UTR3_F)>0) {
        C <- UTR3_F
        C$UTR5 <- C[,'start1'] <= UTR5
        D <- C[C$UTR5 ==T,-c(ncol(C))]
        if(nrow(D)>0) {
          D$class <- c('Truncated Di')
          D$sub_class <- c('delete 3\'UTR Di')
        }
        if(nrow(D)==0) {
          D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
        E <- C[C$UTR5 ==F,-c(ncol(C))]
        if(nrow(E)>0)  {
          E$class <- c('Truncated Di')
          E$sub_class <- c('delete 3\' 5\'UTR Di')
        }
        if(nrow(E)==0) {
          E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
        delete_UTR_Di <- rbind(D,E)
      }
      
      B <- B[B$UTR3 ==T,-c(ncol(B))]
      if(nrow(B)>0) {
        C <- B
        C$UTR5 <- C[,'start1'] < UTR5
        delete_5UTR_Di <- C[C$UTR5 ==F,-c(ncol(C))]
        if(nrow(delete_5UTR_Di)>0) {
          delete_5UTR_Di$class <- c('Truncated Di')
          delete_5UTR_Di$sub_class <- c('delete 5\'UTR Di')
        }
        if(nrow(delete_5UTR_Di)==0) {
          delete_5UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
        }
      }
      
      C <- C[C$UTR5 ==T, -c(ncol(C))]
      C$leader <- C$end1 <= leader
      sgmRNA <- C[C$leader == T, -c(ncol(C))]
      if(nrow(sgmRNA)>0) {
        sgmRNA$class <- c('sgmRNA')
        sgmRNA$sub_class <- c('')
      }
      if(nrow(sgmRNA)==0) {sgmRNA <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      
      DiRNA <- C[C$leader ==F,-c(ncol(C))]
      if(nrow(DiRNA)>0) {
        DiRNA$class <- c('DiRNA')
        DiRNA$sub_class <- c('')
      }
      if(nrow(DiRNA)==0) {
        DiRNA <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(A)+2)))
      }
      class <- rbind(small_gap, delete_UTR_Di, delete_5UTR_Di, DiRNA, sgmRNA)
      
      B <- class
      rownames(B) <- seq(1,nrow(B))
      if(nrow(B[B$class %in% 'sgmRNA', ])>0) {
        A <- B[!(B$class %in% 'sgmRNA'),]
        sgmRNA <- list()
        D <- B[B$class %in% c('sgmRNA'), ]
        if(nrow(D)>0) {
          for(j in seq_len(nrow(gff))) {
            D$sub_class_bulin <- (gff[j,4]-sgmRNA_BA)<D[,paste('start', 2, sep = '')] & D[,paste('start', 2, sep = '')]<(gff[j,4]+sgmRNA_BA)
            E <- D[D$sub_class_bulin==T,]
            E <- E[,-c(ncol(E))]
            if(nrow(E)>0) {
              E$sub_class <- gff[j,9]
            }
            sgmRNA[[j]] <- E
            D <- D[D$sub_class_bulin ==F,]
          }
          D <- D[D$sub_class_bulin ==F,]
          D <- D[,-c(ncol(D))]
          if(nrow(D)>0) {
            D$sub_class <- c('uncanonical sgRNAs')
          }
          sgmRNA[[(length(sgmRNA)+1)]] <- D
          D <- do.call(rbind, sgmRNA)
          B <- rbind(A,D)
          if(isTRUE(isrecom)) {
            B <- B[order(B$freq, decreasing = T),]
          }
        }
      }
      output_list[[i]] <- B
      names(output_list)[[i]] <- names(list)[[i]]
    }
  }
  
  if(frag>2) {
    output_list <- list()
    for(i in seq_along(list)) {
      frag_other <- list[[i]]
      if(nrow(frag_other)>0) {
        B <- frag_other
        dfgap <- as.data.frame(matrix(data = '', nrow = nrow(B), ncol = 1))
        for(j in 2:frag) {
          dfgap[,(j-1)] <- frag_other[,paste('end', j, sep = '')]-frag_other[,paste('start', (j-1), sep = '')]
        }
        dfgap$min <- apply(dfgap, MARGIN = 1,FUN = min)
        B$gap <- dfgap$min
        B$small_gap <- B[,'gap'] <gap
        small_gap <- B[B$small_gap ==T,-c(ncol(B))]
        if(nrow(small_gap)>0) {
          small_gap$class <- c('small_gap')
          small_gap$sub_class <- c('')
        }
        if(nrow(small_gap)==0) {
          small_gap <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        
        B <- B[B$small_gap ==F,-c(ncol(B))]
        B$UTR3 <-  B[, paste('end', frag, sep = '')] > UTR3
        UTR3_F <- B[B$UTR3 ==F,-c(ncol(B))]
        if(nrow(UTR3_F)>0) {
          C <- UTR3_F
          C$UTR5 <- C[,'start1'] <= UTR5
          D <- C[C$UTR5 ==T,-c(ncol(C))]
          if(nrow(D)>0) {
            D$class <- c('Truncated Di')
            D$sub_class <- c('delete 3\'UTR Di')
          }
          if(nrow(D)==0) {
            D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          E <- C[C$UTR5 ==F,-c(ncol(C))]
          if(nrow(E)>0) {
            E$class <- c('Truncated Di')
            E$sub_class <- c('delete 3\' 5\'UTR Di')
          }
          if(nrow(E)==0) {
            E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          delete_UTR_Di <- rbind(D,E)
        }
        if(nrow(UTR3_F)==0) {
          delete_UTR_Di <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        
        UTR3_T <- B[B$UTR3 ==T,-c(ncol(B))]
        if(nrow(UTR3_T)>0) {
          C <- UTR3_T
          C$UTR5 <- C[,'start1'] > UTR5
          D <- C[C$UTR5 ==T,-c(ncol(C))]
          if(nrow(D)>0) {
            D$class <- c('Truncated Di')
            D$sub_class <- c('delete 5\'UTR Di')
          }
          if(nrow(D)==0) { 
            D <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          E <- C[C$UTR5 ==F,-c(ncol(C))]
          if(nrow(E)>0) {
            E$class <- c('DiRNA')
            E$sub_class <- c('')
          }
          if(nrow(E)==0) {
            E <- as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
          }
          delete_UTR5_Di <- rbind(D,E)
        }
        if(nrow(UTR3_T)==0) {
          delete_UTR5_Di <-  as.data.frame(matrix(data = 0, nrow = 0, ncol = (ncol(frag_other)+2)))
        }
        class <- rbind(delete_UTR5_Di, delete_UTR_Di, small_gap)
        if(isTRUE(isrecom)) {class <- class[order(class$freq, decreasing = T),]
        }
      }
      output_list[[i]] <- class
      names(output_list)[[i]] <- names(list)[[i]]
    }
  }
  
  if(write== FALSE| write==c('false')| write==c('F')| write==c('f')| write==c('False')) {
    return(output_list)
  }
  if(write==TRUE | write==c('True')| write==c('T')| write==c('t')| write==c('true')) {
    if(isFALSE(grepl(write_path, pattern = '&/'))){
      write_path <- paste(write_path,'/', sep = '')
    }
    lapply(seq_along(output_list)
           , function(i) 
             write.csv(output_list[[i]]
                       , file = paste0(write_path
                                       ,names(output_list[i]),'_', frag, 'fragment_', out_prefix, ".csv"),
                       row.names = FALSE, quote = F))
  }
  lapply(paste('package:', packages, sep = ''), detach, character.only = TRUE)
}


cut_by_n_bp <- function(list, genome_length, bin=10, frag) {
  index <- as.data.frame(matrix(data = 0, nrow = length(seq(0, genome_length, 10)), ncol = 0))
  index[,1] <- seq(0, genome_length, 10)
  index[nrow(index),1] <- genome_length
  out_list <- list()
  for(j in seq_along(list)) {
    A <- list[[j]]
    A$index <- seq(1, nrow(A))
    colN <- c(paste('start', seq_len(frag), sep = ''), paste('end', seq_len(frag), sep = ''))
    for(k in colN) {
      C <- as.data.frame(matrix(data = '', nrow = 0, ncol = 0))
      for(i in 2:nrow(index)) {
        A[,k] <- as.numeric(A[,k])
        A[,'clu'] <- A[,k] > index[i-1,1] & A[,k] <= index[i,1]
        row <- A$index[A$clu %in% T]
        if(length(row)>0) {
          B <- A[A$index %in% row,-c(ncol(A))]
          B[rownames(A[A$clu %in% T,]),k] <- index[i,1]
          C <- rbind(B,C)
          A <- A[!(A$index %in% row),-c(ncol(A))]
        }
      }
      A <- C
    }
    
    out_list[[j]] <- A
    names(out_list)[j] <- names(list)[j]
  }
  return(out_list)
}

library(ggplot2)
gg_theme= theme(
  axis.title.x = element_text(size = 22),
  axis.text.x = element_text(size = 18),
  axis.title.y = element_text(size = 22),
  axis.text.y = element_text(size = 18),
  title = element_text(size = 26),
  plot.subtitle =  element_text(size = 24),
  plot.caption = element_text(size = 16), 
  legend.text = element_text(size = 20), 
  legend.key.size = unit(2, 'lines'), 
  legend.key=element_blank())

out_list <- list()
union_list <- list()
rpm <- function(X, reads){
  X[,reads] <- as.numeric(X[,reads])
  return(X[,reads]/colSums(as.matrix(X[,reads]))*100)
}

bed <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/')
for(j in seq(2,5)) {
  library(dplyr)
  
  ## plotting
  frag_10bp <- bed %>%
    specific_frag_format_organ(frag = j) %>% 
    specific_reads_sgRNA_identify(gap = 50, length = 31335
                                  , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                  , leader = 72, UTR5 = 210, UTR3 = 31035, frag = j) %>%
    cut_by_n_bp(genome_length = 31335, bin = 10, frag=j)
  
  
  out_list <- list()
  inter_count <- ''
  x=0
  for(i in seq(1, length(frag_10bp), 2)) {
    lib1 <- frag_10bp[[i]]
    # lib1 <- lib1[!(lib1$class %in% 'small_gap'), ]
    lib1[rownames(lib1[lib1$sub_class %in% '',]),'sub_class'] <- ' '
    # A <- lib1[lib1$sub_class %in% 'uncanonical sgRNAs', ]
    # lib1 <- lib1[!(lib1$class %in% 'sgmRNA'), ]
    # lib1 <- rbind(lib1,A)
    lib2 <- frag_10bp[[i+1]]
    # lib2 <- lib2[!(lib2$class %in% 'small_gap'), ]
    lib2[rownames(lib2[lib2$sub_class %in% '', ]),'sub_class'] <- ' '
    # B <- lib2[lib2$sub_class %in% 'uncanonical sgRNAs', ]
    # lib2 <- lib2[!(lib2$class %in% 'sgmRNA'), ]
    # lib2 <- rbind(lib2,B)
    colNames <-  colnames(lib1)[c(2:(j*2+1),(j*2+6), (j*2+7))] # could be any number of column names here
    lib1$recom = do.call(paste, c(lib1[colNames], sep="-"))
    lib1 <- as.data.frame(table(lib1$recom), stringsAsFactors = F)
    lib2$recom = do.call(paste, c(lib2[colNames], sep="-"))
    lib2 <- as.data.frame(table(lib2$recom), stringsAsFactors = F)
    inter <- intersect(lib1$Var1,lib2$Var1)
    lib1_inter <- lib1[lib1$Var1 %in% inter, ]
    lib2_inter <- lib2[lib2$Var1 %in% inter, ]
    A <- cbind(lib1_inter, lib2_inter[,2])
    if(nrow(A) >0) {
      A$Group <- 'Intersection'
      colnames(A) <- c('recom', 'lib1', 'lib2', 'Group')
    }
    if(nrow(A)==0) {
      A <- as.data.frame(matrix('', nrow = 0, ncol = 4))
    }
    lib1_only <- lib1[lib1$Var1 %in% setdiff(lib1$Var1, inter), ]
    lib1_only$lib2 <- 0
    lib1_only$Group <- 'Rep1_only'
    colnames(lib1_only) <- c('recom', 'lib1', 'lib2', 'Group')
    lib2_only <- lib2[lib2$Var1 %in% setdiff(lib2$Var1, inter), ]
    lib2_only$lib2 <- lib2_only$Freq
    lib2_only$Freq <- 0
    lib2_only$Group <- 'Rep2_only'
    colnames(lib2_only) <- c('recom', 'lib1', 'lib2','Group')
    union <- rbind(A, lib1_only, lib2_only)
    if(nrow(union) > 0) {
      x=x+1
      inter_count[x] <- nrow(A)
      A <- do.call(rbind, strsplit(union[,1], split = '-'))
      A <- cbind(A, union[,2:4])
      out <- as.data.frame(do.call(paste, c(A[colnames(A)[c(1:(j*2))]], sep="-")))
      out <- cbind(out, A[,c((j*2+1):ncol(A))])
      colnames(out) <- c('share_recom', 'class','sub_class', 'lib1', 'lib2', 'Group')
      # out <- out[order(out$lib1, decreasing = T),]
      out_list[[x]] <- out
      names(out_list)[x] <- sub(names(frag_10bp)[i]
                                , pattern = '[0-9]', replacement = '')
    }
    # uni <- union(lib1,lib2)
    # union_list[[(i+1)/2]] <- uni
    # names(union_list)[(i+1)/2] <- sub(names(frag_10bp)[i]
    #                                   , pattern = '[0-9]', replacement = '')
  }
  
  
  # library(ggrepel)
  # plot_list <- list()
  # k=0
  # for(i in seq_along(out_list)) {
  #   two <- out_list[[i]]
  #   label <- NA
  #   two$index <- seq(1, nrow(two))
  #   if(nrow(two) >0 ){
  #     k=k+1
  #     if(inter_count[i]>1) {
  #       label <- round(cor.test(two[1:inter_count[i],5], two[1:inter_count[i],4], method = c("kendall"))$estimate, digits = 2)  
  #     }
  #     two$sum <- two[,4]+two[,5]
  #     two <- two[order(two$sum, decreasing = T), ]
  #     if(nrow(two)>20) {
  #       two[1:20,'plot_label'] <- paste(two[1:20,1], '\n', two[1:20,2], '\n', two[1:20,3], sep = '')
  #     }
  #     if(nrow(two)<20) {
  #       two[1:nrow(two),'plot_label'] <- paste(two[1:nrow(two),1], '\n'
  #                                              , two[1:nrow(two),2], '\n', two[1:nrow(two),3], sep = '')
  #     }
  #     color <- c('black','#93E426','#7726E4')
  #     if(length(table(two$Group))!=3) {
  #       color <- c('#93E426','#7726E4')
  #     }
  #     two <- two[order(two$index, decreasing = T), ]
  #     plot <- ggplot(two, aes(x=lib1, y=lib2, label=plot_label, col=Group))+
  #       geom_point()+theme_bw()+gg_theme+
  #       # geom_smooth(method = "lm", se = FALSE, col='Red')+
  #       geom_label_repel(size=4, show.legend = FALSE
  #                        , min.segment.length = 0, seed = 713
  #                        , nudge_x = (max(two$lib1)-0.4*max(two$lib1))
  #                        , nudge_y = (max(two$lib2)-0.4*max(two$lib2)),segment.curvature = -0.1,
  #                        segment.ncp = 3,
  #                        segment.angle = 20)+
  #       labs(
  #         x = 'Frequency in replication 1 (reads percentage)',
  #         y = 'Frequency in replication 2 (reads percentage)',
  #         title = paste(names(out_list)[i], ', frag ', j
  #                       , " reads, 10 bps cluster,"
  #                       , " recombination intersection between replication", sep = '')
  #         , subtitle = paste('Total 2 fragments recombination reads'
  #                            , ', kendall correlation between intersection recombination = '
  #                            , label, sep = '')
  #         , caption = 'small gap is included')+
  #       # scale_x_continuous(breaks = seq(from = 0, to = max(two$lib1), by = 5))+
  #       # scale_y_continuous(breaks = seq(from = 0, to = max(two$lib2), by = 5))+
  #       scale_color_manual(values = color)
  #     plot_list[[k]] <- plot
  #     names(plot_list)[k] <- names(out_list)[i]
  #   }
  # }
  # 
  # for(i in seq_along(plot_list)) {
  #   png(paste(names(plot_list)[i], ' frag_', j,
  #             ' MHV Total reads recombination.png', sep = '')
  #       , width = 1600, height = 900)
  #   print(plot_list[[i]])
  #   dev.off()
  # }
  
  for(i in seq_along(out_list)) {
    A <- out_list[[i]]
    A <- A[order(A$lib1, decreasing = T),]
    if(nrow(A) > 0) {
      write.csv(A
                , file = paste(names(out_list)[i], ' frag_', j
                               , '_MHV_10bp_cluster_recom_class_union.csv', sep = '')
                , quote = F, row.names = F)
    }
  }
}


###########################################################################################
# frag1 reproduction
lib_recom_matrix <- function(lib_list, frag, libs) {
  recom_table_list <- list()
  recom_list <- list()
  inter_count <- ''
  x=0
  for(i in seq_along(lib_list)) {
    lib1 <- lib_list[[i]]
    lib1 <- lib1[!(lib1$class %in% 'small_gap'), ]
    lib1[rownames(lib1[lib1$sub_class %in% '',]),'sub_class'] <- ' '
    if(frag !=1) {
      colNames <-  colnames(lib1)[c(2:(frag*2+1),(frag*2+6), (frag*2+7))]
    }
    if(frag ==1) {
      colNames <-  colnames(lib1)[c(2,3,8,9)]
    }
    lib1$recom <-  do.call(paste, c(lib1[colNames], sep="-"))
    lib1 <- as.data.frame(table(lib1$recom), stringsAsFactors = F)
    recom_list[[i]] <- lib1$Var1
    recom_table_list[[i]] <- lib1
  }
  names(recom_list) <- names(lib_list)
  names(recom_table_list) <- names(lib_list)
  
  all_recom <- Reduce(union, recom_list)
  lib_matrix <- as.data.frame(matrix(data = 0, nrow = length(all_recom), ncol = libs))
  lib_matrix$index <- seq(1, nrow(lib_matrix))
  rownames(lib_matrix) <- all_recom; colnames(lib_matrix) <- names(lib_list)
  
  for(i in seq_along(recom_table_list)) {
    B <- recom_table_list[[i]]
    B <- B[order(B$Var1, decreasing = T), ]
    index <- lib_matrix[rownames(lib_matrix) %in% B$Var1, (ncol(lib_matrix))]
    C <- lib_matrix[lib_matrix[,ncol(lib_matrix)] %in% index,]
    C <- C[order(rownames(C), decreasing = T), ]
    lib_matrix <- lib_matrix[!(lib_matrix[,ncol(lib_matrix)] %in% index), ]
    print(table(rownames(C)== B$Var1))
    C[,i] <- B$Freq
    lib_matrix <- rbind(lib_matrix,C)
  }
  lib_matrix <- lib_matrix[, -c(ncol(lib_matrix))]
}



lin_frag1_identification <- function(list, recom_matrix, gff_path, libs
                                     , UTR5, UTR3, length, sgmRNA_BA=50
                                     , whole_genome_tolerance=20) {
  unc <- list()
  for(j in seq(1,libs,2)) {
    B <- recom_matrix[,j:(j+1)]
    C <- as.data.frame(ifelse(test = B>=2, yes = 1, no = 0))
    C$con <- rowSums(as.matrix(C[,1:2]))
    C <- C[C$con > 0, ]
    C <- C[grep(rownames(C), pattern = 'uncanonical'), ]
    C <- C[C$con >1, ]
    C <- as.data.frame(do.call(rbind, strsplit(rownames(C), split = '-')))
    unc[[j]] <- rep(as.numeric(unique(C$V3)))
    unc[[j+1]] <- rep(as.numeric(unique(C$V3)))
  }
  
  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('plyr', 'dplyr')
  ipak(packages)
  
  gff <- read.csv(file = gff_path, header = F)
  
  output_list <- list()
  for(i in seq_along(list))  {
    f1 <- list[[i]]
    f1$length <- f1[,3]-f1[,2]
    f1$complete <- f1$length > (length-whole_genome_tolerance)
    complete <- f1[f1$complete %in% TRUE,-c(ncol(f1))]
    if(nrow(complete)>0){
      complete$class <- c('complete genome')
    }
    if(nrow(complete)==0){
      complete <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- f1[f1$complete %in% FALSE,-c(ncol(f1))]
    f1$utr3 <- f1$end> UTR3
    delete_utr3 <- f1[f1$utr3 %in% FALSE,-c(ncol(f1))]
    if(nrow(delete_utr3)>0){
      delete_utr3$utr5 <- delete_utr3$start < UTR5
      delete_utr3_utr5 <- delete_utr3[delete_utr3$utr5 %in% FALSE, -c(ncol(delete_utr3))]
      if(nrow(delete_utr3_utr5)>0){
        delete_utr3_utr5$class <- c('delete 3\'UTR and 5\' UTR genome')
      }
      if(nrow(delete_utr3_utr5)==0)
      {delete_utr3_utr5 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
      
      delete_utr3 <- delete_utr3[delete_utr3$utr5 %in% TRUE, -c(ncol(delete_utr3))]
      if(nrow(delete_utr3)>0){
        delete_utr3$class <- c('delete 3\'UTR genome')
      }
      if(nrow(delete_utr3)==0)
      {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    }
    if(nrow(delete_utr3)==0)
    {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- f1[f1$utr3 %in% TRUE,-c(ncol(f1))]
    f1$utr5 <- f1$start < UTR5
    delete_utr5 <- f1[f1$utr5 %in% FALSE,-c(ncol(f1))]
    if(nrow(delete_utr5)>0)
    {
      delete_utr5$class <- c('delete 5\' terminal')
    }
    if(nrow(delete_utr5)==0)
    {delete_utr5 <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    near_complete <- f1[f1$utr5 %in% T,-c(ncol(f1))]
    if(nrow(near_complete)>0)
    {
      near_complete$class <- c('near complete genome')
    }
    if(nrow(near_complete)==0)
    {near_complete <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- rbind(delete_utr3, delete_utr3_utr5, delete_utr5, near_complete)
    output_list[[i]] <- f1
    names(output_list)[[i]] <- names(list)[[i]]
  }
  
  for(k in 1:length(output_list)) {
    B <- output_list[[k]]
    rownames(B) <- seq(1,nrow(B))
    if(nrow(B[B$class %in% 'delete 5\' terminal', ])>0) {
      A <- B[!(B$class %in% 'delete 5\' terminal'),]
      sgmRNA <- list()
      A$sub_class <- c('')
      D <- B[B$class %in% c('delete 5\' terminal'), ]
      if(nrow(D)>0) {
        for(j in seq_len(nrow(gff))){
          D$sub_class_bulin <- (gff[j,4]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(gff[j,4]+sgmRNA_BA)
          E <- D[D$sub_class_bulin==T,]
          E <- E[,-c(ncol(E))]
          if(nrow(E)>0) {
            E$sub_class <- paste('delete leader_', gff[j,9], sep = '')
          }
          sgmRNA[[j]] <- E
          D <- D[D$sub_class_bulin ==F,]
        }
        D <- D[D$sub_class_bulin ==F,]
        D <- D[,-c(ncol(D))]
        if(nrow(D)>0) {
          unc_list <- list()
          unc_cor <- unc[[k]]
          if(length(unc_cor) != 0) {
            for(j in seq_len(length(unc_cor))) {
              D$sub_class_bulin <- (unc_cor[j]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(unc_cor[j]+sgmRNA_BA)
              E <- D[D$sub_class_bulin==T,-c(ncol(D))]
              if(nrow(E)>0) {
                E$sub_class <- paste('delete leader uncanonincal_', unc_cor[j], sep = '')
              }
              unc_list[[j]] <- E
              D <- D[D$sub_class_bulin ==F,]
            }
            D <- D[,-c(ncol(D))]
            if(nrow(D)>0) {
              D$sub_class <- c('delete 5\' UTR genome')
            }
            D <- rbind(do.call(rbind, sgmRNA), do.call(rbind, unc_list), D, A)
            output_list[[k]] <- D 
          }
          else if(length(unc_cor) == 0){
            if(nrow(D)>0) {
              D$sub_class <- c('delete 5\' UTR genome')
            }
            D <- rbind(do.call(rbind, sgmRNA),D, A)
            output_list[[k]] <- D 
          }
        }
      }
    }
    if(nrow(B[B$class %in% 'delete 5\' terminal', ])==0) {
      B$sub_class <- c('')
      B <- B[order(B$freq, decreasing = T), ]
      output_list[[k]] <- B
    }
  }
  return(output_list)
}



# function6 : construct recombination point sort and new seq name
one_frag_recom_freq <- function(list)  {
  recom <- list()
  for(i in seq_along(list))  {
    B <- list[[i]]
    B[B$sub_class %in% '','sub_class'] <- '.'
    colNames <-  colnames(B)[c(2:(1*2+1),(1*2+6), (1*2+7))] # could be any number of column names here
    B$recom <- do.call(paste, c(B[colNames], sep="-"))
    recom_table <- as.data.frame(table(B$recom), stringsAsFactors = F)
    recom_table <- recom_table[order(recom_table$Freq, decreasing = T),]
    frequency <- recom_table$Freq
    out <- as.data.frame(do.call(rbind, strsplit(recom_table$Var1, split = '-')))
    out$name <- recom_table$Var1
    out$freq <- recom_table$Freq
    
    recom_table <- as.data.frame(matrix('U00735.2', nrow = nrow(out), ncol = 1))
    recom_table[,2:(1+6)] <- out[,1:ncol(out)]
    colnames(recom_table)[seq(2, (1+(1*2)), 2)] <- paste('start', 1, sep = '')
    colnames(recom_table)[seq(3, (1+(1*2)), 2)] <- paste('end', 1, sep ='')
    recom_table$score <- 0
    recom_table$strand <- c('+')
    for(j in seq(2, (1*2+1),1))
    {recom_table[,j] <- as.numeric(recom_table[,j])}
    D <- 0
    recom_table$length <- recom_table$end1- recom_table$start1+1
    recom_table$frag <- 1
    recom_table <- recom_table[order(recom_table$freq, decreasing = T), ]
    recom[[i]] <- recom_table
    names(recom)[i] <- names(list)[i]
  }
  return(recom)
}



frag2 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 2) %>% 
  specific_reads_sgRNA_identify(gap = 50, length = 31335
                                , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                , leader = 72, UTR5 = 210, UTR3 = 31035, frag = 2)

lib_matrix <- lib_recom_matrix(lib_list = frag2, frag = 2, libs = 4)

frag1 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 1)  %>%
  lin_frag1_identification(recom_matrix = lib_matrix, libs = 4
                           , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                           , UTR5 = 210, UTR3 = 31035, length = 31335) %>%
  cut_by_n_bp(genome_length = 31335, bin = 10, frag = 1)

for(j in 1) {
  out_list <- list()
  inter_count <- ''
  x=0
  for(i in seq(1, length(frag1), 2)) {
    lib1 <- frag1[[i]]
    # lib1 <- lib1[!(lib1$class %in% 'small_gap'), ]
    lib1[rownames(lib1[lib1$sub_class %in% '',]),'sub_class'] <- ' '
    # A <- lib1[lib1$sub_class %in% 'uncanonical sgRNAs', ]
    # lib1 <- lib1[!(lib1$class %in% 'sgmRNA'), ]
    # lib1 <- rbind(lib1,A)
    lib2 <- frag1[[i+1]]
    # lib2 <- lib2[!(lib2$class %in% 'small_gap'), ]
    lib2[rownames(lib2[lib2$sub_class %in% '', ]),'sub_class'] <- ' '
    # B <- lib2[lib2$sub_class %in% 'uncanonical sgRNAs', ]
    # lib2 <- lib2[!(lib2$class %in% 'sgmRNA'), ]
    # lib2 <- rbind(lib2,B)
    colNames <-  colnames(lib1)[c(2,3,8,9)] # could be any number of column names here
    lib1$recom = do.call(paste, c(lib1[colNames], sep="-"))
    lib1 <- as.data.frame(table(lib1$recom), stringsAsFactors = F)
    lib2$recom = do.call(paste, c(lib2[colNames], sep="-"))
    lib2 <- as.data.frame(table(lib2$recom), stringsAsFactors = F)
    inter <- intersect(lib1$Var1,lib2$Var1)
    lib1_inter <- lib1[lib1$Var1 %in% inter, ]
    lib2_inter <- lib2[lib2$Var1 %in% inter, ]
    A <- cbind(lib1_inter, lib2_inter[,2])
    if(nrow(A) >0) {
      A$Group <- 'Intersection'
      colnames(A) <- c('recom', 'lib1', 'lib2', 'Group')
    }
    if(nrow(A)==0) {
      A <- as.data.frame(matrix('', nrow = 0, ncol = 4))
    }
    lib1_only <- lib1[lib1$Var1 %in% setdiff(lib1$Var1, inter), ]
    lib1_only$lib2 <- 0
    lib1_only$Group <- 'Rep1_only'
    colnames(lib1_only) <- c('recom', 'lib1', 'lib2', 'Group')
    lib2_only <- lib2[lib2$Var1 %in% setdiff(lib2$Var1, inter), ]
    lib2_only$lib2 <- lib2_only$Freq
    lib2_only$Freq <- 0
    lib2_only$Group <- 'Rep2_only'
    colnames(lib2_only) <- c('recom', 'lib1', 'lib2','Group')
    union <- rbind(A, lib1_only, lib2_only)
    if(nrow(union) > 0) {
      x=x+1
      inter_count[x] <- nrow(A)
      A <- do.call(rbind, strsplit(union[,1], split = '-'))
      A <- cbind(A, union[,2:4])
      out <- as.data.frame(do.call(paste, c(A[colnames(A)[c(1:(j*2))]], sep="-")))
      out <- cbind(out, A[,c((j*2+1):ncol(A))])
      colnames(out) <- c('share_recom', 'class','sub_class', 'lib1', 'lib2', 'Group')
      # out <- out[order(out$lib1, decreasing = T),]
      out_list[[x]] <- out
      names(out_list)[x] <- sub(names(frag1)[i]
                                , pattern = '[0-9]', replacement = '')
    }
    # uni <- union(lib1,lib2)
    # union_list[[(i+1)/2]] <- uni
    # names(union_list)[(i+1)/2] <- sub(names(frag1)[i]
    #                                   , pattern = '[0-9]', replacement = '')
  }
  
  
  # library(ggrepel)
  # plot_list <- list()
  # k=0
  # for(i in seq_along(out_list)) {
  #   two <- out_list[[i]]
  #   label <- NA
  #   two$index <- seq(1, nrow(two))
  #   if(nrow(two) >0 ){
  #     k=k+1
  #     if(inter_count[i]>1) {
  #       label <- round(cor.test(two[1:inter_count[i],5], two[1:inter_count[i],4], method = c("kendall"))$estimate, digits = 2)  
  #     }
  #     two$sum <- two[,4]+two[,5]
  #     two <- two[order(two$sum, decreasing = T), ]
  #     if(nrow(two)>20) {
  #       two[1:20,'plot_label'] <- paste(two[1:20,1], '\n', two[1:20,2], '\n', two[1:20,3], sep = '')
  #     }
  #     if(nrow(two)<20) {
  #       two[1:nrow(two),'plot_label'] <- paste(two[1:nrow(two),1], '\n'
  #                                              , two[1:nrow(two),2], '\n', two[1:nrow(two),3], sep = '')
  #     }
  #     color <- c('black','#93E426','#7726E4')
  #     if(length(table(two$Group))!=3) {
  #       color <- c('#93E426','#7726E4')
  #     }
  #     two <- two[order(two$index, decreasing = T), ]
  #     plot <- ggplot(two, aes(x=lib1, y=lib2, label=plot_label, col=Group))+
  #       geom_point()+theme_bw()+gg_theme+
  #       # geom_smooth(method = "lm", se = FALSE, col='Red')+
  #       geom_label_repel(size=4, show.legend = FALSE
  #                        , min.segment.length = 0, seed = 713
  #                        , nudge_x = (max(two$lib1)-0.4*max(two$lib1))
  #                        , nudge_y = (max(two$lib2)-0.4*max(two$lib2)),segment.curvature = -0.1,
  #                        segment.ncp = 3,
  #                        segment.angle = 20)+
  #       labs(
  #         x = 'Frequency in replication 1 (reads percentage)',
  #         y = 'Frequency in replication 2 (reads percentage)',
  #         title = paste(names(out_list)[i], ', frag ', j
  #                       , " reads, 10 bps cluster,"
  #                       , " recombination intersection between replication", sep = '')
  #         , subtitle = paste('Total 2 fragments recombination reads'
  #                            , ', kendall correlation between intersection recombination = '
  #                            , label, sep = '')
  #         , caption = 'small gap is included')+
  #       # scale_x_continuous(breaks = seq(from = 0, to = max(two$lib1), by = 5))+
  #       # scale_y_continuous(breaks = seq(from = 0, to = max(two$lib2), by = 5))+
  #       scale_color_manual(values = color)
  #     plot_list[[k]] <- plot
  #     names(plot_list)[k] <- names(out_list)[i]
  #   }
  # }
  # 
  # 
  # 
  # for(i in seq_along(plot_list)) {
  #   png(paste(names(plot_list)[i], ' frag_', j,
  #             ' BCoV Total reads recombination.png', sep = '')
  #       , width = 1600, height = 900)
  #   print(plot_list[[i]])
  #   dev.off()
  # }
  
  for(i in seq_along(out_list)) {
    A <- out_list[[i]]
    A <- A[order(A$lib1, decreasing = T),]
    if(nrow(A) > 0) {
      write.csv(A
                , file = paste(names(out_list)[i], ' frag_', j
                               , '_MHV_10bp_cluster_recom_class_union.csv', sep = '')
                , quote = F, row.names = F)
    }
  }
}


############################################################################################
## qname extraction
A <-  as.data.frame(matrix(data = '', nrow = 0, ncol = 0))
out_list <- list(A,A,A,A)
bed <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/')
for(j in 2) {
  frag <- bed%>%
    specific_frag_format_organ(frag = j) %>% 
    specific_reads_sgRNA_identify(gap = 50, length = 31335
                                  , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                  , leader = 72, UTR5 = 210, UTR3 = 31035, frag = j)
  names(out_list) <- names(frag)
  for(i in seq_along(frag)) {
    frag_extract <- frag[[i]]
    frag_extract <- frag_extract[frag_extract$class %in% 'sgmRNA', ]
    f2_di <- as.data.frame(frag_extract$V4[!(frag_extract$sub_class %in% 'uncanonical sgRNAs')])
    print(nrow(f2_di))
    A <- out_list[[i]]
    A <- rbind(A, f2_di)
    out_list[[i]] <- A
  }
}


for(i in seq_along(out_list)) {
  A <- out_list[[i]]
  if(nrow(A) > 0) {
    write.table(A
                , file = paste(names(out_list)[i]
                               , '_MHV_sgmRNA_ca_qname', sep = '')
                , quote = F, row.names = F, col.names = F)
  }
}





############################################################################################
## qname extraction
frag2 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 2) %>% 
  specific_reads_sgRNA_identify(gap = 50, length = 31335
                                , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                , leader = 72, UTR5 = 210, UTR3 = 31035, frag = 2)

lib_matrix <- lib_recom_matrix(lib_list = frag2, frag = 2, libs=4)

bed <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/')
frag <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 1)  %>%
  lin_frag1_identification(recom_matrix = lib_matrix, libs = 4
                           , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                           , UTR5 = 210, UTR3 = 31035, length = 31335)

A <-  as.data.frame(matrix(data = '', nrow = 0, ncol = 0))
out_list <- list(A,A,A,A)
names(out_list) <- names(frag)
for(i in seq_along(frag)) {
  frag_extract <- frag[[i]]
  frag_extract <- frag_extract[frag_extract$class %in% 'delete 5\' terminal', ]
  f2_di <- as.data.frame(frag_extract$V4[grep(frag_extract$sub_class, pattern = '^delete leader uncanonincal')])
  # f2_di <- as.data.frame(frag_extract$V4[-c(grep(frag_extract$sub_class, pattern = 'delete leader uncanonincal'))])
  print(nrow(f2_di))
  A <- out_list[[i]]
  A <- rbind(A, f2_di)
  out_list[[i]] <- A
}


for(i in seq_along(out_list)) {
  A <- out_list[[i]]
  if(nrow(A) > 0) {
    write.table(A
                , file = paste(names(out_list)[i]
                               , '_MHV_1f_unc_qname', sep = '')
                , quote = F, row.names = F, col.names = F)
  }
}



############################################################################################
lib_recom_matrix <- function(lib_list, frag, libs) {
  recom_table_list <- list()
  recom_list <- list()
  inter_count <- ''
  x=0
  for(i in seq_along(lib_list)) {
    lib1 <- lib_list[[i]]
    lib1 <- lib1[!(lib1$class %in% 'small_gap'), ]
    lib1[rownames(lib1[lib1$sub_class %in% '',]),'sub_class'] <- ' '
    if(frag !=1) {
      colNames <-  colnames(lib1)[c(2:(frag*2+1),(frag*2+6), (frag*2+7))]
    }
    if(frag ==1) {
      colNames <-  colnames(lib1)[c(2,3,8,9)]
    }
    lib1$recom <-  do.call(paste, c(lib1[colNames], sep="-"))
    lib1 <- as.data.frame(table(lib1$recom), stringsAsFactors = F)
    recom_list[[i]] <- lib1$Var1
    recom_table_list[[i]] <- lib1
  }
  names(recom_list) <- names(lib_list)
  names(recom_table_list) <- names(lib_list)
  
  all_recom <- Reduce(union, recom_list)
  lib_matrix <- as.data.frame(matrix(data = 0, nrow = length(all_recom), ncol = libs))
  lib_matrix$index <- seq(1, nrow(lib_matrix))
  rownames(lib_matrix) <- all_recom; colnames(lib_matrix) <- names(lib_list)
  
  for(i in seq_along(recom_table_list)) {
    B <- recom_table_list[[i]]
    B <- B[order(B$Var1, decreasing = T), ]
    index <- lib_matrix[rownames(lib_matrix) %in% B$Var1, (ncol(lib_matrix))]
    C <- lib_matrix[lib_matrix[,ncol(lib_matrix)] %in% index,]
    C <- C[order(rownames(C), decreasing = T), ]
    lib_matrix <- lib_matrix[!(lib_matrix[,ncol(lib_matrix)] %in% index), ]
    print(table(rownames(C)== B$Var1))
    C[,i] <- B$Freq
    lib_matrix <- rbind(lib_matrix,C)
  }
  lib_matrix <- lib_matrix[, -c(ncol(lib_matrix))]
}



for(j in seq(2,5)) {
  
  lib_matrix <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
    specific_frag_format_organ(frag = j) %>% 
    specific_reads_sgRNA_identify(gap = 50, length = 31335
                                  , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                  , leader = 72, UTR5 = 210, UTR3 = 31035, frag = j) %>%
    cut_by_n_bp(genome_length = 31335, bin = 10, frag=j) %>%
    lib_recom_matrix(frag = j, libs = 4)
  
  
  B <- lib_matrix
  B <- as.data.frame(cbind(rowSums(B[,1:2]), rowSums(B[,3:4])))
  colnames(B) <- c('cell', 'liver')
  
  C <- as.data.frame(do.call(rbind, strsplit(rownames(B), split = '-')))
  C[,colnames(B)] <- B[,1:2]
  C$v5 <- 'intersect'
  C[(C$cell !=0 &C$liver==0),'v5'] <- 'cell'
  C[(C$liver !=0 &C$cell==0),'v5'] <- 'liver'
  C <- C[order(C$cell, decreasing = T), ]
  
  if(nrow(C) > 0) {
    write.csv(C, file = paste('~/Analysis/MHV/result/f1_0515/'
                              , 'cell-liver', 'frag_', j, '.csv'))
  }
}


##### frag1

frag2 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 2) %>% 
  specific_reads_sgRNA_identify(gap = 50, length = 31335
                                , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                                , leader = 72, UTR5 = 210, UTR3 = 31035, frag = 2)

lib_matrix <- lib_recom_matrix(lib_list = frag2, frag = 2, libs = 4)

j=1
lib_matrix_1 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 1)  %>%
  lin_frag1_identification(recom_matrix = lib_matrix, libs = 4
                           , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                           , UTR5 = 210, UTR3 = 31035, length = 31335) %>%
  cut_by_n_bp(genome_length = 31335, bin = 10, frag=j) %>%
  lib_recom_matrix(frag = j, libs = 4)


B <- lib_matrix_1
B <- as.data.frame(cbind(rowSums(B[,1:2]), rowSums(B[,3:4])))
colnames(B) <- c('cell', 'liver')

C <- as.data.frame(do.call(rbind, strsplit(rownames(B), split = '-')))
C[,colnames(B)] <- B[,1:2]
C$v5 <- 'intersect'
C[(C$cell !=0 &C$liver==0),'v5'] <- 'cell'
C[(C$liver !=0 &C$cell==0),'v5'] <- 'liver'
C <- C[order(C$cell, decreasing = T), ]
if(nrow(C) > 0) {
  write.csv(C, file = paste('~'
                            , 'cell-liver', 'frag_', j, '.csv'))
}




#########################################################################################

lin_frag1_identification <- function(list, recom_matrix, gff_path, libs
                                     , UTR5, UTR3, length, sgmRNA_BA=50
                                     , whole_genome_tolerance=20) {
  
  
  B <- recom_matrix[,1:2]
  C <- as.data.frame(ifelse(test = B>=2, yes = 1, no = 0))
  C$con <- rowSums(as.matrix(C[,1:2]))
  C <- C[C$con > 0, ]
  C <- C[grep(rownames(C), pattern = 'uncanonical'), ]
  C <- C[C$con >1, ]
  C <- as.data.frame(do.call(rbind, strsplit(rownames(C), split = '-')))
  unc <- as.numeric(unique(C$V3))

  
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  packages <- c('plyr', 'dplyr')
  ipak(packages)
  
  gff <- read.csv(file = gff_path, header = F)
  
  output_list <- list()
  for(i in seq_along(list))  {
    f1 <- list[[i]]
    f1$length <- f1[,3]-f1[,2]
    f1$complete <- f1$length > (length-whole_genome_tolerance)
    complete <- f1[f1$complete %in% TRUE,-c(ncol(f1))]
    if(nrow(complete)>0){
      complete$class <- c('complete genome')
    }
    if(nrow(complete)==0){
      complete <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- f1[f1$complete %in% FALSE,-c(ncol(f1))]
    f1$utr3 <- f1$end> UTR3
    delete_utr3 <- f1[f1$utr3 %in% FALSE,-c(ncol(f1))]
    if(nrow(delete_utr3)>0){
      delete_utr3$utr5 <- delete_utr3$start < UTR5
      delete_utr3_utr5 <- delete_utr3[delete_utr3$utr5 %in% FALSE, -c(ncol(delete_utr3))]
      if(nrow(delete_utr3_utr5)>0){
        delete_utr3_utr5$class <- c('delete 3\'UTR and 5\' UTR genome')
      }
      if(nrow(delete_utr3_utr5)==0)
      {delete_utr3_utr5 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
      
      delete_utr3 <- delete_utr3[delete_utr3$utr5 %in% TRUE, -c(ncol(delete_utr3))]
      if(nrow(delete_utr3)>0){
        delete_utr3$class <- c('delete 3\'UTR genome')
      }
      if(nrow(delete_utr3)==0)
      {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    }
    if(nrow(delete_utr3)==0)
    {delete_utr3 <- as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- f1[f1$utr3 %in% TRUE,-c(ncol(f1))]
    f1$utr5 <- f1$start < UTR5
    delete_utr5 <- f1[f1$utr5 %in% FALSE,-c(ncol(f1))]
    if(nrow(delete_utr5)>0)
    {
      delete_utr5$class <- c('delete 5\' terminal')
    }
    if(nrow(delete_utr5)==0)
    {delete_utr5 <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    near_complete <- f1[f1$utr5 %in% T,-c(ncol(f1))]
    if(nrow(near_complete)>0)
    {
      near_complete$class <- c('near complete genome')
    }
    if(nrow(near_complete)==0)
    {near_complete <-as.data.frame(matrix(data = '', nrow = 0, ncol = 8))}
    
    f1 <- rbind(delete_utr3, delete_utr3_utr5, delete_utr5, near_complete)
    output_list[[i]] <- f1
    names(output_list)[[i]] <- names(list)[[i]]
  }
  
  for(k in 1:length(output_list)) {
    B <- output_list[[k]]
    rownames(B) <- seq(1,nrow(B))
    if(nrow(B[B$class %in% 'delete 5\' terminal', ])>0) {
      A <- B[!(B$class %in% 'delete 5\' terminal'),]
      sgmRNA <- list()
      A$sub_class <- c('')
      D <- B[B$class %in% c('delete 5\' terminal'), ]
      if(nrow(D)>0) {
        for(j in seq_len(nrow(gff))){
          D$sub_class_bulin <- (gff[j,4]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(gff[j,4]+sgmRNA_BA)
          E <- D[D$sub_class_bulin==T,]
          E <- E[,-c(ncol(E))]
          if(nrow(E)>0) {
            E$sub_class <- paste('delete leader_', gff[j,9], sep = '')
          }
          sgmRNA[[j]] <- E
          D <- D[D$sub_class_bulin ==F,]
        }
        D <- D[D$sub_class_bulin ==F,]
        D <- D[,-c(ncol(D))]
        if(nrow(D)>0) {
          unc_list <- list()
          unc_cor <- unc
          if(length(unc_cor) != 0) {
            for(j in seq_len(length(unc_cor))) {
              D$sub_class_bulin <- (unc_cor[j]-sgmRNA_BA)<D[,'start1'] & D[,'start1']<(unc_cor[j]+sgmRNA_BA)
              E <- D[D$sub_class_bulin==T,-c(ncol(D))]
              if(nrow(E)>0) {
                E$sub_class <- paste('delete leader uncanonincal_', unc_cor[j], sep = '')
              }
              unc_list[[j]] <- E
              D <- D[D$sub_class_bulin ==F,]
            }
            D <- D[,-c(ncol(D))]
            if(nrow(D)>0) {
              D$sub_class <- c('delete 5\' UTR genome')
            }
            D <- rbind(do.call(rbind, sgmRNA), do.call(rbind, unc_list), D, A)
            output_list[[k]] <- D 
          }
          else if(length(unc_cor) == 0){
            if(nrow(D)>0) {
              D$sub_class <- c('delete 5\' UTR genome')
            }
            D <- rbind(do.call(rbind, sgmRNA),D, A)
            output_list[[k]] <- D 
          }
        }
      }
    }
    if(nrow(B[B$class %in% 'delete 5\' terminal', ])==0) {
      B$sub_class <- c('')
      B <- B[order(B$freq, decreasing = T), ]
      output_list[[k]] <- B
    }
  }
  return(output_list)
}


frag1 <- read_as_list(path = '~/Analysis/MHV/bed/cell_liver/v2/') %>%
  specific_frag_format_organ(frag = 1)  %>%
  lin_frag1_identification(recom_matrix = lib_matrix, libs = 4
                           , gff_path = '~/Analysis/reference/MHV_edit_gff.csv'
                           , UTR5 = 210, UTR3 = 31035, length = 31335) %>%
  cut_by_n_bp(genome_length = 31335, bin = 10, frag = 1)

