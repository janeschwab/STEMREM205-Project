#Load required packages

# List of packages
packages <- c("Seurat", "dplyr", "magrittr", "data.table", "Matrix", 
              "devtools", "RcppArmadillo", "Rcpp", "scales", "pheatmap", 
              "gplots", "ggplot2", "cowplot", "tibble", "xlsx")

# Install and load packages
#for (package in packages) {
#  if (!require(package, character.only = TRUE)) {
#    install.packages(package)
#    library(package, character.only = TRUE)
#  } else {
#    library(package, character.only = TRUE)
#  }
#}

# Install and load remotes package
library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(tibble)
library(data.table)

#Install batch correction package harmony
#install_github("immunogenomics/harmony")
#devtools::install_github("immunogenomics/harmony")
#remotes::install_github("immunogenomics/presto")
library(harmony)

#____________________________________________________________________________________________________________________________________________________________________#

# 1. Functions

#Interactome Code

check_genes <- function(genes, database, object_genes) {
  'Check to make sure that wanted genes are in reference you provide and in object
  
   Args:
   genes (chr vector): list of potential wanted genes
   database (data.frame): table with reference ligand/receptor pairs
   object_genes (chr vector): list of genes present in Seurat object
   
   Returns:
   not_found_list (chr list): list of genes not found in database and/or object
  '
  
  database_genes <- as.vector(unlist(database))
  not_found_database <- c()
  not_found_object <- c()
  
  for (i in 1:length(genes)) {
    if (!(genes[i] %in% database_genes)) {
      not_found_database <- c(not_found_database, genes[i])
    }
    
    if (!(genes[i] %in% object_genes)) {
      not_found_object <- c(not_found_object, genes[i])
    }
  }
  
  not_found_list <- list(database=not_found_database, object=not_found_object)
  return(not_found_list)
}

make_LR_pairs <- function(ligands, receptors, database) {
  'Make all LR pairs based on wanted ligands and/or receptors and database provided
  
   Args:
   ligands (chr vector): list of wanted ligands
   receptors (chr vector): list of wanted receptors
   database (data.frame): table with reference ligand/receptor pairs
   
   Returns:
   wanted_LR (data.frame): data.frame with ligand/receptor pairs from wanted ligands and receptors
  '
  
  wanted_LR <- data.frame(Ligand = character(0), Receptor = character(0))
  
  for (ligand in ligands){
    # list of corresponding receptors
    corresponding_receptors <- unique(as.character(database[,2][grep(ligand, database[,1])]))
    
    for (receptor in corresponding_receptors) {
      LR_row <- data.frame(Ligand = ligand, Receptor = receptor)
      wanted_LR <- rbind(wanted_LR, LR_row)
    }
  }
  
  # filter out unwanted receptors
  wanted_LR <- wanted_LR[which(wanted_LR$Receptor %in% wanted_receptors),]
  
  return(wanted_LR)
}

create_LR_table <- function(ligands, receptors, cell_types, LRs, avg0, avg1) {
  'Create table w/ ligand, receptor, source and target cells, average expressions, and IDs
   
   Args:
   ligands (chr vector): list of wanted ligands
   receptors (chr vector): list of wanted receptors
   cell_types (chr vector): list of common cell types between two objects
   LRs (data.frame): table with with wanted ligand/receptor pairs (from make_LR_pairs function)
   avg0 (list of num data.frame): average expression table of object 0
   avg1 (list of num data.frame): average expression table of object 1
   
   Returns:
   LR_table (data.frame): table of potential ligand/receptor pairs
    Contains ligand, receptor, source, target, average expression of ligand and receptors from 
    source and target cells. Also contains IDs (important for Cytoscape). 
  '
  
  LR_table <- data.frame()
  count <- 0
  
  for (ligand in ligands) {
    known_receptors <- LRs$Receptor[which(LRs$Ligand == ligand)]
    
    for (receptor in known_receptors) {
      for (i in c(1:length(cell_types))) {
        for (j in c(1:length(cell_types))) {
          LR_table <- rbind(LR_table,data.frame(ligands = ligand, receptors = receptor, 
                                                source = cell_types[i], target = cell_types[j], 
                                                avg_lig_0 = avg0[ligand, cell_types[i]],
                                                avg_rec_0 = avg0[receptor, cell_types[j]],
                                                avg_lig_1 = avg1[ligand, cell_types[i]],
                                                avg_rec_1 = avg1[receptor, cell_types[j]]))
        }
      }
    }
    
    cat(count, ' ')
    count <- count+1
  }
  
  # Create IDs
  source_ID <- c()
  target_ID <- c()
  
  # Find optimal cell type IDs
  ids_key <- c()
  
  for (i in 1:length(cell_types)) {
    for (j in 1:min(nchar(cell_types))) {
      name_1 <- substring(cell_types[i], 1, j)
      name_list <- c()
      name_list <- c(sapply(cell_types[-i], function(x) name_list <- c(name_list, substring(x, 1, j))))
      
      if(name_1 %in% name_list) {
        next
      }
      else {
        ids_key <- c(ids_key, name_1) 
        break
      }
    }
  }
  
  names(ids_key) <- cell_types
  
  for (i in c(1:length(LR_table[,1]))) {
    letter1 <- as.character(ids_key[LR_table$source[i]])
    letter2 <- as.character(ids_key[LR_table$target[i]])
    n1 <- which(ligands == LR_table$ligands[i])
    n2 <- which(receptors == LR_table$receptors[i])
    
    source_ID <- c(source_ID, paste(letter1, 'L', n1, sep = ''))
    target_ID <- c(target_ID, paste(letter2, 'R', n2, sep = ''))
  }
  
  LR_table <- cbind(LR_table, data.frame(source_ID = source_ID, target_ID = target_ID))
  
  return(LR_table)
}

avg_LR_filt <- function(table, threshold) {
  'Calculates states (ON/OFF --> 1/0) and filter if there is no expression in both groups
  
   Args:
   table (data.frame): table with potential ligand/receptors (from create_LR_table)
   threshold (num): average expression threshold 
   
   Returns:
   table (data.frame): filtered table based on average expression
  '
  
  # Find states of pairs in each group
  LR_states <- data.frame(lig_0 = character(0), rec_0 = character(0),
                          lig_1 = character(0), rec_1 = character(0),
                          is_on = logical(0))
  
  for (i in c(1:length(LR_table$avg_lig_0))) {
    row_states <- data.frame(lig_0 = table$avg_lig_0[i] > threshold, 
                             rec_0 = table$avg_rec_0[i] > threshold,
                             lig_1 = table$avg_lig_1[i] > threshold, 
                             rec_1 = table$avg_rec_1[i] > threshold)
    row_states$is_on <- (row_states$lig_0 & row_states$rec_0) | (row_states$lig_1 & row_states$rec_1)
    
    LR_states <- rbind(LR_states, row_states)
  }
  
  table <- cbind(table, ifelse((LR_states$lig_0 & LR_states$rec_0), 1, 0))
  table <- cbind(table, ifelse((LR_states$lig_1 & LR_states$rec_1), 1, 0))
  
  colnames(table)[11] <- 'is_on_0'
  colnames(table)[12] <- 'is_on_1'
  
  # Filter out pairs if pairs in both group are less than threshold 
  table <- table[which(LR_states$is_on == TRUE),]
  
  return(table)
}

LR_diff <- function(table, data_0, data_1, genes, label, alpha = 0.05) {
  'Calculate Wilcoxon-rank test on ligands and receptors (separately) between groups
    Order of test: data_0 --> data_1
  
    Args:
    table (data.frame): table with potential ligand/receptors (from create_LR_table)
    data_0 (data.frame): table of gene expression from object 0
    data_1 (data.frame): table of gene expression from object 1
    genes (chr vector): list of wanted genes
    label (chr): name of meta.data slot in Seurat object
    alpha (num): alpha level for Wilcoxon-rank test
    
    Returns:
    table (data.frame): filtered table based on Wilcoxon-rank test
  '
  
  table$lig_diff <- rep(0, length(table[,1]))
  table$rec_diff <- rep(0, length(table[,1]))
  
  table$lig_diff_p <- rep(0, length(table[,1]))
  table$rec_diff_p <- rep(0, length(table[,1]))
  
  for (i in 1:length(table[,1])) {
    ligand <- table$ligands[i]
    receptor <- table$receptors[i]
    source <- table$source[i]
    target <- table$target[i]
    
    lig_0_data <- data_0[which(data_0[,label] == source), ligand]
    rec_0_data <- data_0[which(data_0[,label] == target), receptor]
    
    lig_1_data <- data_1[which(data_1[,label] == source), ligand]
    rec_1_data <- data_1[which(data_1[,label] == target), receptor]
    
    lig_wilcox <- wilcox.test(lig_0_data, lig_1_data, exact = F, paired = F)
    rec_wilcox <- wilcox.test(rec_0_data, rec_1_data, exact = F, paired = F)
    
    table$lig_diff_p[i] <- lig_wilcox$p.value
    table$rec_diff_p[i] <- rec_wilcox$p.value
  }
  
  table$lig_diff_p_adj <- p.adjust(table$lig_diff_p, method = 'bonferroni')
  table$rec_diff_p_adj <- p.adjust(table$rec_diff_p, method = 'bonferroni')
  
  # If not significant, then 0
  # If significant, then 1
  for (i in 1:length(table[,1])) {
    table$lig_diff[i] <- ifelse(table$lig_diff_p_adj[i] < alpha, 1, 0)
    table$rec_diff[i] <- ifelse(table$rec_diff_p_adj[i] < alpha, 1, 0)
  }
  
  # # If there is difference, then find if increase/decrease
  # # for ligands
  # for (i in 1:length(table[,1])) {
  #   if (table$lig_diff[i] == 1) {
  #     table$lig_diff[i] <- ifelse(table$avg_lig_0[i] > table$avg_lig_1[i], yes = -1, no = 1)
  #   }
  #   
  #   if (table$rec_diff[i] == 1) {
  #     table$rec_diff[i] <- ifelse(table$avg_rec_0[i] > table$avg_rec_1[i], yes = -1, no = 1)
  #   }
  # }
  
  return(table)
}

generate_supplement <- function(LR_table, cytoscape_nodes) {
  'Generate supplemental information for nodes to input back into Cytoscape
    1) ligand/receptor 2) ligand/receptor name 
  
   Args:
    LR_table (data.frame): table with ligand/receptor pairs
    cytoscape_ids (data.frame): exported cytoscape node table
    
   Returns:
    node_names_ids (data.frame): table with gene names matching IDs and if ligand/receptor
  '
  
  # List of node names
  ids <- c(as.character(LR_table$source_ID), as.character(LR_table$target_ID))
  gene_names <- c(as.character(LR_table$ligands), as.character(LR_table$receptors))
  node_names_ids <- data.frame(ID = ids, names = gene_names, LR = rep(c('L', 'R'),each=length(LR_table$source_ID)),stringsAsFactors = F)
  node_names_ids <- unique(node_names_ids)
  
  # sort supplemental names in order of what is in cytoscape
  cytoscape_ids <- gsub('""', '', cytoscape_nodes$`shared name`)
  order <- match(cytoscape_ids, node_names_ids$ID)
  node_names_ids <- node_names_ids[order,]
  node_names_ids$ID <- paste('"', node_names_ids$ID, '"', sep = '')
  
  return(node_names_ids)
}

# AutomatedClusterMarkerTable returns FindAllMarkers table with extra bits of useful information
# and an educated guess about cluster identity

AutomatedClusterMarkerTable <- function(Seurat_Object){
  library(dplyr)
  library(tibble)
  library(Seurat)
  ClusterList <- list()
  Idents(object = Seurat_Object) <- "seurat_clusters"
  current.cluster.ids <- sort(as.numeric(levels(Seurat_Object@active.ident)))
  new.cluster.ids <- c()
  
  
  for(i in current.cluster.ids){
    List_Position <- i + 1
    ClusterList[[List_Position]] <- FindMarkers(object = Seurat_Object, ident.1 = i, min.pct = 0.25, only.pos = TRUE)
    Positive_Genes <- rownames(ClusterList[[List_Position]])
    Num_Positive_Genes <- length(Positive_Genes)
    
    RPS_Num <- length(grep(pattern = "^RPS", x = Positive_Genes))
    RPL_Num <- length(grep(pattern = "^RPL", x = Positive_Genes))
    RP_Percent <- sum(RPS_Num, RPL_Num)/length(Positive_Genes)*100
    RP_Label <- paste("RP%:", RP_Percent, sep = " ")
    
    Mito_Num <- length(grep(pattern = "^MT-", x = Positive_Genes))
    Mito_Percent <- Mito_Num/length(Positive_Genes)*100
    Mito_Label <- paste("Mito%:", RP_Percent, sep = " ")
    
    ClusterCells <- WhichCells(object = Seurat_Object, idents = i)
    Cell_Barcodes <- unlist(Seurat_Object@assays$RNA@counts@Dimnames[2])
    Cell_Number <- c()
    for(k in 1:length(ClusterCells)){
      Cell_Position <- grep(pattern = ClusterCells[k],x = Cell_Barcodes, value = FALSE)
      Cell_Number <- c(Cell_Number,Cell_Position)
    }
    
    
    S_Score <- Seurat_Object@meta.data$S.Score
    G2M_Score <- Seurat_Object@meta.data$G2M.Score
    Cluster_S_Score <- S_Score[Cell_Number]
    Cluster_G2M_Score <- G2M_Score[Cell_Number]
    Avg_Cluster_S_Score <- mean(Cluster_S_Score)
    Avg_Cluster_G2M_Score <- mean(Cluster_G2M_Score)
    Cluster_S_Score_Range <- range(Cluster_S_Score)
    Cluster_G2M_Score_Range <- range(Cluster_G2M_Score)
    
    nFeature <- Seurat_Object@meta.data$nFeature_RNA
    nCount <- Seurat_Object@meta.data$nCount_RNA
    Mito <- Seurat_Object@meta.data$percent.mt
    
    Cluster_nFeature <- nFeature[Cell_Number]
    Cluster_nCount <- nCount[Cell_Number]
    Cluster_Mito <- Mito[Cell_Number]
    
    Avg_Cluster_nFeature <- as.integer(mean(Cluster_nFeature))
    Avg_Cluster_nCount <- as.integer(mean(Cluster_nCount))
    Max_Cluster_Mito <- max(Cluster_Mito)
    
    Cell_Types <- c("Epi","T Cell","Myeloid","B Cell","Fibroblast","RBC","NK", "Endo","Acinar")
    
    Epi_Markers <- c("KRT7","KRT8","KRT18","KRT19","EPCAM","CDH1")
    T_Cell_Markers <- c("CD3E","CD3G","CD3D","CD4","IL7R","CD8A","LEF1")
    Myeloid_Markers <- c("CD14","ITGAM","MNDA","MPEG1","ITGAX")
    B_Cell_Markers <- c("CD79A","MS4A1","CD19")
    Fibroblast_Markers <- c("CDH11","PDGFRA","PDGFRB","ACTA2")
    RBC_Markers <- c("HBA1","HBB","HBA2")
    NK_Markers <- c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38","KLRC1")
    Endo_Markers <- c("CDH5","PECAM1")
    Acinar_Markers <- c("TRY4","SPINK1","AMY2A")
    All_Markers <- list(Epi_Markers,T_Cell_Markers,Myeloid_Markers,B_Cell_Markers,Fibroblast_Markers,RBC_Markers,NK_Markers,Endo_Markers,Acinar_Markers)
    
    Epi_Score <- 0
    T_Cell_Score <- 0
    Myeloid_Score <- 0
    B_Cell_Score <- 0
    Fibroblast_Score <- 0
    RBC_Score <- 0
    NK_Score <- 0
    Endo_Score <- 0
    Acinar_Score <- 0 
    All_Scores <- list(Epi_Score,T_Cell_Score,Myeloid_Score,B_Cell_Score,Fibroblast_Score,RBC_Score,NK_Score,Endo_Score,Acinar_Score)
    Weighted_Scores <- c()
    Score_Weights <- c(1.85,1.85,2.22,3.7,2.78,3.7,1.85,5.56,3.7) 
    
    for(h in 1:length(All_Markers)){
      Markers_to_Test<- All_Markers[[h]]
      Marker_Row <- h
      for(j in 1:length(Markers_to_Test)){
        Gene_Found <- 0
        Gene_Found <- length(grep(pattern = Markers_to_Test[j], x = Positive_Genes))
        if(Gene_Found > 0 ){
          All_Scores[[Marker_Row]] <- All_Scores[[Marker_Row]]+1
        }
      }
      Weighted_Scores[Marker_Row] <- All_Scores[[Marker_Row]]*Score_Weights[Marker_Row]
    }
    
    ClusterID <- which(Weighted_Scores >= 5.5)
    if(length(ClusterID) > 0){
      if(length(ClusterID) > 1){
        ID <- "Multiple"
      }else{
        ID <- Cell_Types[ClusterID]
      }
    }else{
      ID <- i
    } 
    if(RP_Percent > 30){
      ID <- paste("RP_",ID,sep = "")
    }
    if(Avg_Cluster_S_Score > 0.01 | Avg_Cluster_G2M_Score > 0.01){
      CellCycleID <- "Cycling"
      ID <- paste("Cycling_",ID,sep = "")
    }else{
      CellCycleID <- "N/A"
    }
    if(Avg_Cluster_nCount < 700){
      ID <- paste("G_",ID,sep = "")
    }
    new.cluster.ids <- c(new.cluster.ids,ID)
    
    Label_Row <- length(Positive_Genes) + 1
    Label_Row2 <- length(Positive_Genes) + 2
    Label_Row3 <- length(Positive_Genes) + 3
    Label_Row4 <- length(Positive_Genes) + 4
    Label_Row5 <- length(Positive_Genes) + 5
    
    Label1 <- c("Summary:",paste("Cluster",i, sep = " "), paste("ID:",ID, sep = " "),paste("Mito%:",Mito_Percent, sep = " "),paste("RP%:",RP_Percent, sep = " "))
    Label2 <- c("Immune Summary",paste("T Cell Score:",All_Scores[[2]], sep = " "),paste("Myeloid Score:",All_Scores[[3]], sep = " "),paste("B Cell Score:",All_Scores[[4]], sep = " "),
                paste("NK Score:",All_Scores[[7]], sep = " "))
    Label3 <- c(paste("Epi Score:",All_Scores[[1]], sep = " "),paste("Fib Score:",All_Scores[[5]], sep = " "),paste("Acinar Score:",All_Scores[[9]], sep = " "),
                paste("Endo Score:",All_Scores[[8]], sep = " "),paste("RBC Score:",All_Scores[[6]], sep = " "))
    Label4 <- c("Avg S Score:", Avg_Cluster_S_Score, "Avg G2M Score:", Avg_Cluster_G2M_Score, CellCycleID)
    Label5 <- c("Filter Info",paste("Avg. nGene:", Avg_Cluster_nFeature, sep = " "),paste("Avg. nCounts:", Avg_Cluster_nCount, sep = " ")
                , paste("Highest Mito:",Max_Cluster_Mito, sep = " "), paste("# Cells:",length(Cell_Number), sep = " "))
    
    ClusterList[[List_Position]][Label_Row,] <- Label1
    ClusterList[[List_Position]][Label_Row2,] <- Label2
    ClusterList[[List_Position]][Label_Row3,] <- Label3
    ClusterList[[List_Position]][Label_Row4,] <- Label4
    ClusterList[[List_Position]][Label_Row5,] <- Label5
    
    ClusterList[[List_Position]] <- rownames_to_column(.data = ClusterList[[List_Position]],var = "Gene")
    ClusterList[[List_Position]][Label_Row,"Gene"] <- "Summary1"
    ClusterList[[List_Position]][Label_Row2,"Gene"] <- "Summary2"
    ClusterList[[List_Position]][Label_Row3,"Gene"] <- "Summary3"
    ClusterList[[List_Position]][Label_Row4,"Gene"] <- "Summary4"
    ClusterList[[List_Position]][Label_Row5,"Gene"] <- "Summary5"
    
    
  }
  
  ClusterDataFrame <- bind_rows(ClusterList, .id = "column_label")
  ClusterDataFrame <- ClusterDataFrame[,-1]
  ClusterPackage <- list(ClusterDataFrame, new.cluster.ids)
  return(ClusterPackage)
}

#Circos Code

# CircosFunctions will output the text files required to run Circos
#SetIdents to Circos Labels Prior to Starting
CircosFunctions <- function(InteractomeData, SeuratObject, CellTypeVector, POI, Lig_or_Rec, LR_List, Species, Cutoff_Lig, Cutoff_Rec){
  #InteractomeData: Filtered Interactome Data
  #SeuratObject: Object used to make interactome 
  #CellTypeVector: Vector of cell types to be included in interactome
  #POI: A number representing the cell population of interest to single out as the Ligand in a receptor plot or the receptor in a ligand plot
  #Lig_or_Rec: T = POI as Ligand in receptor plot, F = POI as receptor in ligand plot
  #LR_List: A dataframe of ligands and paired receptors, ex. The Ramilowski List (organized human, mouse, suborganized ligand-receptor)
  #Species: T = Human, F = Mouse
  #Cutoff_Lig: Cutoff ligand expression value for the interactome
  #Cutoff_Rec: Cutoff receptor expression value for the interactome
  #Cutoff Values:
  #Positive = Cutoff determined by summary function ex Min, 1st Q, Median, Mean, 3rd Q, Max (1-6)
  #Negative = Arbitrary cutoff (set -0.05 for a 0.05 cutoff)
  
  #POI Numbers:
  #1 - CD8 T Cells
  #2 - CD4 T Cells
  #3 - T Cells
  #4 - Pericytess
  #5 - iCAF Fibroblasts
  #6 - Fibroblasts
  #7 - Epithelial
  #8 - Acinar
  #9 - Endothelial
  #10 - Mast Cells
  #11 - Granulocytes
  #12 - Macrophages
  #13 - MDSCs 
  #14 - Myeloid
  #15 - NK Cells
  #16 - Dendritic Cells
  #17 - Endocrine
  #18 - Perivascular
  #19 - B Cells
  
  CircosFiles <- list()
  
  CellTypes <- c("CD8TCells","CD4TCells","TCells","myCAFFibroblast","iCAFFibroblast","Fibroblasts","Epithelial","Acinar","Endothelial","MastCells",
                 "Granulocytes","Macrophages","MDSCs","Myeloid","NKCells","DendriticCells","Endocrine","Perivascular","BCells")
  
  CellTypeColors <- c("darkgreen", "limegreen","forestgreen","darkslategray3","darkcyan","darkcyan","red3","deeppink",
                      "mediumvioletred","gold","orange1","darkorange2", "tan2","darkorange2","darkorchid","chocolate4","darkred","lightskyblue","gold1")
  
  #Karyotype Function
  
  Karyotype_File <- data.frame(
    "Chr" = "chr -",
    "Name" = "Epithelial",
    "Label" = "Epithelial",
    "Start" = 0,
    "End" = 0,
    "Color" = "chr0"
  )
  
  LigList_File <- data.frame(
    "Chr" = "A",
    "Start" = 0,
    "End" = 0,
    "Name" = "A"
  )
  
  RecList_File <- data.frame(
    "Chr" = "A",
    "Start" = 0,
    "End" = 0,
    "Name" = "A"
  )
  
  DataVec_Source <- c()
  
  for(i in 1:length(CellTypeVector)){
    CellType <- CellTypeVector[i]
    DataVec_Source <- c(DataVec_Source, which(x = InteractomeData$source == CellType))
  }
  
  Data_Source <- InteractomeData[DataVec_Source,]
  
  DataVec_Target <- c()
  
  for(i in 1:length(CellTypeVector)){
    CellType <- CellTypeVector[i]
    DataVec_Target <- c(DataVec_Target, which(x = Data_Source$target == CellType))
  }
  
  Data_Target <- Data_Source[DataVec_Target,]
  
  Data <- Data_Target[which(Data_Target$source != Data_Target$target),]
  
  if(Cutoff_Lig > 0 & Cutoff_Rec > 0){
    Data <- Data[which(Data$avg_lig_0 > summary(Data$avg_lig_0)[1] | Data$avg_rec_0 > summary(Data$avg_rec_0)[Cutoff_Rec]),]
  }else{
    if(Cutoff_Lig < 0 & Cutoff_Lig < 0){
      Data <- Data[which(Data$avg_lig_0 > -1*-.01),]
      Data <- Data[which(Data$avg_rec_0 > -1*Cutoff_Rec),]
    }else{
      if(Cutoff_Lig > 0){
        Data <- Data[which(Data$avg_lig_0 > summary(Data$avg_lig_0)[Cutoff_Lig]),]
        Data <- Data[which(Data$avg_rec_0 > -1*Cutoff_Rec),]
      }else{
        Data <- Data[which(Data$avg_rec_0 > summary(Data$avg_rec_0)[Cutoff_Rec]),]
        Data <- Data[which(Data$avg_lig_0 > -1*Cutoff_Lig),]
      }
    }
  }
  
  for (i in 1:length(table(Data$source))) {
    #Determine if population is present after filtering by significance
    if( length(which(Data$source == rownames(table(Data$source))[i])) > 0 ){
      #Subset Population of Interest
      POI_Lig_Data <- Data[which(Data$source == rownames(table(Data$source))[i]),]
      POI_Rec_Data <- Data[which(Data$target == rownames(table(Data$source))[i]),]
      POI_Lig <- unique(as.character(POI_Lig_Data$ligands))
      POI_Rec <-  unique(as.character(POI_Rec_Data$receptors))
      POI_Lig_Num <- length(POI_Lig)
      POI_Rec_Num <- length(POI_Rec)
      if(POI_Lig_Num > POI_Rec_Num){
        BandSize <- POI_Lig_Num*10-1
      }else{
        BandSize <- POI_Rec_Num*10-1
      }
      InsertDF <- data.frame(
        "Chr" = "chr -",
        "Name" = as.character(rownames(table(Data$source))[i]),
        "Label" = as.character(rownames(table(Data$source))[i]),
        "Start" = 0,
        "End" = BandSize,
        "Color" = paste("chr",i,sep = "")
      ) 
      Karyotype_File <- merge(x = Karyotype_File, y = InsertDF, by = c("Chr","Name","Label","Start","End","Color"), all.y = T, all.x = T, sort = F)  
      
      if(POI_Lig_Num > POI_Rec_Num){
        
        InsertDF_LigList <- data.frame(
          "Chr" = rep(x = as.character(rownames(table(Data$source))[i]), times = POI_Lig_Num),
          "Start" = seq(from = 0, to = POI_Lig_Num*10-10, by = 10),
          "End" = seq(from = 9, to = POI_Lig_Num*10-1, by = 10),
          "Name" = POI_Lig
        ) 
        
        InsertDF_RecList <- data.frame(
          "Chr" = rep(x = as.character(rownames(table(Data$target))[i]), times = POI_Rec_Num),
          "Start" = ceiling(seq(from = 0, to = POI_Lig_Num*10-(POI_Lig_Num*10-10-0)/POI_Rec_Num, by = (POI_Lig_Num*10-1)/POI_Rec_Num)),
          "End" = floor(seq(from = (POI_Lig_Num*10-1)/POI_Rec_Num, to = POI_Lig_Num*10-1, by = (POI_Lig_Num*10-1)/POI_Rec_Num)),
          "Name" = POI_Rec
        ) 
        
      }else{
        
        InsertDF_LigList <- data.frame(
          "Chr" = rep(x = as.character(rownames(table(Data$source))[i]), times = POI_Lig_Num),
          "Start" = ceiling(seq(from = 0, to = POI_Rec_Num*10-(POI_Rec_Num*10-10-0)/POI_Lig_Num, by = (POI_Rec_Num*10-1)/POI_Lig_Num)),
          "End" = floor(seq(from = (POI_Rec_Num*10-1)/POI_Lig_Num, to = POI_Rec_Num*10-1, by = (POI_Rec_Num*10-1)/POI_Lig_Num)),
          "Name" = POI_Lig
        ) 
        
        InsertDF_RecList <- data.frame(
          "Chr" = rep(x = as.character(rownames(table(Data$target))[i]), times = POI_Rec_Num),
          "Start" = seq(from = 0, to = POI_Rec_Num*10-10, by = 10),
          "End" = seq(from = 9, to = POI_Rec_Num*10-1, by = 10),
          "Name" = POI_Rec
        ) 
      }
      
      Karyotype_File <- merge(x = Karyotype_File, y = InsertDF, by = c("Chr","Name","Label","Start","End","Color"), all.y = T, all.x = T)
      LigList_File <- merge(x = LigList_File, y = InsertDF_LigList, by = c("Chr","Start","End","Name"), all.y = T, all.x = T, sort = F)
      RecList_File <- merge(x = RecList_File, y = InsertDF_RecList, by = c("Chr","Start","End","Name"), all.y = T, all.x = T, sort = F)
      
    } 
  }
  
  Karyotype_File <- Karyotype_File[-1,]
  LigList_File <- LigList_File[-1,]
  RecList_File <- RecList_File[-1,]
  
  sentenceString <- Karyotype_File$Name
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  Karyotype_File$Name <- sentenceString
  Karyotype_File$Label <- sentenceString
  
  sentenceString <- LigList_File$Chr
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  LigList_File$Chr <- sentenceString
  
  
  sentenceString <- RecList_File$Chr
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  RecList_File$Chr <- sentenceString
  
  Band_Color <- c()
  
  for (i in 1:dim(Karyotype_File)[1]) {
    CellType <- as.character(Karyotype_File[i,"Name"])
    Band_Color <- c(Band_Color, CellTypeColors[which(CellTypes == CellType)])
  }
  Karyotype_File[,"Color"] <- Band_Color
  
  LigList <- LigList_File
  RecList <- RecList_File
  
  if(Species == T){
    LigRecList <- LR_List[,1:2]
  }else{
    LigRecList <- LR_List[,3:4]
  }
  
  
  
  
  #Expression Function
  
  LigExpression <- data.frame(
    "Chr"   = "A",
    "Name" = "A",
    "Expression" = 0
  )
  
  RecExpression <- data.frame(
    "Chr"   = "A",
    "Name" = "A",
    "Expression" = 0
  )
  
  Avg_Expression <- AverageExpression(object = SeuratObject, features = c(LigList$Name,RecList$Name))
  Avg_Expression <- Avg_Expression[[1]]
  
  for(i in 1:length(CellTypeVector)){
    
    CellSearch <- CellTypeVector[i]
    sentenceString <- CellSearch
    searchString <- ' '
    replacementString <- ''
    sentenceString = sub(searchString,replacementString,sentenceString)
    sentenceString = sub(searchString,replacementString,sentenceString)
    CellSearch <- sentenceString
    
    ident1 <- CellTypeVector[i]
    CellExpression <- subset.data.frame(x = Avg_Expression, select = ident1)
    
    LigFeatures <- as.character(LigList[which(LigList$Chr == CellSearch),"Name"])
    RecFeatures <- as.character(RecList[which(RecList$Chr == CellSearch),"Name"])
    
    Lig_DE <- data.frame(
      "Exp" = 0,
      "Name" = "a"
    )
    Rec_DE <- data.frame(
      "Exp" = 0,
      "Name" = "a"
    )
    for(k in 1:length(LigFeatures)){
      
      Lig_Exp <- subset.data.frame(x = CellExpression, subset = rownames(Avg_Expression) == LigFeatures[k])
      Lig_Exp[1,2] <- LigFeatures[k]
      colnames(Lig_Exp)[1] <- "Exp"
      colnames(Lig_Exp)[2] <- "Name"
      
      Lig_DE <- merge(x = Lig_DE, y = Lig_Exp, by = c("Exp","Name"), all.y = T, all.x = T, sort = F)
      
    }
    
    Lig_DE <- Lig_DE[-1,]
    
    for(t in 1:length(RecFeatures)){
      
      Rec_Exp <- subset.data.frame(x = CellExpression, subset = rownames(Avg_Expression) == RecFeatures[t])
      Rec_Exp[1,2] <- RecFeatures[t]
      colnames(Rec_Exp)[1] <- "Exp"
      colnames(Rec_Exp)[2] <- "Name"
      
      Rec_DE <- merge(x = Rec_DE, y = Rec_Exp, by = c("Exp","Name"), all.y = T, all.x = T, sort = F)
    }
    
    Rec_DE <- Rec_DE[-1,] 
    
    Lig_DF <- data.frame(
      "Chr" = rep(x = CellTypeVector[i], times = dim(Lig_DE)[1]),
      "Name" = Lig_DE$Name,
      "Expression" = Lig_DE$Exp
    )
    Rec_DF <- data.frame(
      "Chr" = rep(x = CellTypeVector[i], times = dim(Rec_DE)[1]),
      "Name" = Rec_DE$Name,
      "Expression" = Rec_DE$Exp
    )
    
    LigExpression <- merge(x = LigExpression, y = Lig_DF, by = c("Chr","Name","Expression"), all.y = T, all.x = T, sort = F)
    RecExpression <- merge(x = RecExpression, y = Rec_DF, by = c("Chr","Name","Expression"), all.y = T, all.x = T, sort = F)
    
    
  }
  
  LigExpression <- LigExpression[-1,]
  RecExpression <- RecExpression[-1,]
  
  sentenceString <- LigExpression$Chr
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  LigExpression$Chr <- sentenceString
  
  
  sentenceString <- RecExpression$Chr
  searchString <- ' '
  replacementString <- ''
  sentenceString = sub(searchString,replacementString,sentenceString)
  sentenceString = sub(searchString,replacementString,sentenceString)
  RecExpression$Chr <- sentenceString
  
  ExpressionOutput <- merge(x = LigExpression, y = RecExpression, by = c("Chr","Name","Expression"), all.y = T, all.x = T, sort = F)
  
  #Squeeze expression values
  ExpressionInput <- ExpressionOutput
  Normalized_Expression <- log(ExpressionOutput$Expression*1000)
  ScaledExpression <- ((Normalized_Expression - range(Normalized_Expression)[1])/
                         (range(Normalized_Expression)[2]- range(Normalized_Expression)[1]))*4
  ExpressionInput$Expression <- ScaledExpression
  
  #Text and Link Function
  
  CellNumCounter <- 1
  
  if(Lig_or_Rec == T){
    
    if(POI == 1){
      CellType <- CellTypes[1]
      CellTypeColor <- CellTypeColors[1]
    }else{
      while (POI != CellNumCounter) {
        CellNumCounter <- CellNumCounter + 1
      }
      CellType <- CellTypes[CellNumCounter]
      CellTypeColor <- CellTypeColors[CellNumCounter]
    }
    
    POI_Lig_Data <- LigList[which(LigList$Chr == CellType),]
    POI_Rec_Data <- RecList[which(RecList$Chr != CellType),]
    Text_File <- merge(x = POI_Lig_Data, y = POI_Rec_Data, by = c("Chr","Start","End","Name"), all.y = T, all.x = T, sort = F)
    
    Fil_LigRecList_Vec <- c()
    
    for(i in 1:dim(LigRecList)[1]){
      LigRecList_Lig <- as.character(LigRecList[i,1])
      LigRecList_Rec <- as.character(LigRecList[i,2])
      Lig_Present <- grep(pattern = paste("^", LigRecList_Lig,"$", sep = ""), x = POI_Lig_Data$Name)
      Rec_Present <- grep(pattern = paste("^", LigRecList_Rec,"$", sep = ""), x = POI_Rec_Data$Name)
      LigList_Row <- which(LigRecList[,1] == LigRecList_Lig)
      RecList_Row <- which(LigRecList[,2] == LigRecList_Rec)
      if(length(Lig_Present) > 0 & length(Rec_Present) > 0){
        Fil_LigRecList_Vec <- c(Fil_LigRecList_Vec, intersect(LigList_Row,RecList_Row))
      }
    }
    
    LigRecList_Fil <- LigRecList[Fil_LigRecList_Vec,]
    
    Link_File <- data.frame(
      "Chr" = "A",
      "Start" = 0,
      "End" =  0,
      "Chr1" = "A",
      "Start1" = 0,
      "End1" = 0
    )
    
    for (i in 1:length(POI_Lig_Data$Name)) {
      Rec_Data <- data.frame(
        "Chr" = "A",
        "Start" = 0,
        "End" = 0
      )
      
      Lig <- as.character(POI_Lig_Data$Name[i])
      Receptors <- as.character(LigRecList_Fil[,2][grep(x = LigRecList_Fil[,1], pattern = paste("^", Lig,"$", sep = ""),value = F)])
      
      if(length(Receptors) == 1){
        Rec_Data_Insert <- POI_Rec_Data[grep(x = POI_Rec_Data$Name, pattern = paste("^", Receptors,"$", sep = ""),value = F),1:3]
        Rec_Data <- merge(x = Rec_Data, y = Rec_Data_Insert, by = c("Chr","Start","End"), all.y = T, all.x = T, sort = F)
      }else{
        for (l in 1:length(Receptors)) {
          Receptor <- Receptors[l]
          Rec_Data_Insert <- POI_Rec_Data[grep(x = POI_Rec_Data$Name, pattern = paste("^", Receptor,"$", sep = ""),value = F),1:3]
          Rec_Data <- merge(x = Rec_Data, y = Rec_Data_Insert, by = c("Chr","Start","End"), all.y = T, all.x = T, sort = F)
        }
      }
      Rec_Data <- Rec_Data[-1,]
      
      Lig_Data <- POI_Lig_Data[i,1:3]
      
      if(length(Rec_Data$Chr) > 1){
        Lig_Data <- rbind(Lig_Data, Lig_Data[rep(1, length(Rec_Data$Chr)-1), ])
      }
      
      Lig_Rec_Merge <- bind_cols(x = Lig_Data, y = Rec_Data)
      
      Link_File <- rbind(x = Link_File, y = Lig_Rec_Merge)
    }
    
    
    
    Link_File <- Link_File[-1,]
    Link_File[,(dim(Link_File)[2]+1)] <-rep(x = paste("color=",CellTypeColor,sep = ""), times = dim(Link_File)[1])
    
    Expression_File <- Text_File
    Expression_File[,5] <- rep(0, times = dim(Expression_File)[1])
    
    for(k in 1:dim(Expression_File)[1]){
      GeneCell <- Expression_File[k,"Chr"]
      Gene <- as.character(Expression_File[k,"Name"])
      ExpressionPull <- ExpressionInput[which(ExpressionInput$Chr == GeneCell & ExpressionInput$Name == Gene),]
      Expression_File[k,5] <- ExpressionPull$Expression
    }
    
    Expression_File <- Expression_File[,-4]
    
  }else{
    
    if(POI == 1){
      CellType <- CellTypes[1]
    }else{
      while (POI != CellNumCounter) {
        CellNumCounter <- CellNumCounter + 1
      }
      CellType <- CellTypes[CellNumCounter]
      CellTypeColor <- CellTypeColors[CellNumCounter]
    }
    
    POI_Lig_Data <- LigList[which(LigList$Chr != CellType),]
    POI_Rec_Data <- RecList[which(RecList$Chr == CellType),]
    Text_File <- merge(x = POI_Rec_Data, y = POI_Lig_Data, by = c("Chr","Start","End","Name"), all.y = T, all.x = T, sort = F)
    
    Fil_LigRecList_Vec <- c()
    
    for(i in 1:dim(LigRecList)[1]){
      LigRecList_Lig <- as.character(LigRecList[i,1])
      LigRecList_Rec <- as.character(LigRecList[i,2])
      Lig_Present <- grep(pattern = paste("^", LigRecList_Lig,"$", sep = ""), x = POI_Lig_Data$Name)
      Rec_Present <- grep(pattern = paste("^", LigRecList_Rec,"$", sep = ""), x = POI_Rec_Data$Name)
      LigList_Row <- which(LigRecList[,1] == LigRecList_Lig)
      RecList_Row <- which(LigRecList[,2] == LigRecList_Rec)
      if(length(Lig_Present) > 0 & length(Rec_Present) > 0){
        Fil_LigRecList_Vec <- c(Fil_LigRecList_Vec, intersect(LigList_Row,RecList_Row))
      }
    }
    
    LigRecList_Fil <- LigRecList[Fil_LigRecList_Vec,]    
    
    Link_File <- data.frame(
      "Chr" = "A",
      "Start" = 0,
      "End" =  0,
      "Chr1" = "A",
      "Start1" = 0,
      "End1" = 0
    )
    
    
    
    for (i in 1:length(POI_Rec_Data$Name)) {
      Lig_Data <- data.frame(
        "Chr" = "A",
        "Start" = 0,
        "End" = 0
      )
      
      Rec <- as.character(POI_Rec_Data$Name[i])
      Ligands <- as.character(LigRecList_Fil[,1][grep(x = LigRecList_Fil[,2], pattern = paste("^", Rec,"$", sep = ""),value = F)])
      
      if(length(Ligands == 1)){
        Lig_Data_Insert <- POI_Lig_Data[grep(x = POI_Lig_Data$Name, pattern = paste("^", Ligands,"$", sep = ""),value = F),1:3]
        Lig_Data <- merge(x = Lig_Data, y = Lig_Data_Insert, by = c("Chr","Start","End"), all.y = T, all.x = T, sort = F)
      }else{
        for (l in 1:length(Ligands)) {
          Ligand <- Ligands[l]
          Lig_Data_Insert <- POI_Lig_Data[grep(x = POI_Lig_Data$Name, pattern = paste("^", Ligand,"$", sep = ""),value = F),1:3]
          Lig_Data <- merge(x = Lig_Data, y = Lig_Data_Insert, by = c("Chr","Start","End"), all.y = T, all.x = T, sort = F)
        }
      }
      Lig_Data <- Lig_Data[-1,]
      
      Rec_Data <- POI_Rec_Data[i,1:3]
      
      if(length(Lig_Data$Chr) > 1){
        Rec_Data <- rbind(Rec_Data, Rec_Data[rep(1, length(Lig_Data$Chr)-1), ])
      }
      Lig_Rec_Merge <- bind_cols(x = Rec_Data, y = Lig_Data)
      
      Link_File <- rbind(x = Link_File, y = Lig_Rec_Merge)
    }
    
    
    
    Link_File <- Link_File[-1,]
    
    Link_Color <- c()
    
    for (i in 1:dim(Link_File)[1]) {
      CellType <- as.character(Link_File[i,"Chr1"])
      Link_Color <- c(Link_Color, CellTypeColors[which(CellTypes == CellType)])
    }
    Link_Color <- paste("color=",Link_Color,"_a4",sep = "")
    Link_File[,(dim(Link_File)[2]+1)] <- Link_Color
    
    Expression_File <- Text_File
    Expression_File[,5] <- rep(0, times = dim(Expression_File)[1])
    
    for(k in 1:dim(Expression_File)[1]){
      GeneCell <- Expression_File[k,"Chr"]
      Gene <- as.character(Expression_File[k,"Name"])
      ExpressionPull <- ExpressionInput[which(ExpressionInput$Chr == GeneCell & ExpressionInput$Name == Gene),]
      Expression_File[k,5] <- ExpressionPull$Expression
    }
    
    Expression_File <- Expression_File[,-4]
    
  }
  
  CircosFiles[[1]] <- Karyotype_File
  CircosFiles[[2]] <- Text_File
  CircosFiles[[3]] <- Expression_File
  CircosFiles[[4]] <- Link_File 
  
  return(CircosFiles)
}


#____________________________________________________________________________________________________________________________________________________________________#

# 2. Data Pre-processing

#Tissue Object

#Load in all filtered runs using Read10X or Read10X_h5 functions
#   Data.data <- Read10X("filepath/")/Read10X_h5("filepath/file.ext")

PDAC_TISSUE_1.data <- Read10X("~/Desktop/PDAC_TISSUE_1/filtered_feature_bc_matrix/")
PDAC_TISSUE_2.data <- Read10X("~/Desktop/PDAC_TISSUE_2/filtered_feature_bc_matrix/")
PDAC_TISSUE_3.data <- Read10X("~/Desktop/PDAC_TISSUE_3/filtered_feature_bc_matrix/")
PDAC_TISSUE_4.data <- Read10X("~/Desktop/PDAC_TISSUE_4/filtered_feature_bc_matrix/")
PDAC_TISSUE_5.data <- Read10X("~/Desktop/PDAC_TISSUE_5/filtered_feature_bc_matrix/")
PDAC_TISSUE_6.data <- Read10X("~/Desktop/PDAC_TISSUE_6/filtered_feature_bc_matrix/")
PDAC_TISSUE_7.data <- Read10X("~/Desktop/PDAC_TISSUE_7/filtered_feature_bc_matrix/")
PDAC_TISSUE_8.data <- Read10X("~/Desktop/PDAC_TISSUE_8/filtered_feature_bc_matrix/")
PDAC_TISSUE_9.data <- Read10X("~/Desktop/PDAC_TISSUE_9/filtered_feature_bc_matrix/")
PDAC_TISSUE_10.data <- Read10X("~/Desktop/PDAC_TISSUE_10/filtered_feature_bc_matrix/")
PDAC_TISSUE_11A.data <- Read10X("~/Desktop/PDAC_TISSUE_11A/filtered_feature_bc_matrix/")
PDAC_TISSUE_11B.data <- Read10X("~/Desktop/PDAC_TISSUE_11B/filtered_feature_bc_matrix/")
PDAC_TISSUE_12.data <- Read10X("~/Desktop/PDAC_TISSUE_12/filtered_feature_bc_matrix/")
PDAC_TISSUE_13.data <- Read10X("~/Desktop/PDAC_TISSUE_13/filtered_feature_bc_matrix/")
PDAC_TISSUE_14.data <- Read10X_h5("~/Desktop/PDAC_TISSUE_14/filtered_feature_bc_matrix.h5")
PDAC_TISSUE_15.data <- Read10X("~/Desktop/PDAC_TISSUE_15/filtered_feature_bc_matrix/")
PDAC_TISSUE_16.data <- Read10X_h5("~/Desktop/PDAC_TISSUE_16/filtered_gene_bc_matrices_h5.h5")

AdjNorm_TISSUE_1.data <- Read10X("~/Desktop/AdjNorm_TISSUE_1/filtered_feature_bc_matrix/")
AdjNorm_TISSUE_2.data <- Read10X("~/Desktop/AdjNorm_TISSUE_2/filtered_feature_bc_matrix/")
AdjNorm_TISSUE_3.data <- Read10X("~/Desktop/AdjNorm_TISSUE_3/filtered_feature_bc_matrix/")

#Generate Seurat objects using CreateSeuratObject function, min.cells=3, min.features=100
# Data <- CreateSeuratObject(counts = Data.data, project = "PDA", min.cells=3, min.features=100)

PDAC_TISSUE_1 <- CreateSeuratObject(counts = PDAC_TISSUE_1.data, project = 'PDAC_TISSUE_1', min.cells = 3, min.features = 100)
PDAC_TISSUE_2 <- CreateSeuratObject(counts = PDAC_TISSUE_2.data, project = 'PDAC_TISSUE_2', min.cells = 3, min.features = 100)
PDAC_TISSUE_3 <- CreateSeuratObject(counts = PDAC_TISSUE_3.data, project = 'PDAC_TISSUE_3', min.cells = 3, min.features = 100)
PDAC_TISSUE_4 <- CreateSeuratObject(counts = PDAC_TISSUE_4.data, project = 'PDAC_TISSUE_4', min.cells = 3, min.features = 100)
PDAC_TISSUE_5 <- CreateSeuratObject(counts = PDAC_TISSUE_5.data, project = 'PDAC_TISSUE_5', min.cells = 3, min.features = 100)
PDAC_TISSUE_6 <- CreateSeuratObject(counts = PDAC_TISSUE_6.data, project = 'PDAC_TISSUE_6', min.cells = 3, min.features = 100)
PDAC_TISSUE_7 <- CreateSeuratObject(counts = PDAC_TISSUE_7.data, project = 'PDAC_TISSUE_7', min.cells = 3, min.features = 100)
PDAC_TISSUE_8 <- CreateSeuratObject(counts = PDAC_TISSUE_8.data, project = 'PDAC_TISSUE_8', min.cells = 3, min.features = 100)
PDAC_TISSUE_9 <- CreateSeuratObject(counts = PDAC_TISSUE_9.data, project = 'PDAC_TISSUE_9', min.cells = 3, min.features = 100)
PDAC_TISSUE_10 <- CreateSeuratObject(counts = PDAC_TISSUE_10.data, project = 'PDAC_TISSUE_10', min.cells = 3, min.features = 100)
PDAC_TISSUE_11A <- CreateSeuratObject(counts = PDAC_TISSUE_11A.data, project = 'PDAC_TISSUE_11A', min.cells = 3, min.features = 100)
PDAC_TISSUE_11B <- CreateSeuratObject(counts = PDAC_TISSUE_11B.data, project = 'PDAC_TISSUE_11B', min.cells = 3, min.features = 100)
PDAC_TISSUE_12 <- CreateSeuratObject(counts = PDAC_TISSUE_12.data, project = 'PDAC_TISSUE_12', min.cells = 3, min.features = 100)
PDAC_TISSUE_13 <- CreateSeuratObject(counts = PDAC_TISSUE_13.data, project = 'PDAC_TISSUE_13', min.cells = 3, min.features = 100)
PDAC_TISSUE_14 <- CreateSeuratObject(counts = PDAC_TISSUE_14.data, project = 'PDAC_TISSUE_14', min.cells = 3, min.features = 100)
PDAC_TISSUE_15 <- CreateSeuratObject(counts = PDAC_TISSUE_15.data, project = 'PDAC_TISSUE_15', min.cells = 3, min.features = 100)
PDAC_TISSUE_16 <- CreateSeuratObject(counts = PDAC_TISSUE_16.data, project = 'PDAC_TISSUE_16', min.cells = 3, min.features = 100)

AdjNorm_TISSUE_1 <- CreateSeuratObject(counts = AdjNorm_TISSUE_1.data, project = 'AdjNorm_TISSUE_1', min.cells = 3, min.features = 100)
AdjNorm_TISSUE_2 <- CreateSeuratObject(counts = AdjNorm_TISSUE_2.data, project = 'AdjNorm_TISSUE_2', min.cells = 3, min.features = 100)
AdjNorm_TISSUE_3 <- CreateSeuratObject(counts = AdjNorm_TISSUE_3.data, project = 'AdjNorm_TISSUE_3', min.cells = 3, min.features = 100)

#Add Metadata

#ID Metadata
PDAC_TISSUE_16$ID <- "PDAC_TISSUE_16"
PDAC_TISSUE_15$ID <- "PDAC_TISSUE_15"
PDAC_TISSUE_13$ID <- "PDAC_TISSUE_13"
PDAC_TISSUE_14$ID <- "PDAC_TISSUE_14"
PDAC_TISSUE_12$ID <- "PDAC_TISSUE_12"
PDAC_TISSUE_11B$ID <- "PDAC_TISSUE_11B" 
PDAC_TISSUE_11A$ID <- "PDAC_TISSUE_11A" 
PDAC_TISSUE_10$ID <- "PDAC_TISSUE_10"
PDAC_TISSUE_9$ID <- "PDAC_TISSUE_9"
PDAC_TISSUE_8$ID <- "PDAC_TISSUE_8"
PDAC_TISSUE_7$ID <- "PDAC_TISSUE_7"
PDAC_TISSUE_6$ID <- "PDAC_TISSUE_6"
PDAC_TISSUE_5$ID <- "PDAC_TISSUE_5"
PDAC_TISSUE_4$ID <- "PDAC_TISSUE_4"
PDAC_TISSUE_3$ID <- "PDAC_TISSUE_3"
PDAC_TISSUE_2$ID <- "PDAC_TISSUE_2"
PDAC_TISSUE_1$ID <- "PDAC_TISSUE_1"
AdjNorm_TISSUE_1$ID <-"AdjNorm_TISSUE_1"
AdjNorm_TISSUE_3$ID <- "AdjNorm_TISSUE_3"
AdjNorm_TISSUE_2$ID <- "AdjNorm_TISSUE_2"

#Disease State Metadata
PDAC_TISSUE_16$DiseaseState <-"PDAC"
PDAC_TISSUE_15$DiseaseState <-"PDAC"
PDAC_TISSUE_13$DiseaseState <-"PDAC"
PDAC_TISSUE_14$DiseaseState <-"PDAC"
PDAC_TISSUE_12$DiseaseState <-"PDAC"
PDAC_TISSUE_11B$DiseaseState <-"PDAC"  
PDAC_TISSUE_11A$DiseaseState <-"PDAC"
PDAC_TISSUE_10$DiseaseState <-"PDAC"
PDAC_TISSUE_9$DiseaseState <-"PDAC"
PDAC_TISSUE_8$DiseaseState <-"PDAC"
PDAC_TISSUE_7$DiseaseState <-"PDAC"
PDAC_TISSUE_6$DiseaseState <-"PDAC"
PDAC_TISSUE_5$DiseaseState <-"PDAC"
PDAC_TISSUE_4$DiseaseState <-"PDAC"
PDAC_TISSUE_3$DiseaseState <-"PDAC"
PDAC_TISSUE_2$DiseaseState <-"PDAC"
PDAC_TISSUE_1$DiseaseState <-"PDAC"
AdjNorm_TISSUE_1$DiseaseState <- "AdjNorm"
AdjNorm_TISSUE_3$DiseaseState <- "AdjNorm"
AdjNorm_TISSUE_2$DiseaseState <- "AdjNorm"

#Merge all objects
TotalTissue.combined <- merge(PDAC_TISSUE_1, y = c(PDAC_TISSUE_2, PDAC_TISSUE_3, PDAC_TISSUE_4, 
                                                   PDAC_TISSUE_5, PDAC_TISSUE_6, PDAC_TISSUE_7, 
                                                   PDAC_TISSUE_8, PDAC_TISSUE_9, PDAC_TISSUE_10, 
                                                   PDAC_TISSUE_11A, PDAC_TISSUE_11B, PDAC_TISSUE_12,
                                                   PDAC_TISSUE_13, PDAC_TISSUE_14, PDAC_TISSUE_15, 
                                                   PDAC_TISSUE_16))
#Check all objects present
table(TotalTissue.combined$orig.ident)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata
Idents(object = TotalTissue.combined) <- 'ID'
# Check active identity
levels(TotalTissue.combined)

#NORMALIZE DATA
TotalTissue.combined <- NormalizeData(object = TotalTissue.combined, normalization.method = "LogNormalize", 
                                      scale.factor = 10000)
#Percent Mitochondrial Genes
#QC Metric used to remove cells with overabundant Mitochondrial genes, typically associated with nuclear wash out during sequencing
TotalTissue.combined[["percent.mt"]] <- PercentageFeatureSet(TotalTissue.combined, pattern = "^MT-")

#Plot the nFeatures/counts/% Mito to get general idea about the quality of your data
VlnPlot(TotalTissue.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = .0)

#FIND VARIABLE GENES
TotalTissue.combined<- FindVariableFeatures(object = TotalTissue.combined, mean.function = ExpMean, dispersion.function = LogVMR, 
                                            x.low.cutoff = 0.0125, y.cutoff = 0.5)


#Calculate Cell Cycle Score (S-G2M Difference)
#s.genes <- cc.genes$s.genes
#g2m.genes <- cc.genes$g2m.genes
#TotalTissue.combined<- CellCycleScoring(TotalTissue.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#TotalTissue.combined$CC.Difference <- TotalTissue.combined$S.Score - TotalTissue.combined$G2M.Score

#Clear up memory
rm(PDAC_TISSUE_1.data, PDAC_TISSUE_2.data, PDAC_TISSUE_3.data, PDAC_TISSUE_4.data, PDAC_TISSUE_5.data, PDAC_TISSUE_6.data, PDAC_TISSUE_7.data, 
   PDAC_TISSUE_8.data, PDAC_TISSUE_9.data, PDAC_TISSUE_10.data, PDAC_TISSUE_11A.data, PDAC_TISSUE_11B.data, PDAC_TISSUE_12.data,
   PDAC_TISSUE_13.data, PDAC_TISSUE_14.data, PDAC_TISSUE_15.data, PDAC_TISSUE_16.data)
rm(PDAC_TISSUE_1, PDAC_TISSUE_2, PDAC_TISSUE_3, PDAC_TISSUE_4, PDAC_TISSUE_5, PDAC_TISSUE_6, PDAC_TISSUE_7, 
PDAC_TISSUE_8, PDAC_TISSUE_9, PDAC_TISSUE_10, PDAC_TISSUE_11A, PDAC_TISSUE_11B, PDAC_TISSUE_12,
PDAC_TISSUE_13, PDAC_TISSUE_14, PDAC_TISSUE_15, PDAC_TISSUE_16)

#THIS STEP MAY TAKE A VERY LONG TIME
#Scale Data
TotalTissue.combined<- ScaleData(object = TotalTissue.combined, vars.to.regress = "nCount_RNA", features = rownames(TotalTissue.combined))

#It is recommended that you save this object so you do not have to re-run the scaling step more than necessary
save(TotalTissue.combined,file="TotalTissue.combined.RData")
#____________________________________________________________________________________________________________________________________________________________________#

#Data Visualization

#Run PCA and Determine Dimensions for 90% Variance
TotalTissue.combined <- RunPCA(object = TotalTissue.combined, features = VariableFeatures(object = TotalTissue.combined))
stdev <- TotalTissue.combined@reductions$pca@stdev
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)


#Find Neighbors + Find CLusters (without harmony batch correction)
TotalTissue.combined <- FindNeighbors(object = TotalTissue.combined, dims = 1:PCNum)
TotalTissue.combined <- FindClusters(object = TotalTissue.combined, resolution = 1.2)

#Run UMAP and get unlabelled cluster UMAP and violin plot (without harmony batch correction)
TotalTissue.combined <- RunUMAP(object = TotalTissue.combined, dims = 1:PCNum)
DimPlot(object = TotalTissue.combined, reduction = "umap", label = TRUE, pt.size = 0.5)
TotalTissue.combined[["UMAP_Clusters"]] <- Idents(object = TotalTissue.combined)

#UMAP split by group
DimPlot(object = TotalTissue.combined, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "ID", ncol=5)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata - will need to do this to make the right plots
Idents(object = TotalTissue.combined) <- 'UMAP_Clusters'
# Check active identity
levels(TotalTissue.combined)

#Umap with overlaid group
Idents(object = TotalTissue.combined) <- 'ID'
DimPlot(object = TotalTissue.combined, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(TotalTissue.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Find neighbors and clusters WITH harmony batch correction
options(repr.plot.height = 2.5, repr.plot.width = 6)
TotalTissue.combined <- TotalTissue.combined %>% 
  RunHarmony("ID", plot_convergence = TRUE)

TotalTissue.combined.harmony <- FindNeighbors(object = TotalTissue.combined, dims = 1:PCNum, reduction ="harmony")
TotalTissue.combined.harmony <- FindClusters(object = TotalTissue.combined, resolution = 1.2, reduction ="harmony")

#Run UMAP and get unlabelled cluster UMAP and violin plot
TotalTissue.combined.harmony <- RunUMAP(object = TotalTissue.combined.harmony, dims = 1:PCNum, reduction = "harmony")
DimPlot(object = TotalTissue.combined.harmony, reduction = "umap", label = TRUE, pt.size = 0.5)

#UMAP split by group
DimPlot(object = TotalTissue.combined.harmony, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "ID", ncol=5)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata - will need to do this to make the right plots
Idents(object = TotalTissue.combined.harmony) <- 'group'
# Check active identity
levels(TotalTissue.combined.harmony)

#Umap with overlaid group
Idents(object = TotalTissue.combined.harmony) <- 'ID'
DimPlot(object = TotalTissue.combined.harmony, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(TotalTissue.combined.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Make the active identity UMAP clusters again
Idents(object = TotalTissue.combined.harmony) <- 'UMAP_Clusters'
# Check active identity
levels(TotalTissue.combined.harmony)

#MAKE CLUSTER MARKER TABLE

# This function automates the FindMarkers function and uses the list of markers to broadly
# identify cell types based on a preselected list of markers. Markers chosen for human samples.
# The output is a dataframe containing the FindMarkers output.

ClusterMarkerTable <- function(ClustNum,Data,NumMarkers,PCT) {
  ClusterMarks <- FindMarkers(object = Data, ident.1 = ClustNum[1], min.pct = PCT)
  ClusterTable <- as.data.frame(head(x = ClusterMarks, n = NumMarkers))
  MarkerLabelRow <- length(ClusterTable$p_val)+1
  MarkerLabel <- c("Cluster", ClustNum[1], "Markers", "Are", "Above")
  ClusterTable[MarkerLabelRow,] <- MarkerLabel
  
  if(length(ClustNum) > 2) {
    for(i in ClustNum[2:length(ClustNum)]){
      ClusterMarks <- FindMarkers(object = Data, ident.1 = i, min.pct = PCT)
      ClusterTable <- rbind(ClusterTable,head(x = ClusterMarks, n = NumMarkers))
      MarkerLabelRow <- length(ClusterTable$p_val)+1
      MarkerLabel <- c("Cluster", i, "Markers", "Are", "Above")
      ClusterTable[MarkerLabelRow,] <- MarkerLabel
    }
  }
  View(ClusterTable)
  return(ClusterTable)
}

#Clear up memory
rm(TotalTissue.combined)
rm(AdjNorm_TISSUE_1, AdjNorm_TISSUE_1.data, AdjNorm_TISSUE_2, AdjNorm_TISSUE_2.data, AdjNorm_TISSUE_3, AdjNorm_TISSUE_3.data)

#Save data
save(TotalTissue.combined.harmony,file="TotalTissue.combined.harmony.RData")

