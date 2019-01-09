#' Finding ligand-receptor pairs
#'
#' This function loads the highly expressed genes or differentail expressed
#' genes as a dataframe. Significant interactions are found through mapping
#' these genes to our ligand-receptor database.
#'
#' @param data_1 Data used to find the ligand-receptor pairs
#' @param data_2 Second dataset used to find ligand-receptor pairs. If set NULL,
#' paris will be found within data_1. Otherwise, pairs will be found between
#' data_1 and data_2. Default is NULL.
#' @param datatype Type of data used as input. Options are "mean count"
#' and "DEG"
#' @param comm_type Communication type. Available options are "cytokine",
#' "checkpoint", "growth factor", "other"
#' @param database Database used to find ligand-receptor pairs. If set NULL,
#' the build-in database will be used.
#' @import dplyr
#' @references Cytokines, Inflammation and Pain. Zhang et al,2007.
#' @references Cytokines, Chemokines and Their Receptors. Cameron et al, 2000-2013
#' @references Robust prediction of response to immune checkpoint blockade therapy
#' in metastatic melanoma. Auslander et al, 2018.
#' @references A draft network of ligand-receptor-mediated multicellular signalling
#' in human, Jordan A. Ramilowski, Nature Communications, 2015
#' @return A dataframe of the significant interactions
#' @export
FindLR<-function(data_1,data_2=NULL,datatype,comm_type,database=NULL){
  if(is.null(database)){
    data('LR_database')
    database<-db
  }
  database<-database[database$Classification==comm_type,]
  if(datatype=='mean count'){
    gene_list_1<-data_1
    if(is.null(data_2)){
      gene_list_2<-gene_list_1
    }else{
      gene_list_2<-data_2
    }
    ligand_ind<-which(database$Ligand.ApprovedSymbol %in% gene_list_1$gene)
    receptor_ind<-which(database$Receptor.ApprovedSymbol %in% gene_list_2$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    FilterTable_1<-database[ind,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
      left_join(gene_list_1[,c('gene','exprs','cell_type')],by=c('Ligand.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_from_mean_exprs=mean_exprs,cell_from=cell_type) %>%
      left_join(gene_list_2[,c('gene','exprs','cell_type')],by=c('Receptor.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_to_mean_exprs=mean_exprs,cell_to=cell_type)
    ligand_ind<-which(database$Ligand.ApprovedSymbol %in% gene_list_2$gene)
    receptor_ind<-which(database$Receptor.ApprovedSymbol %in% gene_list_1$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    FilterTable_2<-database[ind,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
      left_join(gene_list_2[,c('gene','exprs','cell_type')],by=c('Ligand.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_from_mean_exprs=mean_exprs,cell_from=cell_type) %>%
      left_join(gene_list_1[,c('gene','exprs','cell_type')],by=c('Receptor.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_to_mean_exprs=mean_exprs,cell_to=cell_type)
    FilterTable<-rbind(FilterTable_1,FilterTable_2)
  }else if(datatype=='DEG'){
    gene_list_1<-data_1
    if(is.null(data_2)){
      gene_list_2<-gene_list_1
    }else{
      gene_list_2<-data_2
    }
    ligand_ind<-which(database$Ligand.ApprovedSymbol %in% gene_list_1$gene)
    receptor_ind<-which(database$Receptor.ApprovedSymbol %in% gene_list_2$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    FilterTable_1<-database[ind,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
      left_join(gene_list_1[,c('gene','logFC','q.value','cell_type')],by=c('Ligand.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_from_logFC=logFC,cell_from_q.value=q.value,cell_from=cell_type) %>%
      left_join(gene_list_2[,c('gene','logFC','q.value','cell_type')],by=c('Receptor.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_to_logFC=logFC,cell_to_q.value=q.value,cell_to=cell_type)
    ligand_ind<-which(database$Ligand.ApprovedSymbol %in% gene_list_2$gene)
    receptor_ind<-which(database$Receptor.ApprovedSymbol %in% gene_list_1$gene)
    ind<-intersect(ligand_ind,receptor_ind)
    FilterTable_2<-database[ind,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')] %>%
      left_join(gene_list_2[,c('gene','logFC','q.value','cell_type')],by=c('Ligand.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_from_logFC=logFC,cell_from_q.value=q.value,cell_from=cell_type) %>%
      left_join(gene_list_1[,c('gene','logFC','q.value','cell_type')],by=c('Receptor.ApprovedSymbol'='gene')) %>%
      dplyr::rename(cell_to_logFC=logFC,cell_to_q.value=q.value,cell_to=cell_type)
    FilterTable<-rbind(FilterTable_1,FilterTable_2)
  }else{
    stop('Error: invalid data type')
  }

  FilterTable<-FilterTable[!duplicated(FilterTable),]
  res<-as.data.frame(FilterTable) %>% dplyr::rename(ligand=Ligand.ApprovedSymbol,receptor=Receptor.ApprovedSymbol)
  if(datatype=='DEG'){
    res<-res[!(res$cell_from_logFC==0.0001 & res$cell_to_logFC==0.0001),]
  }
  res<-res %>% mutate(comm_type=comm_type)
  return(res)
}
