#' Plotting ligand-receptor pairs
#'
#' This function loads the significant interactions as a dataframe. A circle
#' plot will be generated using package circlize. The width of the arrow
#' represents the expression level/log fold change of the ligand; while the
#' width of arrow head represents the expression level/log fold change of the
#' receptor. Different color and the type of the arrow stands for whether
#' the ligand and/or receptor are upregulated or downregulated. Users can select
#' the colors represent the cell type by their own or chosen randomly by default.
#'
#' @references Gu, Z. (2014) circlize implements and enhances circular
#' visualization in R. Bioinformatics.
#' @param data A dataframe contains significant ligand-receptor pairs and related
#' information such as expression level/log fold change and cell type
#' @param datatype Type of data. Options are "mean count" and "DEG"
#' @param gene_col Colors used to represent different categories of genes.
#' @param transparency Transparency of link colors, 0 means no transparency and
#' 1 means full transparency. If transparency is already set in col or row.col
#' or column.col, this argument will be ignored. NAalso ignores this argument.
#' @param link.arr.lwd line width of the single line link which is put in the
#' center of the belt.
#' @param link.arr.lty line type of the single line link which is put in the
#' center of the belt.
#' @param link.arr.col color or the single line link which is put in the center
#' of the belt.
#' @param link.arr.width size of the single arrow head link which is put in the
#' center of the belt.
#' @param link.arr.type Type of the arrows, pass to Arrowhead. Default value is
#' triangle. There is an additional option big.arrow
#' @param facing Facing of text.
#' @param cell_col Colors used to represent types of cells. If set NULL, it
#' will be generated randomly
#' @param print.cell Whether or not print the type of cells on the outer layer
#' of the graph.
#' @param track.height_1 height of the cell notation track
#' @param track.height_2 height of the gene notation track
#' @param annotation.height_1 Track height corresponding to values in annotationTrack.
#' @param annotation.height_2 Track height corresponding to values in annotationTrack.
#' @param text.vjust adjustment on ’vertical’ (radical) direction. Besides to set it
#' as numeric values, the value can also be a string contain absoute unit, e.g.
#' "2.1mm", "-1 inche", but only "mm", "cm", "inches"/"inche" are allowed.
#' @import randomcoloR
#' @import graphics
#' @import circlize
#'
#' @return A figure of the significant interactions
#' @export
LRPlot<-function(data,datatype,gene_col=NULL,transparency=0.5,link.arr.lwd=1,link.arr.lty=NULL,link.arr.col=NULL,link.arr.width=NULL,
                 link.arr.type=NULL,facing='clockwise',cell_col=NULL,print.cell=TRUE,track.height_1=uh(2,'mm'),track.height_2=uh(12,'mm'),
                 annotation.height_1=0.01,annotation.height_2=0.01,text.vjust = '0.4cm',...){
  cell_group<-unique(c(data$cell_from,data$cell_to))
  genes<-c(structure(data$ligand,names=data$cell_from),structure(data$receptor,names=data$cell_to))
  genes<-genes[!duplicated(paste(names(genes),genes))]
  genes<-genes[order(names(genes))]
  if(is.null(link.arr.lty)){
    if(datatype=='mean count'){
      link.arr.lty='solid'
    }else if(datatype=='DEG'){
      link.arr.lty=structure(ifelse(data$cell_from_logFC==0.0001,'dashed','solid'),names=paste(data$cell_from,data$receptor))
    }else{
      print('invalid datatype')
    }
  }
  if(is.null(link.arr.col)){
    if(datatype=='mean count'){
      data<-data %>% mutate(link_col='black')
    }else if(datatype=='DEG'){
      data<-data %>% mutate(link_col=ifelse(cell_from_logFC==0.0001,ifelse(cell_to_logFC>0,'#d73027','#00ccff'),
                                            ifelse(cell_to_logFC==0.0001,ifelse(cell_from_logFC>0,'#d73027','#00ccff'),
                                                   ifelse(cell_from_logFC>0,ifelse(cell_to_logFC>0,'#d73027','#dfc27d'),
                                                          ifelse(cell_to_logFC>0,'#9933ff','#00ccff')))))
    }else{
      print('invalid datatype')
    }
  }else{
    data$link_col=link.arr.col
  }
  if(is.null(link.arr.type)){
    if(datatype=='mean count'){
      link.arr.type='triangle'
    }else if(datatype=='DEG'){
      link.arr.type=structure(ifelse(data$cell_to_logFC==0.0001,'ellipse','triangle'),names=paste(data$cell_from,data$receptor))
    }else{
      print('invalid datatype')
    }
  }
  if(is.null(gene_col)){
    comm_col<-structure(c('#99ff99','#99ccff','#ff9999','#ffcc99'),names=c('other','cytokine','checkpoint','growth factor'))
    gene_col<-structure(c(comm_col[data$comm_type],rep('#073c53',length(data$receptor))),names=c(data$ligand,data$receptor))
  }
  if(is.null(cell_col)){
    cell_col<-structure(randomColor(count=length(unique(names(genes))),luminosity='dark'),names=unique(names(genes)))
  }
  if(is.null(link.arr.lwd)){
    data<-data %>% mutate(arr_width=1)
  }else if(max(abs(link.arr.lwd))-min(abs(link.arr.lwd))==0 && all(link.arr.lwd!=0.0001)){
    data<-data %>% mutate(arr_width=ifelse(abs(link.arr.lwd<5),abs(link.arr.lwd),5))
  }else{
    data<-data %>% mutate(arr_width=ifelse(link.arr.lwd==0.0001,2,1+5/(max(abs(link.arr.lwd))-min(abs(link.arr.lwd)))*(abs(link.arr.lwd)-min(abs(link.arr.lwd)))))
  }
  if(length(cell_group)!=1){
    gap.degree <- do.call("c", lapply(table(names(genes)), function(i) c(rep(1, i-1), 8)))
  }else{
    gap.degree <- do.call("c", lapply(table(names(genes)), function(i) c(rep(1, i))))
  }
  circos.par(gap.degree = gap.degree)
  if(length(gene_col)==1){
    grid.col=gene_col
  }else{
    grid.col=gene_col[genes]
    names(grid.col)<-paste(names(genes),genes)
  }
  if(is.null(link.arr.width)){
    data<-data %>% mutate(link.arr.width=data$arr_width/10)
  }else if(max(abs(link.arr.width))-min(abs(link.arr.width))==0 && all(link.arr.width!=0.0001)){
    data<-data %>% mutate(link.arr.width=ifelse(abs(link.arr.width)<0.5,abs(link.arr.width),0.5))
  }else{
    data<-data %>% mutate(link.arr.width=ifelse(link.arr.width==0.0001,0.2,(1+5/(max(abs(link.arr.width))-min(abs(link.arr.width)))*(abs(link.arr.width)-min(abs(link.arr.width))))/10))
  }
  chordDiagram(as.data.frame(cbind(paste(data$cell_from,data$ligand),paste(data$cell_to,data$receptor))), order=paste(names(genes),genes),
               grid.col=grid.col,transparency=transparency,directional=1,direction.type='arrows',link.arr.lwd=data$arr_width,link.arr.lty=link.arr.lty,
               link.arr.type=link.arr.type,link.arr.width=data$link.arr.width,link.arr.col=data$link_col,col='#00000000',annotationTrack=c('grid'),preAllocateTracks = list(
                 list(track.height = track.height_1),list(track.height = track.height_2)),annotationTrackHeight = c(annotation.height_1,annotation.height_2),...)

  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = genes[get.cell.meta.data("sector.numeric.index")]
    circos.text(mean(xlim),mean(ylim),sector.index, col = "black", cex = 0.7, facing = facing, niceFacing = TRUE)
  }, bg.border = 0)

  if(print.cell){
    for(c in unique(names(genes))) {
      gene = as.character(genes[names(genes) == c])
      highlight.sector(sector.index = paste(c,gene), track.index = 1, col = ifelse(length(cell_col)==1,cell_col,cell_col[c]), text = c, text.vjust = text.vjust, niceFacing = TRUE,lwd=1)
    }
  }
  circos.clear()
}
