#' Plotting ligand-receptor pairs
#'
#' This function loads count data as dataframe, ligand, receptor and two interactive
#' cells' names as strings. The plot shows the expression level of ligand and
#' receptor at different time, thus illustrates a dynamic change of a ligand-receptor
#' pairs.
#'
#' @param data A dataframe contains significant ligand-receptor pairs and related
#' information such as expression level/log fold change and cell type
#' @param ligand String as selected ligand
#' @param receptor String as selected receptor
#' @param cell_from The cell type ligand gene belongs to
#' @param cell_to The cell type receptor gene belongs to
#' @param Time Different time points showing on the plot
#' @import tidyr
#' @import ggplot2
#' @return A figure of the paired interactions
#' @export
TimePlot<-function(data,ligand,receptor,cell_from,cell_to,Time=NULL){
  if(is.null(Time)){
    Time=unique(data$time)
  }
  data<-data %>% filter(time %in% Time) %>% select(ligand,receptor,time,cell_type)
  data_long <- gather(data, gene, value, c(ligand,receptor), factor_key=TRUE) %>%
    filter((cell_type==cell_from & gene==ligand) | (cell_type==cell_to & gene==receptor))
  data_long$time<-as.factor(data_long$time)
  g<-ggplot(data_long,aes(x=time,y=value,color=gene))+
    geom_point(position=position_dodge(0.75)) +
    stat_summary(fun.y=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,position=position_dodge(0.75)) +
    stat_summary(fun.ymin=function(x)(mean(x)-sd(x)), fun.ymax=function(x)(mean(x)+sd(x)),geom="errorbar", width=0.1,position=position_dodge(0.75)) +
    theme_minimal()+theme(axis.line=element_line(),plot.title=element_text(size=14,face='bold',hjust=0.5))+ylab('gene expression')+xlab('time')+ggtitle(paste0(ligand,'-',receptor))
  g
  return(g)
}
