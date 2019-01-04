#' Network Viewing of cell-cell communication
#'
#' This function loads the significant interactions as a dataframe, and colors
#' represent different types of cells as a structure. The width of edges represent
#' the strength of the communication. Labels on the edges show exactly how many
#' interactions exist between two types of cells.
#'
#' @references Csardi G, Nepusz T: The igraph software package for complex network
#' research, InterJournal, Complex Systems 1695. 2006.
#' http://igraph.org
#' @param data A dataframe containing ligand-receptor pairs and corresponding
#' cell typesused to do the plotting
#' @param col Colors used to represent different cell types
#' @param label Whether or not shows the label of edges (number of connections
#' between different cell types)
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' This can be a logical or a numeric vector or scalar.
#' First the vector is replicated to have the same length as the number of
#' edges in the graph. Then it is interpreted for each edge separately.
#' A numeric value specifies the curvature of the edge; zero curvature means
#' straight edges, negative values means the edge bends clockwise, positive
#' values the opposite. TRUE means curvature 0.5, FALSE means curvature zero
#' @param shape The shape of the vertex, currently “circle”, “square”,
#' “csquare”, “rectangle”, “crectangle”, “vrectangle”, “pie” (see
#' vertex.shape.pie), ‘sphere’, and “none” are supported, and only by the
#' plot.igraph command. “none” does not draw the vertices at all, although
#' vertex label are plotted (if given). See shapes for details about vertex
#' shapes and vertex.shape.pie for using pie charts as vertices.
#' @param layout The layout specification. It must be a call to a layout
#' specification function.
#' @param vertex.size The size of vertex
#' @param margin The amount of empty space below, over, at the left and right
#'  of the plot, it is a numeric vector of length four. Usually values between
#'  0 and 0.5 are meaningful, but negative values are also possible, that will
#'  make the plot zoom in to a part of the graph. If it is shorter than four
#'  then it is recycled.
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param arrow.width The width of arrows
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.max.width The maximum arrow size
#' @import network
#' @import igraph
#' @return A network graph of the significant interactions
#' @export
NetView<-function(data,col,label=TRUE,edge.curved=0.5,shape='circle',layout=nicely(),vertex.size=20,margin=0.2,
                  vertex.label.cex=1.5,vertex.label.color='black',arrow.width=1.5,edge.label.color='black',edge.label.cex=1,edge.max.width=10){
  net<-data %>% group_by(cell_from,cell_to) %>% dplyr::summarize(n=n())
  net<-as.data.frame(net,stringsAsFactors=FALSE)
  g<-graph.data.frame(net,directed=TRUE)
  edge.start <- ends(g, es=E(g), names=FALSE)
  coords<-layout_(g,layout)
  coords_scale<-scale(coords)
  loop.angle<-ifelse(coords_scale[V(g),1]>0,-atan(coords_scale[V(g),2]/coords_scale[V(g),1]),pi-atan(coords_scale[V(g),2]/coords_scale[V(g),1]))
  V(g)$size<-vertex.size
  V(g)$color<-col[V(g)]
  V(g)$label.color<-vertex.label.color
  V(g)$label.cex<-vertex.label.cex
  if(label){
    E(g)$label<-E(g)$n
  }
  E(g)$width<-1+edge.max.width/(max(E(g)$n)-min(E(g)$n))*(E(g)$n-min(E(g)$n))
  E(g)$arrow.width<-arrow.width
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
  E(g)$color<-V(g)$color[edge.start[,1]]
  E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin)
  return(g)
}
