####Version1.2.0###
###Author Wangningqi####
#' Conduct Network analysis
#' @description A convenient and fast network analysis function, with output results suitable for cytoscape and gephi
#' @param input Input dataframe with otu/gene/taxa in row and sample ID in column,at least 5 replicates(more than 8 replicates are recommened).
#' @param inputtype Input dataframe type
#'
#' 1:dataframe with first column of OTUID and last column of taxonomy
#'
#' 2:dataframe with first column of OTUID/taxonomy
#'
#' 3:dataframe of all numeric
#'
#' @param input2 A second input data frame with otu/gene/taxa in row and sample ID in column. Default:NULL
#' @param input2type The second input data frame type. Details the same as above. Default:NULL
#' @param n  Only keep otu/gene/taxa appearing in n sample size
#' @param threshold Threshold of correlation r value
#' @param method A character string indicating which correlation coefficient  is to be computed. One of "pearson" or "spearman"
#' @param display If display a preview plot of network based on igraph. F for the first attempt is recommended in case of too many vertices and edges.
#'
#' @details
#'
#' 1. We had optimized the correlation algorithm to achieve a faster running speed. It takes less than 2 minute to calculate dataframe correlation and p value which more than 400 samples and 10000 OTUs for computer with dual Core i5 processor.
#' However, too many vertices(>2000) or links(>10000) may slow the statistical process and visualization,so we recommend that in your first attempt,set `display` paramter as `F` to have a preview.
#' Then you can adjust your n/threshold/method paramter to generate a suitable visualization network
#'
#' 2. We display a preview plot so as to adjusting your network. Generally a global figure (like we show in examples) with less than 1000 vertices and 5000 edges/links
#' is recommended. Further more,we recommend you to output the statistics and adjacency table and use software like cytoscape or gephi for better visualization.
#'
#' 3.\code{\link{left_join}}  from \code{\link{dplyr}} is available to match otu/gene/taxa annotaion with output results for further analysis.
#'
#' @note
#' 1. Replicates should be at least 5,more than 8 is recommend. See details in \code{\link{rcorr}}.
#'
#' 2. In case of too many edges/links or not a global network plot, you can stop the process immediately to provent wasting too much time.
#'
#' 3. Package magrittr,tidyr,\code{\link{Hmisc}}(correlation analysis),\code{\link{fdrtools}}(p value adjustment),\code{\link{igraph}}(create network graph for statistics and visualization) is required in this function.
#' @author  Wang Ningqi <2434066068@qq.com>
#' @return  One list contains a statistics table of network vertices/nodes and an adjacency table. One preview plot of network in the plot interface and an igraph object(named `igraph1`) in global environment.
#' @export
#'
#' @examples
#' ###data prepration###
#' data(testotu)
#' rownames(testotu)<-testotu[,1]
#' inputotu<-testotu[,-c(1,ncol(testotu))]
#' head(inputotu)
#'
#' ###one input network analysis###
#' network_result<-network_analysis(inputotu,3,10,0.9,"spearman",T)
#'
#' network_stat<- as.data.frame(network_result[1])  ##vertices table
#' head(network_stat)
#' network_adjacency<-as.data.frame(network_result[2]) ##adjacency table
#' head(network_adjacency)
#' network_matrix<-as.data.frame(network_result[3])  ##compelete adjacency matrix
#' network_matrix[1:10,1:10]
#' network_stat(igraph1) ##In case you want to see statistics again or do other analysis based on igraph.
#'
#' ###two inputs network analysis###
#' inputotu1<- inputotu[1:456,]
#' inputotu2<- inputotu[524:975,]
#' network_result<-network_analysis(input=inputotu1,inputtype=3,input2=inputotu2,input2type=3,n=10,threshold=0.85,method="spearman",display=T)
#'
#'
#'####incorrect demonstration!!##
#'###WARNING:it may takes too long to creat the graph!!###
#'
#'network_result<-network_analysis(inputotu,3,3,0.8,"spearman",T)
#'
#'#Total edges/links: 10199
#'#Total vertices: 826
#'##too many edges and not a global network####
#'
network_analysis<-function(input,inputtype,n,threshold,method,display,input2,input2type){
  require(igraph);require(fdrtool);require(tidyr);require(Hmisc);require(magrittr)
  if(inputtype==1){input1<-input[,-c(1,ncol(input))];rownames(input1)<-input[,1]}else
  if(inputtype==2){input1<-input[,-1];rownames(input1)<-input[,1]}else
  if(inputtype==3){input1<-input}else{print("Please choose correct inputtype(1,2,3)")}
  zero_count=function(input){length(which(input==0)) %>% return()}
  zerocount=apply(input1,1,zero_count)
  input1=input1[which(zerocount<=(ncol(input1)-n)),]
  input1=input1[which(rowSums(input1)>0),]
  if(missing(input2)){
    corr=rcorr(as.matrix(t(input1)),type=method)
    cor.p=corr$P;cor.p[is.na(cor.p)]<- 0
    fdr=fdrtool(as.numeric(cor.p), statistic="pvalue",plot=F,verbose = F) ##Global fdr correlation##
    cor.r=corr$r
    cor.q=matrix(fdr$qval,ncol=ncol(cor.p),nrow=nrow(cor.p));cor.q[is.nan(cor.q)]<- 0
    cor.r[cor.q>0.05|abs(cor.r)<threshold] = 0 ##fliter via threshold##
    cor.r[is.na(cor.r)]<-0 ##NA变0##
    cor.r[which(cor.r>0.9999)]<-0 ###对角线的1变成0##
    cor.r1=cor.r;cor.r1[which(cor.r1!=0)]<-1 ##无方向矩阵用于igraph##
    cutoff=which(rowSums(cor.r1)== 0)
    if(length(cutoff)==0){cor.r=cor.r;cor.r1=cor.r1}else{
    cor.r=cor.r[-c(cutoff),-c(cutoff)];cor.r1=cor.r1[-c(cutoff),-c(cutoff)]}
    cor.r1[which(cor.r>0)]<-1;cor.r1[which(cor.r<0)]<- -1;
    adj_matrix=cor.r1
    cor.r1[lower.tri(cor.r1,diag = T)]<-NA}else{
#####input&input2####
  if(input2type==1){input3<-input2[,-c(1,ncol(input2))];rownames(input3)<-input2[,1]}else
  if(input2type==2){input3<-input2[,-1];rownames(input3)<-input2[,1]}else
  if(input2type==3){input3<-input2}else{print("Please choose correct inputtype(1,2,3)")}
  zero_count=function(input){length(which(input==0)) %>% return()}
  zerocount=apply(input3,1,zero_count)
  input3=input3[which(zerocount<=(ncol(input3)-n)),]
  input3=input3[which(rowSums(input3)>0),]
  corr=rcorr(x=as.matrix(t(input1)),y=as.matrix(t(input3)),type=method)
  cor.p=corr$P[1:nrow(input1),(nrow(input1)+1):(nrow(input1)+nrow(input3))];cor.p[is.na(cor.p)]<- 0
  fdr=fdrtool(as.numeric(cor.p), statistic="pvalue",plot=F,verbose = F) ##Global fdr correlation##
  cor.q=matrix(fdr$qval,ncol=ncol(cor.p),nrow=nrow(cor.p));cor.q[is.nan(cor.q)]<- 0
  cor.r=corr$r[1:nrow(input1),(nrow(input1)+1):(nrow(input1)+nrow(input3))]
  cor.r[cor.q>0.05|abs(cor.r)<threshold] = 0 ##fliter via threshold##
  cor.r[is.na(cor.r)]<-0 ##NA变0##
  cor.r[which(cor.r>0.9999)]<-0 ###1变成0##
  cor.r1=cor.r;cor.r1[which(cor.r1!=0)]<-1 ##无方向矩阵用于igraph##
  cutoff1=which(rowSums(cor.r1)== 0);cutoff2=which(colSums(cor.r1)== 0)
  cor.r=cor.r[-c(cutoff1),-c(cutoff2)];cor.r1=cor.r1[-c(cutoff1),-c(cutoff2)]
  cor.r1[which(cor.r>0)]<-1;cor.r1[which(cor.r<0)]<- -1
  adj_matrix=cor.r1
    }
  ###创造邻接表###
  source<-rownames(cor.r1)
  adjacency<-cor.r1%>%cbind(source,.) %>%as.data.frame %>% gather("target","value",-source) %>% subset(.,.$value!=0)
  adjacency$value=as.numeric(adjacency$value)
  igraph1<<-graph_from_data_frame(adjacency,directed = F)##全局输出##
  network_stat(igraph1)
  if(length(E(igraph1))>10000){cat("

Warning:too many edges/links!Better STOP the process")}
  if(display==T){
    plot(igraph1,main="Co-occurrence network",vertex.frame.color=NA,vertex.label=NA,edge.width=1,
         vertex.size=5,edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))}
  V(igraph1)$degree<-igraph::degree(igraph1)
  set.seed(999)
  V(igraph1)$modularity <- membership(cluster_fast_greedy(igraph1))%>%as.numeric()
  nodes_list <- data.frame(nodes_id = V(igraph1)$name, degree = V(igraph1)$degree, modularity = V(igraph1)$modularity)
  nodes_list<-nodes_list[order(nodes_list$nodes_id,decreasing = F),]
  if(missing(input2)){
    cor.tempr=cor.r[order(rownames(cor.r),decreasing = F),];cor.tempr=cor.tempr[,order(colnames(cor.tempr),decreasing = F)]
    rownames(nodes_list)=rownames(cor.tempr)
    zi_pi8 <- zi.pi(nodes_list, cor.tempr, degree = 'degree', modularity_class = 'modularity')}else
    {cor.tempr=matrix(rep(0,nrow(nodes_list)^2),nrow=nrow(nodes_list),ncol=nrow(nodes_list),byrow=T)
    cor.tempr[1:nrow(cor.r),1:ncol(cor.r)]=cor.r
    rownames(cor.tempr)=c(rownames(cor.r),colnames(cor.r))
    colnames(cor.tempr)=c(colnames(cor.r),rownames(cor.r))
    cor.tempr=cor.tempr[order(rownames(cor.tempr),decreasing = F),];cor.tempr=cor.tempr[,order(colnames(cor.tempr),decreasing = F)]
    rownames(nodes_list)=rownames(cor.tempr)
    zi_pi8 <- zi.pi(nodes_list, cor.tempr, degree = 'degree', modularity_class = 'modularity')
    }

  output=data.frame(nodes_id = V(igraph1)$name, node_degree = V(igraph1)$degree, node_betw=betweenness(igraph1),
                    node_evcent=evcent(igraph1,scale = F)$vector,Clustering_coefficient=transitivity(igraph1,type="local"),
                    No.module = V(igraph1)$modularity,Zi=zi_pi8$Zi,Pi=zi_pi8$Pi)
  outlist=c(list(output),list(adjacency),list(adj_matrix)) %>%return()
}
