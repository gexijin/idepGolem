#' fct_analysis_random.R
#'
#'
#' @section fct_analysis_random.R functions:
#' \code{gene_group_heatmap} heatmap with color bar define gene groups
#'
#' 
#' @name fct_analysis_random.R
NULL

##### work in process, want to rewrite other code that is relays on first
gene_group_heatmap <- function (x,bar=NULL,n=-1,mycolor=1,clusterNames=NULL,sideColors=NULL ) {
	# number of genes to show
	ngenes = as.character( table(bar))
	if(length(bar) >n && n != -1) {ix = sort( sample(1:length(bar),n) ); bar = bar[ix]; x = x[ix,]  }
	if(! is.null(bar) )
		if(is.null(sideColors) ) 
			sideColors = mycolors

	# this will cutoff very large values, which could skew the color 
	x=as.matrix(x)-apply(x,1,mean)
	cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
	x[x>cutoff] <- cutoff
	cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
	x[x< cutoff] <- cutoff
	#colnames(x)= detectGroups(colnames(x))
	if(is.null(bar)) # no side colors
		heatmap.2(x,  Rowv =F,Colv=F, dendrogram ="none",
			col=heatColors[as.integer(mycolor),], density.info="none", trace="none", scale="none", keysize=.3
			,key=F, labRow = F
			#,RowSideColors = mycolors[bar]
			,margins = c(8, 24)
			,srtCol=45
		) else
		heatmap.2(x,  Rowv =F,Colv=F, dendrogram ="none",
			col=heatColors[as.integer(mycolor),], density.info="none", trace="none", scale="none", keysize=.3
			,key=F, labRow = F
			,RowSideColors = sideColors[bar]
			,margins = c(8, 24)
			,srtCol=45
		)
		
	if(!is.null(bar)) { 

		legend.text = paste("Cluster ", toupper(letters)[unique(bar)], " (N=", ngenes,")", sep="") 
		if( !is.null( clusterNames ) && length(clusterNames)>= length( unique(bar) ) )  
			legend.text = paste(clusterNames[ 1:length( unique(bar) )  ], " (N=", ngenes,")", sep="") 
		
		par(lend = 1)           # square line ends for the color legend
		legend("topright",      # location of the legend on the heatmap plot
		legend = legend.text, # category labels
		col = sideColors,  # color key
		lty= 1,             # line style
		lwd = 10 )           # line width
		}
}