#' mod_02_pre_process.R
#'
#' 
#' @section mod_02_pre_process.R functions:
#'
#'
#' @name mod_02_pre_process.R
NULL

##### work in process, want to rewrite other code that is relays on first
# plot the expression of one or more genes in the preprocess tab
plotGenes <- function(converted_data, all_gene_info, readSampleInfo, geneSearch, genePlotBox, useSD, selectOrg){
   symbols <- rownames(converted_data)
   if (selectOrg != "NEW" &&  ncol(all_gene_info) != 1) {
     ix <- match(symbols, all_gene_info[, 1])
     # symbol really exists? 
     if (sum(is.na(all_gene_info$symbol)) != nrow(all_gene_info)) {  
       Symbols = as.character(allGeneInfo$symbol[ix] )
       Symbols[which( nchar(Symbols) <= 2 ) ] <- rownames(x) [which( nchar(Symbols) <= 2 ) ] 
     }
   }
   x = as.data.frame(x)
   x$Genes = Symbols
   
   # matching from the beginning of symbol
   searchWord = gsub("^ ", "", geneSearch )
   ix = which(regexpr(  paste("^" , toupper(searchWord),sep="")   ,toupper(x$Genes)) > 0)
   if(grepl(" $", searchWord)  )  # if there is space character at the end, do exact match
     ix = match(gsub(" ","", toupper(searchWord)), toupper(x$Genes) )
   
   if(grepl(",|;", searchWord)  ) { # if there is comma or semicolon, split into multiple words
     Words <- unlist( strsplit(searchWord,",|;") ) # split words
     Words <- gsub(" ", "", Words)
     ix = match( toupper(Words), toupper(x$Genes) )
   }
   ix = ix[!is.na(ix)] # remove NAs
   # too few or too many genes found
   if(length(ix) == 0 | length(ix) > 50 ) return(NULL)
   # no genes found
   
   mdf = melt(x[ix,],id.vars="Genes", value.name="value", variable.name="samples")
   # bar plot of individual samples
   p1 <- ggplot(data=mdf, aes(x=samples, y=value, group = Genes, shape=Genes, colour = Genes)) +
     geom_line() +
     geom_point( size=5,  fill="white")+ #shape=21  circle
     #theme(axis.text.x = element_text(size=16,angle = 45, hjust = 1)) +
     labs(y="Transformed expression level") +
     coord_cartesian(ylim = c(0, max(mdf$value)))
   p1 <- p1 + theme(plot.title = element_text(size = 16,hjust = 0.5)) + # theme(aspect.ratio=1) +
     theme(axis.text.x = element_text(angle=45, size = 16, hjust=1),
           axis.text.y = element_text( size = 16),
           axis.title.x = element_blank(),
           axis.title.y = element_text( size = 16) ) +
     theme(legend.text=element_text(size=12))	
   
   
   #ggplotly(p) %>% layout(margin = list(b = 250,l=100))  # prevent cutoff of sample names
   
   # Barplot with error bars
   mdf$count = 1
   g = detectGroups(mdf$samples, readSampleInfo)
   mdf$g = g	
   
   options(dplyr.summarise.inform = FALSE)
   #calculate mean, SD, N, per gene per condition
   summarized <- mdf %>% 
     group_by(g, Genes) %>%  
     summarise(Mean = mean(value), SD = sd(value), N = sum(count))
   colnames(summarized)= c("Samples","Genes","Mean","SD","N")
   summarized$SE = summarized$SD / sqrt(summarized$N)	
   
   if(grepl(",|;", searchWord)  ) { # re-order according to user input, not alphabetically
     levels <- unique(summarized$Genes)
     iy <- match(toupper(Words), toupper(levels) )
     levels <- levels[iy]
     summarized$Genes <- factor( summarized$Genes, levels = levels) 
   }
   
   #http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization
   p2 <- ggplot(summarized, aes(x=Genes, y=Mean,fill=Samples) ) + # data & aesthetic mapping
     geom_bar(stat="identity", position=position_dodge()) + # bars represent average
     geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.2,position=position_dodge(.9)) +
     labs(y="Expression Level") 
   if(useSD == 1) { 
     p2 <- ggplot(summarized, aes(x=Genes, y=Mean,fill=Samples) ) + # data & aesthetic mapping
       geom_bar(stat="identity", position=position_dodge()) + # bars represent average
       geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2,position=position_dodge(.9)) +
       labs(y="Expression Level") 
   }
   
   p2 <- p2 +  theme(plot.title = element_text(size = 16,hjust = 0.5)) + # theme(aspect.ratio=1) +
     theme(axis.text.x = element_text(angle=45, size = 16, hjust=1),
           axis.text.y = element_text( size = 16),
           axis.title.x = element_blank(),
           axis.title.y = element_text( size = 16) ) +
     theme(legend.text=element_text(size=16))
   
   if( genePlotBox == 1)  
     return(p1) else 
     return(p2)
   
 }