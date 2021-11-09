#' fct_07_genome.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_07_genome.R functions:
#'
#'
#' @name fct_07_genome.R
NULL

#' MAIN PLOTLY FUNCTION
chromosome_plotly <- function(
  limma,
  select_contrast,
  all_gene_info,
  ignore_non_coding,
  limma_p_val_viz,
  limma_fc_viz,
  label_gene_symbol
) {
  # Default plot
  fake <- data.frame(a = 1:3, b = 1:3)
  p <- ggplot2::ggplot(fake, ggplot2::aes(x = a, y = b)) +
    ggplot2::geom_blank() +
    ggplot2::ggtitle("No genes with position info.") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )

  if(length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]  
  } else {
	top <- limma$top_genes
	ix <- match(input$select_contrast, names(top))
	if(is.na(ix)) {
      return(plotly::ggplotly(p))
    }
	top_1 <- top[[ix]] 
  }
  if(dim(top_1)[1] == 0) {
    return(plotly::ggplotly(p))
  }
  colnames(top_1) <- c("Fold", "FDR")
  
  x <- merge(
    top_1,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )
  
  colnames(x)[which(colnames(x) == "Row.names")] <- "ensembl_gene_id"

  # only coding genes? 
  if(ignore_non_coding) {
    x <- subset(x, gene_biotype == "protein_coding")
  }
  
  # If no chromosomes found. For example if user do not convert gene IDs.
  if(dim(x)[1] > 5) {
    x <- x[order(x$chromosome_name, x$start_position), ]
    x$ensembl_gene_id <- as.character(x$ensembl_gene_id)
    # If symbol is missing use Ensembl id
    x$symbol <- as.character(x$symbol)  
    ix <- which(is.na(x$symbol))
    ix2 <- which(nchar(as.character(x$symbol)) <= 2)
    ix3 <- which(duplicated(x$symbol))
    ix <- unique(c(ix, ix2, ix3))
    x$symbol[ix] <- x$ensembl_gene_id[ix] 
    
    x <- x[!is.na(x$chromosome_name), ]
    x <- x[!is.na(x$start_position), ]
            
    tem <- sort(table(x$chromosome_name), decreasing = T)
    # ch with less than 100 genes are excluded
    ch <- names(tem[tem >= 1])  
    if(length(ch) > 50) {
      # At most 50 ch
      ch <- ch[1:50]
      # ch. name less than 10 characters
      ch <- ch[nchar(ch) <= 12]
      ch <- ch[order(as.numeric(ch) ) ]
      tem <- ch
      # The numbers are continous from 1 to length(ch)
      ch <- 1:(length(ch))
      # The names are real chr. names
      names(ch) <- tem
      x <- x[which(x$chromosome_name %in% names(ch)), ]
      x <- droplevels(x)
      # Numeric encoding
      x$chNum <- 1 
      x$chNum <- ch[x$chromosome_name]
      
      # Use max position as chr. length before filtering
      ch_length_table <- aggregate(start_position ~ chromosome_name, data = x, max)
      # Add chr. numer 
      ch_length_table$chNum <- ch[ch_length_table$chromosome_name]
      ch_length_table <- ch_length_table[!is.na(ch_length_table$chNum), ]
      ch_length_table <- ch_length_table[order(ch_length_table$chNum), c(3, 2)]
      ch_length_table <- ch_length_table[order(ch_length_table$chNum), ]
      ch_length_table$start_position <- ch_length_table$start_position / 1e6

      # Only keep significant genes
      ix <- which(
        (x$FDR< as.numeric(limma_p_val_viz)) &
        (abs(x$Fold) > log2(as.numeric(limma_fc_viz)))
      )

      if(length(ix) > 5) { 
        # Remove non-significant / not selected genes
        # Keep a copy
        x0 <- x
        x <- x[ix, ]
        # Prepare coordinates
        # Mbp
        x$start_position <- x$start_position / 1000000
        # Distance between chs.
        chD <- 30
        # Max log2 fold 
        fold_cutoff <- 4
        
        # log2 fold within -5 to 5
        x$Fold[which(x$Fold > fold_cutoff)] <- fold_cutoff
        x$Fold[which(x$Fold < -1 * fold_cutoff)] <- -1 * fold_cutoff 
        x$Fold <- 4 * x$Fold
    
        x$y <- x$chNum * chD + x$Fold
        ch_total <- dim(ch_length_table)[1] 
        x$R <- as.factor(sign(x$Fold))
        
        colnames(x)[which(colnames(x) == "start_position")] <- "x"
        
        # Don't define x and y, so that we could plot use two datasets
        p <- ggplot2::ggplot() +
          ggplot2::geom_point(
            data = x,
            ggplot2::aes(x = x, y = y, colour = R, text = symbol),
            shape = 20,
            size = 0.2
          ) 

        if(label_gene_symbol) {
          p <- p + ggplot2::geom_text(
            data = x,
            ggplot2::aes(x = x, y = y, label = symbol),
            check_overlap = FALSE,
            angle = 45,
            size = 2,
            vjust = 0,
            nudge_y = 4 
          )
        }
        # Label y with ch names
        p <- p + ggplot2::scale_y_continuous(
          labels = paste(
            "chr",
            names(ch[ch_length_table$chNum]), sep = ""
          ), 
          breaks = chD * (1:ch_total), 
          limits = c(0, chD * (ch_total + 1) + 5)
        )
        
        # Draw horizontal lines for each ch.
        for(i in 1:dim(ch_length_table)[1]) {
          p <- p + ggplot2::annotate(
            "segment",
            x = 0,
            xend = ch_length_table$start_position[i],
            y = ch_length_table$chNum[i] * chD,
            yend = ch_length_table$chNum[i] * chD
          )
        }
             # change legend  http://ggplot2.tidyverse.org/reference/scale_manual.html
             p <- p + scale_colour_manual(name="",   # customize legend text
                                          values=c("red", "blue"),
                                          breaks=c("1","-1"),
                                          labels=c("Up", "Dn")) 
             p <- p + xlab("Position on chrs. (Mbp)") +  theme(axis.title.y=element_blank())      
             p <- p + theme(legend.position="none")

             incProgress(0.5)
             # add trend lines------------------------------------------
             x0 <- x0[x0$chromosome_name %in% unique(x$chromosome_name), ]
             x0$chNum <- 1 # numeric encoding
             x0$chNum <- ch[ x0$chromosome_name ]
             x0$start_position = x0$start_position/1e6 # Mbp             

             windowSize = as.numeric( input$MAwindowSize )#Mb            
             steps = as.numeric( input$MAwindowSteps ) # step size is then windowSize / steps       
             cutoff <- as.numeric(input$MAwindowCutoff) 

             x0$Fold <-  x0$Fold - mean(x0$Fold) # centering             

             incProgress(0.6)
                             
             for(i in 0:(steps-1)) {
               #step size is  windowSize/steps   
               # If windowSize=10 and steps = 2; then step size is 5Mb
               # 1.3 becomes 5, 11.2 -> 15 for step 1
               # 1.3 -> -5
               x0$x <- ( floor((x0$start_position - i * windowSize / steps)/ windowSize )  
                        + 0.5 + i / steps ) * windowSize

               movingAverage1 <- x0 %>%
                 select(chNum, x, Fold) %>%
                 filter( x >= 0) %>%   # beginning bin can be negative for first bin in the 2nd step
                 group_by(chNum, x) %>%
                 summarize( ma = mean(Fold),
                            n = n(),
                            pval = ifelse( n() >= 3 && sd(Fold) > 0, t.test(Fold)$p.value, 0 ) ) %>%
                 filter(!is.na(pval)) # na when only 1 data point?

               if(i == 0) {
                 movingAverage <- movingAverage1
               } else {
                 movingAverage <- rbind(movingAverage, movingAverage1)  
               }        
             }
 
            
             # translate fold to y coordinates
             movingAverage <- movingAverage %>%
                filter(n >= 3) %>%
#                mutate( pval = p.adjust(pval, method = "fdr", n = length(pval[pval < 0.9]) ) )
                mutate( pval = p.adjust(pval, method = "fdr" ) ) %>%
                filter( pval < as.numeric(input$chRegionPval) ) %>%
                mutate( y = ifelse(ma > 0, 1, -1)) %>% # upper bound
                mutate(y = chNum * chD + 3 * y) %>%
                mutate( ma = ifelse(ma > 0, 1, -1)) %>%
                mutate( ma = as.factor(ma))

              # significant regions are marked as horizontal error bars 
             if(dim(movingAverage)[1] > 0) {
               p <- p +
                 geom_errorbarh(data = movingAverage, aes(x = x, 
                                                          y = y, 
                                                          xmin = x -windowSize/2, 
                                                          xmax = x + windowSize/2,
                                                          colour = ma), 
                                 size = 2, 
                                 height = 15 )

                 # label significant regions
                 sigCh <- sort(table(movingAverage$chNum), decreasing = TRUE)
                 sigCh <- names(ch)[ as.numeric(names(sigCh)) ]
                 if(length(sigCh) <= 5) { # more than 5 just show 5
                   sigCh <- paste0("chr", sigCh, collapse = ", ")
                 } else {
                   sigCh <- sigCh[1:5]
                   sigCh <- paste0("chr", sigCh, collapse = ", ")                  
                   sigCh <- paste0(sigCh,", ...")
                 }

                 sigCh <- paste(dim(movingAverage)[1], 
                                " enriched regions \n(",
                                round( sum(chLengthTable$start_position)/ windowSize * steps * as.numeric(input$chRegionPval), 2),
                                          " expected)  detected on:\n ", sigCh)
                 
               p <- p + annotate(geom = "text", 
                          x = max(x$x) * 0.70,
                          y = max(x$y) * 0.90,
                          label = sigCh)


             }

         } # have genes after filter
			
      }  # have 5+ genes to begin with
              incProgress(1)
	  ggplotly(p)
    
}