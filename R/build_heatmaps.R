build_heatmaps <-
    function(MEM_values,
             cluster.MEM = "both",
             cluster.medians = "none",
             cluster.IQRs = "none",
             display.thresh = 1,
             output.files = FALSE,
             labels = FALSE,
             only.MEMheatmap = FALSE) {
        #heatmap clustering variables
        dendro_var_MEM <- cluster.MEM
        dendro_var_med <- cluster.medians
        dendro_var_IQR <- cluster.IQRs
        
        ##Build MEM heatmap with natural language output
        #Get MEM values and the max and min MEM value
        heatmap_data <- (MEM_values[[5]])[[1]]
        output.dir <- "./output files/"
        if (length(which(apply(heatmap_data, 1, function(row) any(row < 0)))) == 0){ #referenceless MEM (no negative values)
            scale_max <- 10
            #Initialize heatmap color palette
            heat_palette_MEM <- colorRampPalette(c("black", "yellow", "#FAF7C9"))
            #splitting points for binning data into colors
            pairs.breaks_MEM <-
                c(seq(0, scale_max / 6.6, by = 0.1),
                  seq(scale_max / 6.6, scale_max / 3.3, by = 0.1),
                  seq(scale_max / 3.3, scale_max, by = 0.1)
                )
            pairs.breaks_MEM <- unique(pairs.breaks_MEM)
        }else{ #regular MEM (positive and negative MEM scores)
            scale_max <- 10
            scale_min <- (-10)
            #Initialize heatmap color palette
            heat_palette_MEM <- colorRampPalette(c("#A8C3F4", "#17499B", "black", "yellow", "#FAF7C9"))
            #splitting points for binning data into colors
            pairs.breaks_MEM <-
                c(
                    seq(scale_min, scale_min / 3.3, by = 0.2),
                    seq(scale_min / 3.3, scale_min / 6.6, by = 0.2),
                    seq(scale_min / 6.6, 0, by = 0.2),
                    seq(0, scale_max / 6.6, by = 0.2),
                    seq(scale_max / 6.6, scale_max / 3.3, by = 0.2),
                    seq(scale_max / 3.3, scale_max, by = 0.2)
                )
            pairs.breaks_MEM <- unique(pairs.breaks_MEM)}

        #Initialize heatmap for medians
        medians_MEM_values <- (MEM_values[[1]])[[1]]
        scale_max_med <- max(medians_MEM_values)
        heat_palette_med <- colorRampPalette(c("black", "yellow", "#FAF7C9"))
        pairs.breaks_med <- c(seq(0, scale_max_med / 6.6, by = 0.1),
                              seq(scale_max_med / 6.6, scale_max_med / 3.3, by = 0.1),
                              seq(scale_max_med / 3.3, scale_max_med, by = 0.1))
        pairs.breaks_med <- unique(pairs.breaks_med)

        #Initialize heatmap for IQR
        IQR_MEM_values <- (MEM_values[[3]])[[1]]
        if (max(IQR_MEM_values) <= 2.8) {
          scale_max_IQR <- 2.81
        }else{
          scale_max_IQR <- max(IQR_MEM_values)
        }
        heat_palette_IQR <- colorRampPalette(c("black", "yellow", "#FAF7C9"))
        pairs.breaks_IQR <- c(seq(0.5, (scale_max_IQR + 0.5) / 6.6, by = 0.1),
                              seq((scale_max_IQR + 0.5) / 6.6, (scale_max_IQR + 0.5) / 3.3, by = 0.1),
                              seq((scale_max_IQR + 0.5) / 3.3, scale_max_IQR, by = 0.1))
        pairs.breaks_IQR <- unique(pairs.breaks_IQR)

        #Round MEM enrichment vals
        MEM_vals_scale <- as.matrix(round(heatmap_data, 0))
        
        #Initialize variables
        #Call create.labels
        new_rownames <- create.labels(MEM_vals_scale, display.thresh, heatmap_data)
        new_rownames_txt <- create.labels.txt(MEM_vals_scale, display.thresh, heatmap_data)
        
        #Specify heatmap clustering parameters
        if (dendro_var_MEM == "both") { #cluster MEM rows and columns
            Colv_var_MEM = TRUE
            Rowv_var_MEM = TRUE
        }
        if (dendro_var_MEM == "row") { #only cluster MEM rows
            Colv_var_MEM = FALSE
            Rowv_var_MEM = TRUE
        }
        if (dendro_var_MEM == "col") { #only cluster MEM columns
            Colv_var_MEM = TRUE
            Rowv_var_MEM = FALSE
        }
        if (dendro_var_MEM == "none") { #no clustering
            Colv_var_MEM = FALSE
            Rowv_var_MEM = FALSE
        }
        if (dendro_var_med == "both") { #cluster median rows and columns
            Colv_var_med = TRUE
            Rowv_var_med = TRUE
        }
        if (dendro_var_med == "row") { #only cluster median rows
            Colv_var_med = FALSE
            Rowv_var_med = TRUE
        }
        if (dendro_var_med == "col") { #only cluster median columns
            Colv_var_med = TRUE
            Rowv_var_med = FALSE
        }
        if (dendro_var_med == "none") { #no clustering
            Colv_var_med = FALSE
            Rowv_var_med = FALSE
        }
        if (dendro_var_IQR == "both") { #cluster IQR rows and columns
            Colv_var_IQR = TRUE
            Rowv_var_IQR = TRUE
        }
        if (dendro_var_IQR == "row") { #only cluster IQR rows
            Colv_var_IQR = FALSE
            Rowv_var_IQR = TRUE
        }
        if (dendro_var_IQR == "col") { #only cluster IQR columns
            Colv_var_IQR = TRUE
            Rowv_var_IQR = FALSE
        }
        if (dendro_var_IQR == "none") { #no clustering
            Colv_var_IQR = FALSE
            Rowv_var_IQR = FALSE
        }


        # Print MEM heatmap according to cluster spec
        if (length(which(apply(heatmap_data, 1, function(row) any(row < 0)))) == 0){ #referenceless MEM (no negative values)
            title_MEM <- "   MEM* Heatmap" 
        }else{ #regular referenced MEM (positive and negative values)
            title_MEM <- "   MEM Heatmap"}

        #MEM heatmap arguments
        args_heatmap_MEM <- list(
            heatmap_data,
            main = title_MEM,
            dendrogram = dendro_var_MEM,
            Rowv = Rowv_var_MEM,
            Colv = Colv_var_MEM,
            breaks = pairs.breaks_MEM,
            revC = FALSE,
            symm = FALSE,
            symkey = FALSE,
            symbreaks = FALSE,
            scale = "none",
            cexRow = 0.8,
            cexCol = 0.8,
            col = heat_palette_MEM,
            key = TRUE,
            key.title = NA, key.xlab = NA, key.ylab = NA,
            lhei = c(0.5, 1.5),
            lwid = c(0.5, 2.2),
            trace = "none"
        )
        
        #add labels to MEM arguments if specified
        if (labels == TRUE) {
            args_heatmap_MEM <- c(
                args_heatmap_MEM,
                list(
                    labRow = new_rownames,
                    margins = c(5, 15)
                )
            )
        }else{
            args_heatmap_MEM <- c(
                args_heatmap_MEM,
                list(
                    margins = c(5, 10)
                )
            )
        }
        
        
        #create the MEM heatmap
        if (nrow(MEM_values[[1]][[1]])==1){ #check for one cluster
          markers=colnames(heatmap_data) #get marker names
          cluster=rownames(heatmap_data) #get cluster name
          values=heatmap_data[1,] #get MEM scores
          df <- data.frame(x=markers, y=cluster, z=values)
          
          print(ggplot(df, aes(markers, cluster)) + geom_tile(aes(fill= values)) + scale_fill_gradientn(colors = heat_palette_MEM(100), limits=c(0, 10)) + labs(title=title_MEM, x=element_blank(), y=element_blank()) + theme_classic() + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.3), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank()))
          
          #set some variables for later
          table <- list(rowInd=NA, colInd=NA)
          table$rowInd <- 1
          table$colInd <- c(1:ncol(MEM_values[[1]][[1]]))
        }else{
          table <-do.call(heatmap.2, args_heatmap_MEM)
        }

        clustered_matrix <- heatmap_data[rev(table$rowInd), table$colInd]
        matrix.test<- as.matrix(new_rownames)
        matrix.test_txt <- as.matrix(new_rownames_txt)
        
        #generate the actual MEM labels
        if (length(which(apply(heatmap_data, 1, function(row) any(row < 0)))) >0){
            enrichment_score_ordered_txt <- matrix.test_txt[rev(table$rowInd), ]
        }else{
            enrichment_score_ordered_txt <- matrix.test_txt[rev(table$rowInd), ]}

        if (length(which(apply(heatmap_data, 1, function(row) any(row < 0)))) >0){
            enrichment_score_ordered <- matrix.test[rev(table$rowInd), ]
        }else{
            enrichment_score_ordered <- matrix.test[rev(table$rowInd), ]}


        # Print medians heatmap according to cluster spec
        reorder_medians <- as.matrix(MEM_values[[1]][[1]])[rev(table$rowInd), table$colInd]
        # Print IQR heatmap according to cluster spec
        reorder_IQR <- as.matrix(MEM_values[[3]][[1]])[rev(table$rowInd), table$colInd]
        #transpose if only one cluster
        if (nrow(MEM_values[[1]][[1]])==1){
          reorder_medians <- t(reorder_medians)
          rownames(reorder_medians) <- rownames(MEM_values[[1]][[1]])
          
          reorder_IQR <- t(reorder_IQR)
          rownames(reorder_IQR) <- rownames(MEM_values[[1]][[1]])
        }

        if (only.MEMheatmap == FALSE) {
            #median heatmap arguments
            args_heatmap_MEDs <- list(
                main = "Median Heatmap",
                breaks = pairs.breaks_med,
                revC = FALSE,
                symm = FALSE,
                symkey = FALSE,
                symbreaks = FALSE,
                scale = "none",
                cexRow = 0.8,  
                cexCol = 0.8,
                col = heat_palette_med,
                key = TRUE,
                key.title = NA, key.xlab = NA, key.ylab = NA,
                lhei = c(0.5, 1.5),
                lwid = c(0.5, 2.2),
                trace = "none"
            )

            if (cluster.medians != "none") { #cluster the median heatmap
                args_heatmap_MEDs <- c(
                    list(as.matrix(MEM_values[[1]][[1]])),
                    args_heatmap_MEDs,
                    list(dendrogram = dendro_var_med,
                         Rowv = Rowv_var_med,
                         Colv = Colv_var_med,
                         margins = c(5, 10))
                )
                #create heatmap
                if (nrow(MEM_values[[1]][[1]])==1){ #check if one cluster
                  # image(z=t(MEM_values[[1]][[1]]), main = "Median Heatmap", col = pal) #console.output
                  values=MEM_values[[1]][[1]]
                  df <- data.frame(x=markers, y=cluster, z=values)
                  print(ggplot(df, aes(markers, cluster)) + geom_tile(aes(fill= values)) + scale_fill_gradientn(colors = heat_palette_MEM(100), limits=c(0, scale_max_med)) + labs(title="Median Heatmap", x=element_blank(), y=element_blank()) + theme_classic() + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.3), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank()))
                  medians_table <- list(rowInd=NA, colInd=NA)
                  medians_table$rowInd <- 1
                  medians_table$colInd <- c(1:ncol(MEM_values[[1]][[1]]))
                }else{
                  medians_table <- do.call(heatmap.2, args_heatmap_MEDs)
                  reorder_medians <- as.matrix(MEM_values[[1]][[1]])[rev(medians_table$rowInd), medians_table$colInd]
                }
            }else{ #don't cluster the median heatmap
                args_heatmap_MEDs <- c(
                    list(reorder_medians),
                    args_heatmap_MEDs,
                    list(dendrogram = "none",
                         Rowv = FALSE,
                         Colv = FALSE,
                         margins = c(5, 15))
                )
                #create heatmap
                if (nrow(MEM_values[[1]][[1]])==1){ #check if one cluster
                  # image(z=t(MEM_values[[1]][[1]]), main = "Median Heatmap", col = pal) #console output
                  values=MEM_values[[1]][[1]]
                  df <- data.frame(x=markers, y=cluster, z=values)
                  print(ggplot(df, aes(markers, cluster)) + geom_tile(aes(fill= values)) + scale_fill_gradientn(colors = heat_palette_MEM(100), limits=c(0, scale_max_med)) + labs(title="Median Heatmap", x=element_blank(), y=element_blank()) + theme_classic() + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.3), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank()))
                  table2 <- list(rowInd=NA, colInd=NA)
                  table2$rowInd <- 1
                  table2$colInd <- rev(c(1:ncol(MEM_values[[1]][[1]])))
                }else{
                  table2 <- do.call(heatmap.2, args_heatmap_MEDs)
                }
            }

            #IQR heatmap arguments
            args_heatmap_IQRs <- list(
                main = "IQR Heatmap",
                breaks = pairs.breaks_IQR,
                revC = FALSE,
                symm = FALSE,
                symkey = FALSE,
                symbreaks = FALSE,
                scale = "none",
                cexRow = 0.8,  
                cexCol = 0.8,
                col = heat_palette_IQR,
                key = TRUE,
                key.title = NA, key.xlab = NA, key.ylab = NA,
                lhei = c(0.5, 1.5),
                lwid = c(0.5, 2.2),
                trace = "none"
            )

            if (cluster.IQRs != "none") { #cluster the IQR heatmap
                args_heatmap_IQRs <- c(
                    list(as.matrix(MEM_values[[3]][[1]])),
                    args_heatmap_IQRs,
                    list(dendrogram = dendro_var_IQR,
                         Rowv = Rowv_var_IQR,
                         Colv = Colv_var_IQR,
                         margins = c(5, 10))
                )
                #create IQR heatmap
                if (nrow(MEM_values[[1]][[1]])==1){ #check if one cluster
                  # image(z=t(MEM_values[[3]][[1]]), main = "IQR Heatmap", col = pal) #console output
                  values=MEM_values[[3]][[1]]
                  df <- data.frame(x=markers, y=cluster, z=values)
                  print(ggplot(df, aes(markers, cluster)) + geom_tile(aes(fill= values)) + scale_fill_gradientn(colors = heat_palette_MEM(100), limits=c(0, scale_max_IQR)) + labs(title="IQR Heatmap", x=element_blank(), y=element_blank()) + theme_classic() + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.3), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank()))
                  IQR_table <- list(rowInd=NA, colInd=NA)
                  IQR_table$rowInd <- 1
                  IQR_table$colInd <- c(1:ncol(MEM_values[[1]][[1]]))
                  reorder_IQR <- as.matrix(t(MEM_values[[3]][[1]])[rev(IQR_table$rowInd), IQR_table$colInd])
                }else{
                  IQR_table <- do.call(heatmap.2, args_heatmap_IQRs)
                  reorder_IQR <- as.matrix(MEM_values[[3]][[1]])[rev(IQR_table$rowInd), IQR_table$colInd]
                }
            }else{ #don't cluster the IQR heatmap
                args_heatmap_IQRs <- c(
                    list(reorder_IQR),
                    args_heatmap_IQRs,
                    list(dendrogram = "none",
                         Rowv = FALSE,
                         Colv = FALSE,
                         margins = c(5, 15))
                )
                #create IQR heatmap
                if (nrow(MEM_values[[1]][[1]])==1){ #check if one cluster
                  # image(z=t(MEM_values[[3]][[1]]), main = "IQR Heatmap", col = pal) #console output
                  values=MEM_values[[3]][[1]]
                  df <- data.frame(x=markers, y=cluster, z=values)
                  print(ggplot(df, aes(markers, cluster)) + geom_tile(aes(fill= values)) + scale_fill_gradientn(colors = heat_palette_MEM(100), limits=c(0, scale_max_IQR)) + labs(title="IQR Heatmap", x=element_blank(), y=element_blank()) + theme_classic() + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.3), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank()))
                  table3 <- list(rowInd=NA, colInd=NA)
                  table3$rowInd <- 1
                  table3$colInd <- c(1:ncol(MEM_values[[1]][[1]]))
                }else{
                  table3 <- do.call(heatmap.2, args_heatmap_IQRs)
                }
            }
        }

        dir.create(file.path(getwd(), output.dir), showWarnings = FALSE)

        if (output.files == TRUE) {

            # file attributes
            file_prefix <- file.path(output.dir, strftime(Sys.time(), "%Y-%m-%d_%H%M%S"))

            # ----- pdf output ----------------------------
            cairo_pdf(paste0(file_prefix, " MEM heatmap.pdf"), width = 15, onefile = TRUE)

            # MEM heatmap plot
            if (nrow(MEM_values[[1]][[1]])==1){ #check if one cluster
              # image(z=t(heatmap_data), main = title_MEM, col = pal) 
              values= heatmap_data[1,]
              df <- data.frame(x=markers, y=cluster, z=values)
              print(ggplot(df, aes(markers, cluster)) + geom_tile(aes(fill= values)) + scale_fill_gradientn(colors = heat_palette_MEM(100), limits=c(0, 10)) + labs(title=title_MEM, x=element_blank(), y=element_blank()) + theme_classic() + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.3), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank()))
            }else{
              do.call(heatmap.2, args_heatmap_MEM)
            }

            if (only.MEMheatmap == FALSE) {
              if (nrow(MEM_values[[1]][[1]])==1){ #check if one cluster
                # Plot median heatmap
                # image(z=t(MEM_values[[1]][[1]]), main = "Median Heatmap", col = pal)
                values=MEM_values[[1]][[1]] 
                df <- data.frame(x=markers, y=cluster, z=values)
                print(ggplot(df, aes(markers, cluster)) + geom_tile(aes(fill= values)) + scale_fill_gradientn(colors = heat_palette_MEM(100), limits=c(0, scale_max_med)) + labs(title="Median Heatmap", x=element_blank(), y=element_blank()) + theme_classic() + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.3), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank()))
                # Plot IQR heatmap
                # image(z=t(MEM_values[[3]][[1]]), main = "IQR Heatmap", col = pal) 
                values=MEM_values[[3]][[1]] 
                df <- data.frame(x=markers, y=cluster, z=values)
                print(ggplot(df, aes(markers, cluster)) + geom_tile(aes(fill= values)) + scale_fill_gradientn(colors = heat_palette_MEM(100), limits=c(0, scale_max_IQR)) + labs(title="IQR Heatmap", x=element_blank(), y=element_blank()) + theme_classic() + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.3), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank()))
              }else{
                # Plot median heatmap
                do.call(heatmap.2, args_heatmap_MEDs)
                # Plot IQR heatmap
                do.call(heatmap.2, args_heatmap_IQRs)
              }
            }

            dev.off()
            # ----- end of pdf output ---------------------

            # Write data to text files
            #MEM matrix
            write.table(
                as.matrix(clustered_matrix),
                paste0(file_prefix, " MEM matrix.txt"),
                sep = "\t",
                row.names = TRUE)
            #median matrix
            write.table(
                as.matrix(reorder_medians),
                paste0(file_prefix, " Medians matrix.txt"),
                sep = "\t",
                row.names = TRUE)
            #IQR matrix
            write.table(
                as.matrix(reorder_IQR),
                paste0(file_prefix, " IQRs matrix.txt"),
                sep = "\t",
                row.names = TRUE)

            #MEM label txt file
            if (((MEM_values[[6]])[[1]]) == 0){ #check for multiple files
                write.table(
                    enrichment_score_ordered,
                    paste0(file_prefix, " enrichment score-rownames.txt"),
                    sep = "\t")
            }else{ #multiple files for the clusters
                filenames <- unlist(MEM_values[[6]])
                matrix.filenames <- as.matrix(filenames)
                filenames_ordered <- matrix.filenames[rev(table$rowInd), ]
                new_rownames_filenames <-
                    cbind(filenames_ordered, enrichment_score_ordered)
                colnames(new_rownames_filenames) <- c("File", "MEM label")
                write.table(
                    new_rownames_filenames,
                    paste0(file_prefix, " enrichment score-rownames.txt"),
                    sep = "\t"
                )
            }
        }
        if (((MEM_values[[6]])[[1]]) == 0){
            show(paste("MEM label for cluster",enrichment_score_ordered))
        }else{
            filenames <- unlist(MEM_values[[6]])
            matrix.filenames <- as.matrix(filenames)
            filenames_ordered <- matrix.filenames[rev(table$rowInd), ]
            show(filenames_ordered)
            show(paste("MEM label ",enrichment_score_ordered))
        }
        
        #return the MEM labels
        return(enrichment_score_ordered)
    }
