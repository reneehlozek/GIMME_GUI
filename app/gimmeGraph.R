########################################################################
########################################################################
# This function is written to provide individual-level graphics 
# representing contemporaneous and lagged relationships. ###############
########################################################################
########################################################################
# Author: Stephanie Lane; slane@unc.edu ################################
# Date of creation: May 24, 2016 #######################################
########################################################################
########################################################################
# This code depends on the readxl package, the R.matlab package, and the 
# qgraph package. ######################################################
########################################################################
########################################################################
gimmeGraph <- function(input, 
                       filetype, 
                       sep, 
                       header, 
                       weighted){
  # create directory containing plots as subdirectory of input directory
  dir.create(file.path(input, "plots"), showWarnings = FALSE)
  # list all input matrices
  all_mats  <- list.files(input, full.names = TRUE)
  # make sure to not include any subdirectories
  all_mats  <- all_mats[tools::file_ext(all_mats) != ""]
  # get path of matrix containing group-level connections
  group_mat <- all_mats[grep("group.", all_mats)]
  # set all_ind to be all non-group matrices
  all_ind   <- all_mats[-grep("group.", all_mats)]
  # get base filenames
  filenames <- tools::file_path_sans_ext(basename(all_ind))  
  
  if (weighted == TRUE){
    filenames <- paste0(filenames, "_weighted")
  } else {
    filenames <- paste0(filenames, "_unweighted")
  }
  # read in xlsx and xls files
  if (filetype == "xlsx" | filetype == "xls"){
    group <- readxl::read_excel(group_mat, col_names = header)
    n_col <- ncol(group)
    n_row <- nrow(group)
  }
  # read in txt files
  if (filetype == "txt"){
      group <- read.table(group_mat, header = header, sep = sep)
    n_col <- ncol(group)
    n_row <- nrow(group)
  }
  if (filetype == "csv"){
      group <- read.table(group_mat, header = header, sep = ",")
    n_col <- ncol(group)
    n_row <- nrow(group)
  }
  if (filetype == "mat"){
    group <- R.matlab::readMat(group_mat)[[1]]
    n_col <- ncol(group)
    n_row <- nrow(group)
  }
  
  rois    <- n_col/2

  for (i in 1:length(all_ind)){
    # create matrix for colors
    col_mat <- matrix(NA, ncol = n_col, nrow = n_row)
    # read in individual matrix
    if (filetype == "xlsx" | filetype == "xls"){
      ind_mat <- readxl::read_excel(all_ind[i], col_names = header)
    }
    # read in txt files
    if (filetype == "txt"){
        ind_mat <- read.table(all_ind[i], header = header, sep = sep)
    }
    if (filetype == "csv"){
        ind_mat <- read.table(all_ind[i], header = header, sep = ",")
    }
    if (filetype == "mat"){
      ind_mat <- R.matlab::readMat(all_ind[i])[[1]]
    }
    
    # set element of individual-level matrix to .01 
    # if estimate is rounded to zero but a group-level path 
    # does exist. otherwise, colors/values are misaligned
    ind_mat[ind_mat == 0 & group == 1] <- .01
    # set elements of color matrix to grey where 
    # there exists a nonzero element in the individual matrix
    col_mat[ind_mat > 0] <- "red1"
    col_mat[ind_mat < 0] <- "mediumblue"
    # set elements of color matrix to black where
    # there exists a 1 in the group matrix
    col_mat[group == 1]   <- "black"
    # if weighted plots are not desired, overwrite beta value
    # with thick line for group and thin line for individual 
    if (weighted == FALSE) {
      ind_mat[col_mat == "black"]   <- 2
      ind_mat[col_mat == "grey75"]  <- .5
    }
    
    if (header == TRUE){
      plot_names <- colnames(ind_mat)[1:rois]
    } else {
      plot_names <- paste0("V", seq(1:rois))
    }

    # functions to restructure matrix to edge list
    W2E   <- function(x) cbind(which(x!=0, arr.ind = TRUE), x[x != 0])
    W2Ena <- function(x) cbind(which(x!=0, arr.ind = TRUE), x[!is.na(x)])
    # set up edge lists for qgraph
    ind_mat_t              <- t(abs(ind_mat))
    Lagged                 <- ind_mat_t[1:(rois),         (rois+1):(rois*2)]
    Contemporaneous        <- ind_mat_t[(rois+1):(rois*2),(rois+1):(rois*2)]
    eLagged                <- W2E(Lagged)
    eContemporaneous       <- W2E(Contemporaneous)
    # set up edge list for colors
    col_mat_t              <- t(col_mat)
    colLagged              <- W2Ena(col_mat_t[1:(rois),          (rois+1):(rois*2)])
    colContemporaneous     <- W2Ena(col_mat_t[(rois+1):(rois*2), (rois+1):(rois*2)])
    color_list             <- c(colLagged[,3], colContemporaneous[,3])
    # denote lagged and contemporaneous
    isLagged               <- c(rep(TRUE,  nrow(eLagged)), 
                                rep(FALSE, nrow(eContemporaneous)))
    # set up curve argument to alter curve of paths
    # if there is a lagged and contemporaneous
    # otherwise only one will be visible
    curve       <- rep(1, length(isLagged))
    final_edges <- rbind(eLagged, eContemporaneous)
    # if an element is repeated, change its curve value
    curve[which(duplicated(final_edges[,1:2]))] <- .5
    
    rownames(final_edges) <- NULL
    test <- data.frame(final_edges, color_list, curve, isLagged, stringsAsFactors = F)
    test <- test[order(test$color_list, decreasing = T), ]
    
    final_edges <- test[,1:3]
    color_list  <- test[,4]
    curve       <- test[,5]
    isLagged    <- test[,6]
    
    rownames(final_edges) <- NULL
    
    as.matrix(final_edges)
    
    parallelEdge <- ifelse(color_list == "black", TRUE, FALSE)
    curve        <- ifelse(color_list == "black", 0, curve)
    fade         <- ifelse(color_list == "black", FALSE, TRUE)
    trans        <- ifelse(color_list == "black", 0, .5)
    
    # set up individual file path
    plotind              <- file.path(input, "plots", paste0(filenames[i], ".pdf"))
    pdf(plotind)
    tryCatch(qgraph::qgraph(as.matrix(final_edges),
                    layout              = "circle",
                    lty                 = ifelse(isLagged, 2, 1),
                    edge.labels         = FALSE,
                    curve               = 1.5,
                    parallelEdge        = parallelEdge,
                    fade                = FALSE, 
                    edge.color          = color_list,
                    mode                = ifelse(weighted == FALSE, "direct", "strength"), 
                    labels              = plot_names,
                    label.cex           = 1.5,
                    edge.label.cex      = 1,
                    edge.label.position = .3),
             error = function(e) e)
    dev.off()
  }
}

