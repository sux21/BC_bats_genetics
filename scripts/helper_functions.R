# helper functions for calling allelic peaks from OSIRIS tab delimited analysis result file.  

#'@param x a character string separated by comma or missing value (NA)
#'@retrun a numeric vector
chr_to_num <- function(x) {
  if (is.na(x)) {
    out <- NA
  } else {
    out <- as.numeric(unlist(strsplit(x, ",")))
  }
  return(out)
}

#'@param x a numeric vector or missing value (NA)
#'@return a character string separated by comma
num_to_chr <- function(x) {
  if (unique(is.na(x))) {
    out <- NA
  } else {
    out <- paste(x, collapse = ",")
  }
  return(out)
}

#'@param x a numeric vector or missing value (NA)
#'@param d the maximum difference between two values to be called adjacent values
#'@return integer index of adjacent values if there are adjacent values. Return 0 if no adjacent values. 
find_adjacent_values <- function(x, d) {
  index <- which(diff(x) <= d)
  
  if (length(index) == 0) { 
    index2 <- NA
  } else {
    index2 <- unique(sort(c(index,index+1)))
  }
  
  return(index2)
} 

#'@param x a numeric vector or missing value (NA)
#'@param d the maximum difference between two values to be called adjacent values
#'@return adjacent values which difference is less than or equal to d 
get_adjacent_values <- function(x, d) {
  index <- find_adjacent_values(x, d)
  
  if (length(index) == 0) {
    adj_val <- NA
  } else {
    adj_val <- x[index]
  }
  
  return(adj_val)
}

#'@param x a numeric vector
#'@param d the maximum difference between two values to be called adjacent values
#'@return adjacent values which difference is less than or equal to d 
#'@return integer vector of group id for adjacent values
group_adj_values <- function(x, d) {
  stopifnot(!is.na(x))
  
  group <- c() 
  
  id = 1
  
  adj_val <- get_adjacent_values(x, d)
  
  for (i in 1:length(adj_val)) {
    
    if (i == 1 ) { # first value is always put in group 1
      group <- append(group, 1)
    }
    
    if (i > 1 && i < length(adj_val)) { # middle values
      if (adj_val[i] - adj_val[i-1] <= d) {
        group <- append(group, id)
      } else {
        id = id + 1
        group <- append(group, id) 
      }
    }
    
    if (i == length(adj_val)) { # last value
      if (adj_val[i] - adj_val[i-1] <= d) {
        group <- append(group, id)
      } else {
        id = id + 1
        group <- append(group, id) 
      }
    }
  }
  
  out <- list(adjacent_values = adj_val, group = group)
  return(out)
}

#'@param x a numeric vector
#'@param i a numeric vector contained in x
#'@return integer index of where is i in x
get_value_position <- function(x, i) {
  stopifnot(!is.na(x))
  
  index <- which(x %in% i)
  return(index)
}

#'@param x a numeric vector 
#'@param y a numeric vector with the same number of elements as x
#'@param d the maximum difference between two values to be called adjacent values
#'@return integer vector of group id for adjacent x values
#'@return integer index of maximum y values in each group 
find_adj_x_with_max_y <- function(x, y, d) {
  stopifnot(!is.na(x))
  
  adj_x <- group_adj_values(x, d)
  
  adj_y <- y[find_adjacent_values(x, d)]
    
  max_y_by_group <- aggregate(adj_y ~ adj_x[[2]], FUN = max)
  names(max_y_by_group) <- c("group", "max_y")
  
  out <- list(adj_x = adj_x[[1]]
              , adj_y = adj_y
              , max_y = max_y_by_group$max_y)
  
  return(out)
} 

#'@param x a numeric vector
#'@param y a numeric vector with the same number of elements as x
#'@param d the maximum difference between two values to be called adjacent values
#'@return numeric vector of adjacent x after removing x with maximum y values  
#'@return numeric vector of corresponding y for x above
find_adj_x_with_non_max_y <- function(x, y, d) {
  stopifnot(!is.na(x))
  
  max_y <- find_adj_x_with_max_y(x, y, d)
  
  max_index <- get_value_position(max_y[["adj_y"]], max_y[["max_y"]])
  
  adj_x_non_max <- max_y[["adj_x"]][-max_index]
  
  adj_y_non_max <- max_y[["adj_y"]][-max_index]
  
  out <- list(adj_x_non_max = adj_x_non_max
              , adj_y_non_max = adj_y_non_max)
  
  return(out)
}

#'@param f OSIRIS tab-delimited file
#'@param p patterns of negative control sample name, e.g. p=c("negative", "NEG").Set p=NA if no negative control.
#'@return logical vector indicating which rows correspond to negative control sample
find_negative_ctrl <- function(f, p) {
  neg_sample <- grepl(pattern = paste(p, collapse='|'), f$File.Name)
  return(neg_sample)
} 

#'@param f OSIRIS tab-delimited file
#'@param r row number
#'@param d the maximum size difference between two peaks to be called adjacent peaks
#'@return a character string with adjacent fragments with non-maximum intensity removed
remove_non_max_adjacent_frag <- function (f, r, d) {
  frag_size <- f$Allele[r] |> chr_to_num() # get fragment size and intensity
  
  frag_intensity <- f$RFU[r] |> chr_to_num()
  
  adj_frag_size <- get_adjacent_values(frag_size, d) # get adjacent fragment sizes
  
  adj_frag_intensity <- frag_intensity[find_adjacent_values(frag_size, d)]
  
  if (unique(is.na(adj_frag_size))) { # no adjacent fragments
    
    new_frag_size <- num_to_chr(frag_size) 
    
    new_frag_intensity <- num_to_chr(frag_intensity)
    
    out <- list(new_frag_size, frag_intensity)
    return(out)
    
  } else { # remove adjacent fragments with non-maximum intensity
    
    non_max_adj_frag <- find_adj_x_with_non_max_y(frag_size, frag_intensity, d) # find adjacent fragments with non-maximum intensity
    
    new_frag_size <- frag_size[! frag_size %in% non_max_adj_frag$adj_x_non_max] |> num_to_chr() # remove adjacent fragments with non-maximum intensity
    
    new_frag_intensity <- frag_intensity[! frag_intensity %in% non_max_adj_frag$adj_y_non_max] |> num_to_chr() # remove corresponding intensity of these fragments
    
    out <- list(new_frag_size, new_frag_intensity)
    return(out)
  }  
}



