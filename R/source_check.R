check_njp <- function(x, fit.uncond) {
  
  nm <- "njp"
  if (length(x)) {
    len <- length(fit.uncond$FitList) - 1
    check_integer(x, nm, valid=NULL, min=0, max=len)
  }
  NULL 
}

check_fit.uncond <- function(x, name="fit.uncond") {

  if (!("joinpoint" %in% class(x))) {
    stop(paste0("ERROR: ", name, " must have class 'joinpoint'"))
  }
  if (!is.list(x)) stop(paste0("ERROR: ", name, " must be a list of class 'joinpoint'"))
  obj <- x[["fullpredicted", exact=TRUE]]
  if (is.null(obj)) {
    stop(paste0("ERROR: ", name, " must contain the object 'fullpredicted'"))
  }  

  NULL
}

check_intervals <- function(x1, x2, name1="start.intervals", name2="end.intervals") {

  check_intvec(x1, name1, len=NULL, minvec=1)
  len1 <- length(x1)
  check_intvec(x2, name2, len=len1)
  tmp <- x1 >= x2
  if (any(tmp)) {
    stop(paste0("ERROR: ", name1, " must be < ", name2))
  }
  NULL
}

check_intvec <- function(x, nm, len=NULL, minvec=NULL) {

  lenx <- length(x)
  if (!lenx) stop(paste0("ERROR: ", nm, " has length 0"))
  if (!is.numeric(x) || !is.vector(x)) {
    stop(paste0("ERROR: ", nm, " must be a vector of integers"))
  }
  if (any(x != floor(x))) stop(paste0("ERROR: ", nm, " must be a vector of integers"))
  if (!is.null(len) && (len != lenx)) {
    stop(paste0("ERROR: ", nm, " must have length ", len))
  } 
  minlen <- length(minvec)
  if (minlen == 1) {
    minvec <- rep(minvec, lenx)
    tmp    <- x < minvec
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) {
      stop(paste0("ERROR: ", nm, " must be >= ", minvec[1]))
    }
  }
  NULL
}

check_numvec <- function(x, nm, len=NULL) {

  lenx <- length(x)
  if (!lenx) stop(paste0("ERROR: ", nm, " has length 0"))
  if (!is.numeric(x) || !is.vector(x)) {
    stop(paste0("ERROR: ", nm, " must be a numeric vector"))
  }
  if (!is.null(len) && (len != lenx)) {
    stop(paste0("ERROR: ", nm, " must have length ", len))
  } 
  if (any(!is.finite(x))) {
    stop(paste0("ERROR: ", nm, " must contain finite values"))
  }
  
  NULL
}

check_logical <- function(x, nm) {
  if (!length(x)) stop(paste0("ERROR: ", nm, " must be TRUE or FALSE"))
  if ((x != TRUE) && (x != FALSE)) stop(paste0("ERROR: ", nm, " must be TRUE or FALSE"))
  NULL
}

check_integer <- function(x, nm, valid=NULL, min=NULL, max=NULL) {

  len <- length(x)
  if (len != 1) stop(paste0("ERROR: ", nm, " must be a single integer"))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be a single integer"))
  if (!is.null(valid) && !(x %in% valid)) {
    stop(paste0("ERROR: ", nm, "=", x, " is not valid"))
  }
  if (!is.null(min) &&  (x < min)) {
    stop(paste0("ERROR: ", nm, " must be >= ", min))
  }
  if (!is.null(max) &&  (x > max)) {
    stop(paste0("ERROR: ", nm, " must be <= ", max))
  }

  NULL

}

check_numeric <- function(x, nm, min=NULL, max=NULL) {

  len <- length(x)
  if (len != 1) stop(paste0("ERROR: ", nm, " must be a single numeric value"))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be a numeric value"))
  if (!is.null(min) &&  (x < min)) {
    stop(paste0("ERROR: ", nm, " must be >= ", min))
  }
  if (!is.null(max) &&  (x > max)) {
    stop(paste0("ERROR: ", nm, " must be <= ", max))
  }

  NULL

}

check_formula <- function(x, valid, nm="model.form") {

  if (is.null(x)) return(NULL)
  if (!("formula" %in% class(x))) stop(paste0("ERROR: ", nm, " must be a formula"))
  vars <- all.vars(x)
  if (length(vars)) {
    tmp  <- !(vars %in% valid)
    miss <- vars[tmp]
    str  <- paste0(miss, collapse=", ")
    msg  <- paste0("ERROR in ", nm, ": the variables ", str, " were not found in data")
    stop(msg) 
  }
  NULL
}

check_dataframe <- function(x, nm="data") {

  if (!is.data.frame(x)) stop(paste0("ERROR: ", nm, " must be a data frame"))
  if (nrow(x) < 2) stop(paste0("ERROR: ", nm, " contains too few rows"))
  if (!ncol(x)) stop(paste0("ERROR: ", nm, " contains no columns"))
  NULL

}

check_data <- function(x, nm="data") {
  if (is.data.frame(x)) {
    if (nrow(x) < 1) stop(paste0("ERROR: ", nm, " contains no rows"))
    if (!ncol(x)) stop(paste0("ERROR: ", nm, " contains no columns"))
  } else if (isString(x)) {
    x   <- removeWhiteSpace(x)
    f1  <- paste0(x, c(".dic", ".dic.gz"))
    tmp <- file.exists(f1) 
    if (!any(tmp)) stop(paste0("ERROR: file ", f1[1], " was not found"))
    f2  <- paste0(x, c(".txt", ".txt.gz"))
    tmp <- file.exists(f2) 
    if (!any(tmp)) stop(paste0("ERROR: file ", f2[1], " was not found"))
  } else {
    stop(paste0("ERROR: ", nm, " must be a path to the data or a data frame"))
  }
  x 
}

check_file <- function(x, nm) {

  if (!isString(x)) stop(paste0("ERROR: ", nm, " must be the path to a file"))
  if (!file.exists(x)) stop(paste0("ERROR: ", x, " not found"))
  NULL

}

check_subset <- function(x, len, nm="subset") {

  if (is.null(x)) return(NULL)
  if (is.logical(x)) {
    if (length(x) != len) stop(paste0("ERROR: logical vector ", nm, " must have length ", len, "=nrow(data)"))
    return(NULL)
  }
  if (!isString(x)) stop(paste0("ERROR: ", nm, " must be NULL or a character string defining the subset of observations to include"))
  NULL

}

check_var <- function(x, nm, valid, allow.null=0) {

  if (is.null(x)) {
    if (allow.null) return(NULL)
    stop(paste0("ERROR: ", nm, " must be a column name in data"))
  }

  if (!isString(x)) stop(paste0("ERROR: ", nm, " must be a column name in data"))
  if (!(x %in% valid)) stop(paste0("ERROR: the column ", nm, "='", x, "' not found in data"))

  NULL
}

check_dataVar <- function(data, var, nm, num=1, allow.miss=0) {

  if (!isString(var)) stop(paste0("ERROR: ", nm, " must be a column name in data"))
  estr <- paste0(nm, "=", var)
  if (!(var %in% colnames(data))) stop(paste0("ERROR: ", estr, " not found in data"))
  vec <- getDataVec(data, var)
  if (num && !is.numeric(vec)) stop(paste0("ERROR: ", estr, " must be a numeric column in data"))
  if (!allow.miss && any(!is.finite(vec))) {
   stop(paste0("ERROR: ", estr, " contains missing and/or non-finite values"))
  }
  NULL

}

check_op <- function(op, nm="op") {
 
  if (!length(op)) op <- list()
  if (!is.list(op)) stop(paste0("ERROR: ", nm, " must be NULL or a list"))
  valid <- c("numbetwn", "numfromstart", "numtoend", "print", "DEBUG")
  def   <- list(2, 3, 4, TRUE, FALSE)
  nms   <- names(op)
  if (length(nms)) {
    tmp <- !(nms %in% valid)
    if (any(tmp)) {
      str <- paste0(nms[tmp], collapse=", ")
      msg <- paste0("ERROR: the object(s) '", str, "' in options list ", nm, " are not valid")
      stop(msg)
    }
  }
  op  <- default.list(op, valid, def) 
  for (i in 1:3) check_integer(op[[valid[i]]], paste0(nm, "$", valid[i]), min=0) 
  check_logical(op$print, paste0(nm, "$print"))
  op$conditional <- TRUE # For joinpoint function

  op
}

check_end.interval <- function(x, start.interval, data, interval) {

  vec    <- getDataVec(data, interval)
  maxint <- max(vec, na.rm=TRUE)
  if (length(x)) {
    check_integer(x, "end.interval", valid=NULL, min=start.interval+1, max=maxint)
  } else {
    x <- maxint
  }
  x
}

check_max.cutpoint <- function(x, data, interval) {

  vec    <- getDataVec(data, interval)
  maxint <- max(vec, na.rm=TRUE)
  if (length(x)) {
    check_integer(x, "max.cutpoint", valid=NULL, min=1, max=maxint-1)
  } else {
    x <- maxint-1
  }
  x
}

check_jp.relaxPropObj <- function(x, name="obj") {

  if (!length(x)) stop(paste0("ERROR: length(", name, ") = 0"))
  if (!("jp.relaxProp" %in% class(x))) {
    stop(paste0("ERROR: ", name, " must have class 'jp.relaxProp'"))
  }
  if (!is.list(x)) stop(paste0("ERROR: ", name, " must be a list of class 'jp.relaxProp'"))
  tmp <- x[["all.results", exact=TRUE]]
  if (!length(tmp)) stop(paste0("ERROR: all.results not found in ", name))
  if (!is.list(tmp)) stop("ERROR: all.results must be a list")
  NULL
}

check_add.data.cols <- function(x, name="add.data.cols") {

  if (length(x)) {
    if (!is.character(x)) stop(paste0("ERROR: ", name, " must be character"))
    tmp <- "_ALL_" %in% x
    if (any(tmp)) x <- "_ALL_"
  } else {
    x <- NULL
  }
  x
}

check_legend.pos <- function(x, name="legend.pos") {
  
  if (length(x) != 1) stop(paste0("ERROR: ", name, " must be a string"))
  if (!is.character(x)) stop(paste0("ERROR: ", name, " must be a string"))

  valid <- c("bottomright", "bottom", "bottomleft", 
             "left", "topleft", "top", "topright", 
             "right", "center")
  tmp   <- x %in% valid
  if (!tmp) {
    str <- paste0("'", valid, "'")
    str <- paste0(str, collapse=", ")
    msg <- paste0("ERROR: ", name, " must be one of ", str)
    stop(msg)
  }
  NULL
}
