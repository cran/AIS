######################################################
#       
# is.discrete
#       
# Attempts to evaluate whether vector contains discrete values
#
# Parameters:
#
# x - a vector
# 
######################################################

is.discrete<-function(x) {
	if (is.data.frame(x) || is.list(x)) {
		return(sapply(x,is.discrete))
	}
	x=as.vector(x);
	if (is.integer(x)) {
		return (TRUE);
	}  else if (is.factor(x)) {
		return (TRUE);
  } else if (is.logical(x)) {
    return (TRUE);
  } else if (is.character(x)) {
    return(TRUE);
	} else  if (sum(as.integer(x)!=x)==0) {
		return (TRUE);
	}
	return (FALSE);
}

is.nominal<-function(x) {
  if (is.data.frame(x) || is.list(x)) {
		return(sapply(x,is.nominal))
	}
	if (is.character(x) || is.logical(x) || (is.factor(x) && !is.ordered(x))) {
	 return(TRUE)
  } 
  return(FALSE)
}