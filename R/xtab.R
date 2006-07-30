
DWIMxtabs<-function (formula, data = parent.env(), var.labels=NULL, ...) 
{
    res = list()
    
    mf <- model.frame(formula=formula,data=data)
    contin <- which(!sapply(mf,is.discrete))
    resp <- attr(terms(formula,data=data),"response")

    contin <- setdiff(contin,resp)
    ncat <- sapply(mf[contin], function(x)length(unique(x)))
    ncuts <- ifelse(ncat<5, ncat, pmin(25,round(5+(ncat-5)/10)))
    if (length(contin)>0) {
       for (i in 1:length(contin)) {
          mf[contin[i]] <- cut(mf[[contin[i]]],ncuts[i])
        }
    }
    res[[1]] = do.call("xtabs", list(formula = formula, data = mf, 
        ...))
    if (is.null(var.labels)) {
      # apply variable labels from attribute if present
       var.labels=attr(data,"var.labels")
    } 
    if (!is.null(var.labels)) {
	tmpdim=attr(res[[1]],"dimnames")
	names(tmpdim)=paste(names(tmpdim),
		var.labels[match(names(attr(res[[1]],"dimnames")),names(data))]
		,sep=": ")
	#attr(res[[1]],"dimnames")=tmpdim
	attr(res[[1]],"varNameWLbl")<-names(tmpdim)
    }
    attr(res, "ftable") = ft = ftable(res[[1]])
    attr(res, "stats") = summary(res[[1]])
    
	if (length(dim(res[[1]]))== 0) {
		attr(res, "rowTotals")    <- apply(res[[1]], 1, sum)
		attr(res, "colTotals")    <- apply(res[[1]], 2, sum)
		attr(res, "total")        <- sum(attr(res, "colTotals"))
		attr(res, "rowPercent")   <- prop.table(res[[1]], 1)
		attr(res, "colPercent")   <- prop.table(res[[1]], 2)
		attr(res, "totalPercent") <- (res[[1]] / attr(res,"total"))
	} else if (length(dim(res[[1]])) >= 2){
		attr(res, "rowTotals")    <- apply(ft, 1, sum)
		attr(res, "colTotals")    <- apply(ft, 2, sum)
		attr(res, "total")        <- sum(attr(res, "colTotals"))
		attr(res, "rowPercent")   <- (prop.table(ft, 1))
		attr(res, "colPercent")   <- (prop.table(ft, 2))
		attr(res, "totalPercent") <- (ft / attr(res,"total"))

	}
    class(res) = "DWIMxtabs"
    return(res)
}

print.DWIMxtabs<-function(x,...,
    wantPercentages=T, wantTotals=T, wantStats=T,nameLength=15) {
      ft = attr(x,"ftable")

      tmpat = attributes(ft)
      # handle long labels
      longnames=which(nchar(names(tmpat$col.vars))>nameLength)
      names(tmpat$col.vars)[longnames]=
        paste(substr(names(tmpat$col.vars)[longnames],1,10),"",sep="...")
      longnames=which(nchar(names(tmpat$row.vars))>nameLength)
      names(tmpat$row.vars)[longnames]=
        paste(substr(names(tmpat$row.vars)[longnames],1,10),"",sep="...")

		nan2zero <- function(v){if (is.nan(v)) {0} else {v}}

		makeLabels <- function(lst) {
			lens <- sapply(lst, length)
			cplensU <- c(1, cumprod(lens))
			cplensD <- rev(c(1, cumprod(rev(lens))))
			y <- NULL
			for (i in rev(seq(along = lst))) {
				ind <- 1 + seq(from = 0, to = lens[i] - 1) * cplensD[i + 1]
				tmp <- character(length = cplensD[i])
				tmp[ind] <- lst[[i]]
				y <- cbind(rep(tmp, times = cplensU[i]), y)
			}
			y
		}

		makeNames <- function(x) {
			nmx <- names(x)
			if(is.null(nmx)){
				nmx <- rep("", length.out = length(x))
			}
			nmx
		}

	if (wantPercentages) {
		
		fmtttl <- paste(sapply(round(attr(x,"totalPercent")*100,digits=1),nan2zero),"",sep="")
		
		fmtrow <- paste(sapply(round(attr(x,"rowPercent")*100,digits=1),nan2zero),"",sep="")
		fmtcol <- paste(sapply(round(attr(x,"colPercent")*100,digits=1),nan2zero),"",sep="")
		ftfmt  <- paste(ft, paste("(", paste(fmtttl,fmtrow,fmtcol,sep=", "),")",sep=""), sep="")
		
		dim(ftfmt)=dim(ft)
		

	} else {
		ftfmt =   unclass(ft)
	}

	if (wantTotals) {
		rwp <- attr(x,"rowTotals")/ attr(x,"total")
		clp <- attr(x,"colTotals")/ attr(x,"total")
		if (wantPercentages) {
		   mprow <- paste("(",round(attr(x,"rowTotals")/ attr(x,"total")*100,digits=1),"%)",sep="")
		   mpcol <- paste("(",round(attr(x,"colTotals")/ attr(x,"total")*100,digits=1),"%)",sep="")  
       	mgr <- paste(attr(x,"rowTotals"), mprow, sep="")
		   mgc <- paste(attr(x,"colTotals"), mpcol, sep="") 
		}   else   {
        mgr <- attr(x,"rowTotals")
        mgc <- attr(x,"colTotals")
    }  

    if (wantPercentages) {
		   ftfmt <- rbind(cbind(ftfmt, unclass(mgr)), c(unclass(mgc), paste(attr(x,"total"),'(100%)',sep="")))
    } else {
      ftfmt <- rbind(cbind(ftfmt, unclass(mgr)), c(unclass(mgc), paste(attr(x,"total"),sep="")))
    }
		#ftfmt <- rbind(cbind(ftfmt, RowTotal=mgr), ColTotal=c(mgc, paste(attr(x,"total"),'(100%)',sep="")))
	}

	xrv <- attr(ft, "row.vars")
	xcv <- attr(ft, "col.vars")

	mraw <- rbind(matrix(data="", nrow=1, ncol=(1+dim(ftfmt)[2])), cbind(matrix(data="", nrow=dim(ftfmt)[1], ncol=1), ftfmt))
	lhsm <- rbind(matrix("", nr = length(xcv), nc = length(xrv)), makeNames(xrv),makeLabels(xrv),if (wantTotals) {c("ColTotal", rep("",times= (length(xrv)-1)))})
	rhstop<-cbind(c(makeNames(xcv)), if(length(xcv)) { t(makeLabels(xcv))}, if (wantTotals) {c("RowTotal")})
	rhsm<-rbind(rhstop, mraw)

	m<-cbind(apply(lhsm, 2, format, justify = "left"), apply(rhsm, 2, format, justify = "right"))
	if (!is.null(attr(x[[1]],"varNameWLbl"))) {
		varSet<-paste(attr(x[[1]],"varNameWLbl"),collapse=" and ",sep="")
	} else {
		varSet<-paste(names(dimnames(x[[1]])), collapse=", ")
	}
	cat("\nDWIMxtabs among ",varSet,"\n\n", sep="")

	cat(t(m),sep = c(rep(" ", ncol(m) - 1), "\n"))
	cat("\nDWIMxtab Cell Legend: Frequency (total %, row %, column %)\n\n", sep="")
      
      if (wantStats) {
         tmpc=x[[1]]; attr(tmpc,"call")=NULL
         print(summary(tmpc),...)
      }
      
      return(invisible(ft))
  }

  HTML.DWIMxtabs<-function(x,...,
    wantPercentages=T, wantTotals=T, wantStats=T, nameLength=15) {
      ft = attr(x,"ftable")

      tmpat = attributes(ft)

      # handle long labels
      longnames=which(nchar(names(tmpat$col.vars))>nameLength)
      names(tmpat$col.vars)[longnames]=
        paste(substr(names(tmpat$col.vars)[longnames],1,10),"",sep="...")
      longnames=which(nchar(names(tmpat$row.vars))>nameLength)
      names(tmpat$row.vars)[longnames]=
        paste(substr(names(tmpat$row.vars)[longnames],1,10),"",sep="...")

	nan2zero <- function(v){if (is.nan(v)) {0} else {v}}

	makeLabels <- function(lst) {
		lens <- sapply(lst, length)
		cplensU <- c(1, cumprod(lens))
		cplensD <- rev(c(1, cumprod(rev(lens))))
		y <- NULL
		for (i in rev(seq(along = lst))) {
			ind <- 1 + seq(from = 0, to = lens[i] - 1) * cplensD[i + 1]
			tmp <- character(length = cplensD[i])
			tmp[ind] <- lst[[i]]
			y <- cbind(rep(tmp, times = cplensU[i]), y)
		}
		y
	}

	makeNames <- function(x) {
		nmx <- names(x)
		if(is.null(nmx)){
			nmx <- rep("", length.out = length(x))
		}
		nmx
	}


	if (wantPercentages) {
	
		fmtttl <-      paste("<span class='DWIMrowper'>(",sapply(round(attr(x,"totalPercent")*100,digits=1),nan2zero),"%)</span>",sep="")
		fmtrow <- paste("<br/><span class='DWIMper'>(",sapply(round(attr(x,"rowPercent")*100,digits=1),nan2zero),"%)</span>",sep="")
		fmtcol <- paste("<br/><span class='DWIMrowper'>(",sapply(round(attr(x,"colPercent")*100,digits=1),nan2zero),"%)</span>",sep="")

		ftfmt<-paste(unclass(ft), fmtttl, fmtrow, fmtcol)

		dim(ftfmt)<-dim(ft)
	} else {
		ftfmt <- unclass(ft)
	}

	if (wantTotals) {
		
		rwp <- attr(x,"rowTotals")/ attr(x,"total")
		clp <- attr(x,"colTotals")/ attr(x,"total")
		mprow <- paste("<span class='DWIMrowper'>(",round(attr(x,"rowTotals")/ attr(x,"total")*100,digits=1),"%)</span>",sep="")
		mpcol <- paste("<span class='DWIMcolper'>(",round(attr(x,"colTotals")/ attr(x,"total")*100,digits=1),"%)</span>",sep="")
		mgr <- paste(attr(x,"rowTotals"), mprow, sep="")
		mgc <- paste(attr(x,"colTotals"), mpcol, sep="")

		#ftfmt <- cbind(ftfmt, RowTotal=mgr)
		#ftfmt <- rbind(cbind(ftfmt, RowTotal=mgr), ColTotal=c(mgc, paste(attr(x,"total"),'(100%)',sep="")))
		ftfmt <- rbind(cbind(ftfmt, unclass(mgr)), c(unclass(mgc), paste(attr(x,"total"),'(100%)',sep="")))

	}
	

	xrv <- attr(ft, "row.vars")
	xcv <- attr(ft, "col.vars")
	
	
	if (length(dim(x[[1]])) == 2){
		m <- rbind(matrix(data="", nrow=1, ncol=(1+dim(ftfmt)[2])),cbind(matrix(data="", nrow=dim(ftfmt)[1], ncol=1), ftfmt))
		dmnm<-list();
		dmnm[[1]]<-c(names(dimnames(x[[1]]))[1], dimnames(x[[1]])[[1]], if (wantTotals) {"ColTotal"})
		dmnm[[2]]<-c(names(dimnames(x[[1]]))[2], dimnames(x[[1]])[[2]], if (wantTotals) {"RowTotal"})
		dimnames(m)<-dmnm
	} else if (length(dim(x[[1]])) > 2){
		mraw <- rbind(matrix(data="", nrow=1, ncol=(1+dim(ftfmt)[2])), cbind(matrix(data="", nrow=dim(ftfmt)[1], ncol=1), ftfmt))
		lhsm <- rbind(matrix("", nr = length(xcv), nc = length(xrv)), makeNames(xrv),makeLabels(xrv),if (wantTotals) {c("ColTotal", rep("",times= (length(xrv)-1)))})
		rhstop<-cbind(c(makeNames(xcv)), if(length(xcv)) { t(makeLabels(xcv))}, if (wantTotals) {c("RowTotal")})
		rhsm<-rbind(rhstop, mraw)
		m<-cbind(lhsm,rhsm)
	}

	if (!is.null(attr(x[[1]],"varNameWLbl"))) {
		varSet<-paste(attr(x[[1]],"varNameWLbl"),collapse=" and ",sep="")
	} else {
		varSet<-paste(names(dimnames(x[[1]])), collapse=", ")
	}

	tblcaption<-paste("<h3>DWIMxtabs among ",varSet,"</h3>\n", sep="")
	tbllegend<-paste("<blockquote>DWIMxtab Cell Legend: Frequency (total %, row %, column %)</blockquote>\n", sep="")
	
	HTML(tblcaption,...)
	
	HTML(m,...)
	HTML(tbllegend, ...)
	HTML("<br/><br/>",...)

      if (wantStats) {
         tmpc=x[[1]]; attr(tmpc,"call")=NULL
         HTML(summary(tmpc),...)
      }

      return(invisible(ft))
  }
  
