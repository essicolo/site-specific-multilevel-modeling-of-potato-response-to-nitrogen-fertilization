ilrDefinition <- function(sbp, side="+-", sep.elem = ",", sep.bal = " | ", sep.left = "[", sep.right = "]") {
  
  if (nrow(sbp) != (ncol(sbp)-1)) stop("SBP not valid")
  
  ilrDef <- vector()
  for (n in 1:nrow(sbp)) {
    pos <- names(sbp[n,][which(sbp[n,] == 1)])
    neg <- names(sbp[n,][which(sbp[n,] == -1)])
    if (side=="-+") {
      pos <- rev(pos)
      neg <- rev(neg)
    }
    pos.group <- character()
    neg.group <- character()
    for (i in 1:length(pos)) {
      if (i == 1) {
        pos.group <- paste(pos.group,pos[i], sep="")
      } else {
        pos.group <- paste(pos.group,pos[i], sep=sep.elem)
      }
    }
    for (i in 1:length(neg)) {
      if (i == 1) {
        neg.group <- paste(neg.group,neg[i], sep="")
      } else {
        neg.group <- paste(neg.group,neg[i], sep=sep.elem)
      }
    }
    if (side=="+-") {
      ilrDef[n] <- paste(sep.left,pos.group, sep.bal,neg.group, sep.right, sep="")
    } else if (side=="-+") {
      ilrDef[n] <- paste(sep.left,neg.group, sep.bal,pos.group, sep.right, sep="")
    }
    
  }
  
  ilrDef
}
