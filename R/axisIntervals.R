axisIntervals <- function(side, atTick1=monthBegin.5, atTick2=monthEnd.5,
              atLabels=monthMid, labels=month.abb, cex.axis=0.9, ...){
# 1.  Interval start 
  axis(side, at=atTick1, labels=FALSE, ...)
# 2.  Interval end
  if(any(!is.na(atTick2)))axis(side, at=atTick2, labels=FALSE, ...)
# 3.  Interval labels
  axis(side, at=atLabels, labels=labels, tick=FALSE,
       cex.axis=cex.axis, ...)
}


