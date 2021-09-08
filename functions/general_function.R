# column-wise zero to one scaling of elements in a matrix
zero_one_scale<-function(datt){
  for (d in 1:ncol(datt)){
    dn<-datt[,d]
      maxx<-max(dn,na.rm = T)
      minn<-min(dn,na.rm = T)
      delta<-maxx-minn
      datt[,d]<-(datt[,d]-minn)/delta
    
  } #end for
  return(datt)
} # end function
