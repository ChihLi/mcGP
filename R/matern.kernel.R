matern.kernel <- function(X,para,nu,x=NULL, 
                          para.arg=c("theta", "rho")[1]){
  
  if(para.arg=="theta"){
    if(is.null(x)){
      r <- sqrt(distance(t(t(X)/para)))
    }else{
      r <- sqrt(distance(t(t(X)/para), t(t(x)/para)))
    }
  }else{
    if(is.null(x)){
      r <- sqrt(distance(t(t(X)*sqrt(-log(para)))))
    }else{
      r <- sqrt(distance(t(t(X)*sqrt(-log(para))), t(t(x)*sqrt(-log(para)))))
    }
  }
  
  if(is.null(nu)){ #Gaussian kernel
    out <- exp(-r^2)
  }else if(nu==1/2){ #nu=0.5
    out <- exp(-r)
  }else if(nu==3/2){ #nu=1.5
    out <- (1+r*sqrt(3)) * exp(-r*sqrt(3))
  }else if(nu==5/2){ #nu=2.5
    out <- (1+r*sqrt(5)+5*r^2/3) * exp(-r*sqrt(5))
  }else if(nu==7/2){ #nu=3.5
    out <- (1+r*sqrt(7)+2.8*r^2+7/15*sqrt(7)*r^3) * exp(-r*sqrt(7))
  }else if(nu==9/2){ #nu=4.5
    out <- (1+r*sqrt(9)+27*r^2/7+18/7*r^3+27/35*r^4) * exp(-r*3)
  }
  
  
  return(out)
}