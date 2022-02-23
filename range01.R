range01 = function(x){
  if (all(x == 0)){
    return(x)
  }else{
    return((x-min(x))/(max(x)-min(x)))
  }
}