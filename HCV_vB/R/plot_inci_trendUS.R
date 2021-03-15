plot_inci_trendUS=function(){
  library(ggplot2)
  # incidence function , note put in =51
  x=rep(0,81)
  y=rep(0,81)
  
  for ( t in 0:80){
    if (t>=0 && t<=63){
      scaleI = 1
    }else if (t==64){
      scaleI = 1.6;
    }else if (t==65){
      scaleI = 1.8;
    }else if (t==66){
      scaleI = 1.9;
    }else if (t==67){
      scaleI = 2.2;
    }else if (t==68){
      scaleI = 2.3;
    }else if (t==69){
      scaleI = 2.3;
    }else if (t==70){
      scaleI = 2.3;
    }else if (t==71){
      scaleI = 2.2;    
    }else if (t==72){
      scaleI = 2.1;   
    }else if (t==73){
      scaleI = 2.0; 
    }else if (t==74){
      scaleI = 1.9; 
    }else if (t==75){
      scaleI = 1.8;  
    }else if (t==76){
      scaleI = 1.7;  
    }else if (t==77){
      scaleI = 1.6;  
    }else if (t==78){
      scaleI = 1.5;  
    }else if (t==79){
      scaleI = 1.4;  
    }else if (t==80){
      scaleI = 1.3;  
    }else if (t==81){
      scaleI = 1.2;      
    }else{
      scaleI=1;
    }
    x[t+1]=1950+t
    y[t+1]=scaleI
  }
  df = data.frame(x=x,y=y)
  p=ggplot(data=df, aes(x=x, y=y, group=1)) +
    geom_line(color="red",size=2)+
    scale_x_continuous(name="Year", breaks=seq(1950,2030,10))+
    scale_y_continuous(name="Weight", breaks=seq(0,3,0.5),limits=c(0,3))+
    ggtitle("Incidence weighting trend") +
    theme(plot.title = element_text(margin = margin(b = -20)))
  return(p)
}