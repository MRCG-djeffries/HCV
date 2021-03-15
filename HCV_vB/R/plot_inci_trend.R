plot_inci_trend=function(){
  library(ggplot2)
  # incidence function , note put in =51
  x=rep(0,81)
  y=rep(0,81)

  for ( t in 0:80){
      if (t>=0 && t<=51){
        scaleI = 1+1.5*t/51;
      }else if (t>51 && t<=56){
        scaleI = 2.5-1.5*(t-51)/5;
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
    ggtitle("Incidence weighting due to drug market activity") +
    theme(plot.title = element_text(margin = margin(b = -20)))
  return(p)
}