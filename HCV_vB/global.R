# Library in packages used in this application
library(shiny)
library(DT)
library(RSQLite)
library(DBI)
library(shinyjs)
library(shinycssloaders)
library(shinyWidgets)
library(lubridate)
library(shinyFeedback)
library(dplyr)
library(dbplyr)
library(shinydashboard) 
source('R/parameters_module.R')
source('R/fit_module.R')
source('R/fit_to_data_module.R')
source('R/simple_intervention_module.R')
source('R/runmodel.R')
source('R/runmodel_intervention.R')
conn <- DBI::dbConnect(
  RSQLite::SQLite(),
  dbname = 'www/HCV.db'
)

shiny::onStop(function() {
  dbDisconnect(conn)
})




# Turn off scientific notation
options(scipen = 999)

# Set spinner type (for loading)
options(spinner.type = 8)

# Create 'names_map' dataframe to convert variable names ('names') to clean
# column names ('display_names') in table (i.e. capitalized words, spaces, etc.)
names_map <- data.frame(
  names = c('Vnum', 'Variable', 'Stratum', 'Parameter', 'Units', 'Type', 'Value', 'Ref'),
  display_names = c('Vnum', 'Variable', 'Stratum', 'Parameter', 'Units', 'Type', 'Value', 'Ref'),
  stringsAsFactors = FALSE
)


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

coeffsof_v5=function (){
  # compnames - root name sof the 20 compartments
  # chronic_nams - names of ij ordering
  # chronic_numsij - compartment index for all chronically infected with HCV
  
  
  # Gives the index numbers for categories
  # There are 20*9*4 compartments
  # 20 represents the X notation
  # element 1 to 5 S0,S1,S2,S3,S4
  # infection na?ve or previously achieving spontaneous clearance or SVR through treatment from liver fibrosis stage F0 to F4 respectively
  # elemment 6 is A acute stage
  # element 7 to 11 F0,F1,F2,F3,F4 - the fibrosis stages, chronic
  # element 12 DC stage, chronic
  # element 13 HCC stage, chronic
  # element 14 and 15 LT1 and LT2 - liver transplant stages, chronic
  # element 16 to 20 T0 to T4 
  # chronically infected and in treatment achieving sustained viral response (SVR) (T0 to T4-treated from liver fibrosis stage F0 to F4 respectively)
  compnames=c('S0','S1','S2','S3','S4','A','F0','F1','F2','F3','F4','DC','HCC','LT1','LT2', 'T0','T1','T2','T3','T4')
  chronic_nams=c('00=former,never','01=former,failed','10=current,never','11=current,failed')
  
  # each elemet has three subscripts i,j,k
  # i = 0 formwerPWID, i = 1 current PWID
  # j = 0 never failed treatment, j = 1 had treatment failure
  # k = 1 to 9 the 9 age groups
  
  #for the 20*9*4 compartments
  # ordering is 
  # (i=0,j=0,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) former, never failed
  # (i=0,j=1,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) former, failed
  # (i=1,j=0,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) current, never failed
  # (i=1,j=1,k = 1);(i=0,j=0,k = 2);...;(i=0,j=0,k = 9) current, failed
  
  chronic_numsA =7:15
  chronic_nums00 = chronic_numsA
  chronic_nums01 = 180+chronic_numsA
  chronic_nums10 = 360+chronic_numsA
  chronic_nums11 = 540+chronic_numsA
  
  for (i in 2 : 9){
    chronic_nums00 = c(chronic_nums00 ,20*(i-1)+chronic_numsA)
    chronic_nums01 = c(chronic_nums01 ,180+20*(i-1)+chronic_numsA)
    chronic_nums10 = c(chronic_nums10 ,360+20*(i-1)+chronic_numsA)
    chronic_nums11 = c(chronic_nums11 ,540+20*(i-1)+chronic_numsA)
  }
  
  A_comps_00=seq(6,180,20)
  A_comps_01=A_comps_00+180
  A_comps_10=A_comps_00+360
  A_comps_11=A_comps_00+540
  A_comps=c(A_comps_00,A_comps_01,A_comps_10,A_comps_11)
  
  F0_comps_00 = seq(7,180,20)
  F0_comps_01 = F0_comps_00+180
  F0_comps_10 = F0_comps_00+360
  F0_comps_11 = F0_comps_00+540
  F0_comps=c(F0_comps_00,F0_comps_01,F0_comps_10,F0_comps_11)
  
  F1_comps_00 = seq(8,180,20)
  F1_comps_01 = F1_comps_00+180
  F1_comps_10 = F1_comps_00+360
  F1_comps_11 = F1_comps_00+540
  F1_comps=c(F1_comps_00,F1_comps_01,F1_comps_10,F1_comps_11)
  
  F2_comps_00 = seq(9,180,20)
  F2_comps_01 = F2_comps_00+180
  F2_comps_10 = F2_comps_00+360
  F2_comps_11 = F2_comps_00+540
  F2_comps=c(F2_comps_00,F2_comps_01,F2_comps_10,F2_comps_11)
  
  F3_comps_00 = seq(10,180,20)
  F3_comps_01 = F3_comps_00+180
  F3_comps_10 = F3_comps_00+360
  F3_comps_11 = F3_comps_00+540
  F3_comps=c(F3_comps_00,F3_comps_01,F3_comps_10,F3_comps_11)
  
  F4_comps_00 = seq(11,180,20)
  F4_comps_01 = F4_comps_00+180
  F4_comps_10 = F4_comps_00+360
  F4_comps_11 = F4_comps_00+540
  F4_comps=c(F4_comps_00,F4_comps_01,F4_comps_10,F4_comps_11)
  
  DC_comps_00 = seq(12,180,20)
  DC_comps_01 = DC_comps_00+180
  DC_comps_10 = DC_comps_00+360
  DC_comps_11 = DC_comps_00+540
  DC_comps=c(DC_comps_00,DC_comps_01,DC_comps_10,DC_comps_11)
  
  HCC_comps_00 = seq(13,180,20)
  HCC_comps_01 = HCC_comps_00+180
  HCC_comps_10 = HCC_comps_00+360
  HCC_comps_11 = HCC_comps_00+540
  HCC_comps=c(HCC_comps_00,HCC_comps_01,HCC_comps_10,HCC_comps_11)
  
  LT1_comps_00 = seq(14,180,20)
  LT1_comps_01 = LT1_comps_00+180
  LT1_comps_10 = LT1_comps_00+360
  LT1_comps_11 = LT1_comps_00+540
  LT1_comps=c(LT1_comps_00,LT1_comps_01,LT1_comps_10,LT1_comps_11)
  
  LT2_comps_00 = seq(15,180,20)
  LT2_comps_01 = LT2_comps_00+180
  LT2_comps_10 = LT2_comps_00+360
  LT2_comps_11 = LT2_comps_00+540
  LT2_comps=c(LT2_comps_00,LT2_comps_01,LT2_comps_10,LT2_comps_11)
  
  
  S0_comps_00 = seq(1,180,20)
  S0_comps_01 = S0_comps_00+180
  S0_comps_10 = S0_comps_00+360
  S0_comps_11 = S0_comps_00+540
  S0_comps=c(S0_comps_00,S0_comps_01,S0_comps_10,S0_comps_11)
  
  S1_comps_00 = seq(2,180,20)
  S1_comps_01 = S1_comps_00+180
  S1_comps_10 = S1_comps_00+360
  S1_comps_11 = S1_comps_00+540
  S1_comps=c(S1_comps_00,S1_comps_01,S1_comps_10,S1_comps_11)
  
  S2_comps_00 = seq(3,180,20)
  S2_comps_01 = S2_comps_00+180
  S2_comps_10 = S2_comps_00+360
  S2_comps_11 = S2_comps_00+540
  S2_comps=c(S2_comps_00,S2_comps_01,S2_comps_10,S2_comps_11)
  
  S3_comps_00 = seq(4,180,20)
  S3_comps_01 = S3_comps_00+180
  S3_comps_10 = S3_comps_00+360
  S3_comps_11 = S3_comps_00+540
  S3_comps=c(S3_comps_00,S3_comps_01,S3_comps_10,S3_comps_11)
  
  S4_comps_00 = seq(5,180,20)
  S4_comps_01 = S4_comps_00+180
  S4_comps_10 = S4_comps_00+360
  S4_comps_11 = S4_comps_00+540
  S4_comps=c(S4_comps_00,S4_comps_01,S4_comps_10,S4_comps_11)
  
  T0_comps_00 = seq(16,180,20)
  T0_comps_01 = T0_comps_00+180
  T0_comps_10 = T0_comps_00+360
  T0_comps_11 = T0_comps_00+540
  T0_comps=c(T0_comps_00,T0_comps_01,T0_comps_10,T0_comps_11)
  
  T1_comps_00 = seq(17,180,20)
  T1_comps_01 = T1_comps_00+180
  T1_comps_10 = T1_comps_00+360
  T1_comps_11 = T1_comps_00+540
  T1_comps=c(T1_comps_00,T1_comps_01,T1_comps_10,T1_comps_11)
  
  T2_comps_00 = seq(18,180,20)
  T2_comps_01 = T2_comps_00+180
  T2_comps_10 = T2_comps_00+360
  T2_comps_11 = T2_comps_00+540
  T2_comps=c(T2_comps_00,T2_comps_01,T2_comps_10,T2_comps_11)
  
  T3_comps_00 = seq(19,180,20)
  T3_comps_01 = T3_comps_00+180
  T3_comps_10 = T3_comps_00+360
  T3_comps_11 = T3_comps_00+540
  T3_comps=c(T3_comps_00,T3_comps_01,T3_comps_10,T3_comps_11)
  
  T4_comps_00 = seq(20,180,20)
  T4_comps_01 = T4_comps_00+180
  T4_comps_10 = T4_comps_00+360
  T4_comps_11 = T4_comps_00+540
  T4_comps=c(T4_comps_00,T4_comps_01,T4_comps_10,T4_comps_11)
  # S0, A,F0,F1,F2,F3,F4,DC,HCC,LT1,LT2 - all other comps are empty
  agemat_current=rbind(  # current never failed
    c(361, 366:375),
    c(381 ,386:395),
    c(401 ,406:415),
    c(421 ,426:435),
    c(441 ,446:455),
    c(461 ,466:475),
    c(481 ,486:495),
    c(501 ,506:515),
    c(521 ,526:535))
  
  agemat_former=rbind( # former never failed
    c(361, 366:375),
    c(381 ,386:395),
    c(401 ,406:415),
    c(421 ,426:435),
    c(441 ,446:455),
    c(461 ,466:475),
    c(481 ,486:495),
    c(501 ,506:515),
    c(521 ,526:535))-360
  
  
  return(list(comps=compnames,states=chronic_nams,
              chronic00=chronic_nums00,chronic01=chronic_nums01,chronic10=chronic_nums10,chronic11=chronic_nums11,
              A_comps=A_comps,A_comps_00=A_comps_00,A_comps_01=A_comps_01,A_comps_10=A_comps_10,A_comps_11=A_comps_11,
              F0_comps=F0_comps,F0_comps_00=F0_comps_00,F0_comps_01=F0_comps_01,F0_comps_10=F0_comps_10,F0_comps_11=F0_comps_11,
              F1_comps=F1_comps,F1_comps_00=F1_comps_00,F1_comps_01=F1_comps_01,F1_comps_10=F1_comps_10,F1_comps_11=F1_comps_11,
              F2_comps=F2_comps,F2_comps_00=F2_comps_00,F2_comps_01=F2_comps_01,F2_comps_10=F2_comps_10,F2_comps_11=F2_comps_11,
              F3_comps=F3_comps,F3_comps_00=F3_comps_00,F3_comps_01=F3_comps_01,F3_comps_10=F3_comps_10,F3_comps_11=F3_comps_11,
              F4_comps=F4_comps,F4_comps_00=F4_comps_00,F4_comps_01=F4_comps_01,F4_comps_10=F4_comps_10,F4_comps_11=F4_comps_11,
              DC_comps=DC_comps,DC_comps_00=DC_comps_00,DC_comps_01=DC_comps_01,DC_comps_10=DC_comps_10,DC_comps_11=DC_comps_11,
              HCC_comps=HCC_comps,HCC_comps_00=HCC_comps_00,HCC_comps_01=HCC_comps_01,HCC_comps_10=HCC_comps_10,HCC_comps_11=HCC_comps_11,
              LT1_comps=LT1_comps,LT1_comps_00=LT1_comps_00,LT1_comps_01=LT1_comps_01,LT1_comps_10=LT1_comps_10,LT1_comps_11=LT1_comps_11,
              LT2_comps=LT2_comps,LT2_comps_00=LT2_comps_00,LT2_comps_01=LT2_comps_01,LT2_comps_10=LT2_comps_10,LT2_comps_11=LT2_comps_11,
              S0_comps=S0_comps,S0_comps_00=S0_comps_00,S0_comps_01=S0_comps_01,S0_comps_10=S0_comps_10,S0_comps_11=S0_comps_11,
              S1_comps=S1_comps,S1_comps_00=S1_comps_00,S1_comps_01=S1_comps_01,S1_comps_10=S1_comps_10,S1_comps_11=S1_comps_11,
              S2_comps=S2_comps,S2_comps_00=S2_comps_00,S2_comps_01=S2_comps_01,S2_comps_10=S2_comps_10,S2_comps_11=S2_comps_11,
              S3_comps=S3_comps,S3_comps_00=S3_comps_00,S3_comps_01=S3_comps_01,S3_comps_10=S3_comps_10,S3_comps_11=S3_comps_11,
              S4_comps=S4_comps,S4_comps_00=S4_comps_00,S4_comps_01=S4_comps_01,S4_comps_10=S4_comps_10,S4_comps_11=S4_comps_11,
              T0_comps=T0_comps,T0_comps_00=T0_comps_00,T0_comps_01=T0_comps_01,T0_comps_10=T0_comps_10,T0_comps_11=T0_comps_11,
              T1_comps=T1_comps,T1_comps_00=T1_comps_00,T1_comps_01=T1_comps_01,T1_comps_10=T1_comps_10,T1_comps_11=T1_comps_11,
              T2_comps=T2_comps,T2_comps_00=T2_comps_00,T2_comps_01=T2_comps_01,T2_comps_10=T2_comps_10,T2_comps_11=T2_comps_11,
              T3_comps=T3_comps,T3_comps_00=T3_comps_00,T3_comps_01=T3_comps_01,T3_comps_10=T3_comps_10,T3_comps_11=T3_comps_11,
              T4_comps=T4_comps,T4_comps_00=T4_comps_00,T4_comps_01=T4_comps_01,T4_comps_10=T4_comps_10,T4_comps_11=T4_comps_11,
              agemat_current=agemat_current,agemat_former=agemat_former))
}

age_weighted=function(X,agemat_current,agemat_former){
  # X is the compartments matrix 720 *times
  # agemat_current and agemat_former
  # these are 9 by 11 matrices
  # rows are age group compartment numbers
  # cols are S0,A,F0,F1,F2,F3,F4,DC,HCC,LT1,LT2
  # Assumes 25% of current are in first age group in year 66 (i.e. 2015 assuming 1950:2031 for model)
  # current
  bot = colSums(X[361:540,]); # currently infected
  age_weights_current=rep(9,1)
  for (i in 1 : 9){
    agetot=colSums(X[agemat_current[i,],])
    if (i == 1){
      age_weights_current[i] = 0.25*bot[66]/agetot[66]
    }else{
      age_weights_current[i] = (0.75/8)*bot[66]/agetot[66]
    }
  }
  # former
  bot = colSums(X[1:360,]); # currently infected
  age_weights_former=rep(9,1)
  for (i in 1 : 9){
    agetot=colSums(X[agemat_former[i,],])
    if (i == 1){
      age_weights_former[i] = 0.125*bot[66]/agetot[66]
    }else{
      age_weights_former[i] = (0.875/8)*bot[66]/agetot[66]
    }
  }
  return(list(age_weights_current=age_weights_current,age_weights_former=age_weights_former))
  
  
}
getallpop=function(X,L,age_weights_current,age_weights_former,nt,age_scale00,age_scale10){
  #X is the 720 * 82 matrix
  #L is the list of all the compartment index
  Aformer = matrix(rep(age_weights_former,nt),nrow=9,ncol=nt)
  Apwid = matrix(rep(age_weights_current,nt),nrow=9,ncol=nt)
  with(L,{
    lt2=age_scale00*colSums(Aformer*X[LT2_comps_00,])+age_scale10*colSums(Apwid*X[LT2_comps_10,])+
      age_scale00*colSums(Aformer*X[LT2_comps_01,])+age_scale10*colSums(Apwid*X[LT2_comps_11,])
    lt1=age_scale00*colSums(Aformer*X[LT1_comps_00,])+age_scale10*colSums(Apwid*X[LT1_comps_10,])+
      age_scale00*colSums(Aformer*X[LT1_comps_01,])+age_scale10*colSums(Apwid*X[LT1_comps_11,])
    lt=lt1+lt2;
    dc=age_scale00*colSums(Aformer*X[DC_comps_00,])+age_scale10*colSums(Apwid*X[DC_comps_10,])+
      age_scale00*colSums(Aformer*X[DC_comps_01,])+age_scale10*colSums(Apwid*X[DC_comps_11,])
    hcc=age_scale00*colSums(Aformer*X[HCC_comps_00,])+age_scale10*colSums(Apwid*X[HCC_comps_10,])+
      age_scale00*colSums(Aformer*X[HCC_comps_01,])+age_scale10*colSums(Apwid*X[HCC_comps_11,])
    f4=age_scale00*colSums(Aformer*X[F4_comps_00,])+age_scale10*colSums(Apwid*X[F4_comps_10,])+
      age_scale00*colSums(Aformer*X[F4_comps_01,])+age_scale10*colSums(Apwid*X[F4_comps_11,])
    f3=age_scale00*colSums(Aformer*X[F3_comps_00,])+age_scale10*colSums(Apwid*X[F3_comps_10,])+
      age_scale00*colSums(Aformer*X[F3_comps_01,])+age_scale10*colSums(Apwid*X[F3_comps_11,])
    f2=age_scale00*colSums(Aformer*X[F2_comps_00,])+age_scale10*colSums(Apwid*X[F2_comps_10,])+
      age_scale00*colSums(Aformer*X[F2_comps_01,])+age_scale10*colSums(Apwid*X[F2_comps_11,])
    f1=age_scale00*colSums(Aformer*X[F1_comps_00,])+age_scale10*colSums(Apwid*X[F1_comps_10,])+
      age_scale00*colSums(Aformer*X[F1_comps_01,])+age_scale10*colSums(Apwid*X[F1_comps_11,])
    f0=age_scale00*colSums(Aformer*X[F0_comps_00,])+age_scale10*colSums(Apwid*X[F0_comps_10,])+
      age_scale00*colSums(Aformer*X[F0_comps_01,])+age_scale10*colSums(Apwid*X[F0_comps_11,])
    
    # s compartments
    s4=age_scale00*colSums(Aformer*X[S4_comps_00,])+age_scale10*colSums(Apwid*X[S4_comps_10,])+
      age_scale00*colSums(Aformer*X[S4_comps_01,])+age_scale10*colSums(Apwid*X[S4_comps_11,])
    s3=age_scale00*colSums(Aformer*X[S3_comps_00,])+age_scale10*colSums(Apwid*X[S3_comps_10,])+
      age_scale00*colSums(Aformer*X[S3_comps_01,])+age_scale10*colSums(Apwid*X[S3_comps_11,])
    s2=age_scale00*colSums(Aformer*X[S2_comps_00,])+age_scale10*colSums(Apwid*X[S2_comps_10,])+
      age_scale00*colSums(Aformer*X[S2_comps_01,])+age_scale10*colSums(Apwid*X[S2_comps_11,])
    s1=age_scale00*colSums(Aformer*X[S1_comps_00,])+age_scale10*colSums(Apwid*X[S1_comps_10,])+
      age_scale00*colSums(Aformer*X[S1_comps_01,])+age_scale10*colSums(Apwid*X[S1_comps_11,])
    s0=age_scale00*colSums(Aformer*X[S0_comps_00,])+age_scale10*colSums(Apwid*X[S0_comps_10,])+
      age_scale00*colSums(Aformer*X[S0_comps_01,])+age_scale10*colSums(Apwid*X[S0_comps_11,])
    a0=age_scale00*colSums(Aformer*X[A_comps_00,])+age_scale10*colSums(Apwid*X[A_comps_10,])+
      age_scale00*colSums(Aformer*X[A_comps_01,])+age_scale10*colSums(Apwid*X[A_comps_11,])
    
    t4=age_scale00*colSums(Aformer*X[T4_comps_00,])+age_scale10*colSums(Apwid*X[T4_comps_10,])+
      age_scale00*colSums(Aformer*X[T4_comps_01,])+age_scale10*colSums(Apwid*X[T4_comps_11,])
    t3=age_scale00*colSums(Aformer*X[T3_comps_00,])+age_scale10*colSums(Apwid*X[T3_comps_10,])+
      age_scale00*colSums(Aformer*X[T3_comps_01,])+age_scale10*colSums(Apwid*X[T3_comps_11,])
    t2=age_scale00*colSums(Aformer*X[T2_comps_00,])+age_scale10*colSums(Apwid*X[T2_comps_10,])+
      age_scale00*colSums(Aformer*X[T2_comps_01,])+age_scale10*colSums(Apwid*X[T2_comps_11,])
    t1=age_scale00*colSums(Aformer*X[T1_comps_00,])+age_scale10*colSums(Apwid*X[T1_comps_10,])+
      age_scale00*colSums(Aformer*X[T1_comps_01,])+age_scale10*colSums(Apwid*X[T1_comps_11,])
    t0=age_scale00*colSums(Aformer*X[T0_comps_00,])+age_scale10*colSums(Apwid*X[T0_comps_10,])+
      age_scale00*colSums(Aformer*X[T0_comps_01,])+age_scale10*colSums(Apwid*X[T0_comps_11,])
    
    Tformer=age_scale00*(colSums(Aformer*X[LT2_comps_00,])+colSums(Aformer*X[LT2_comps_01,])+
                           colSums(Aformer*X[LT1_comps_00,])+colSums(Aformer*X[LT1_comps_01,])+
                           colSums(Aformer*X[DC_comps_00,]) +colSums(Aformer*X[DC_comps_01,]) + 
                           colSums(Aformer*X[HCC_comps_00,])+colSums(Aformer*X[HCC_comps_01,])+
                           colSums(Aformer*X[F4_comps_00,]) +colSums(Aformer*X[F4_comps_01,])+
                           colSums(Aformer*X[F3_comps_00,]) +colSums(Aformer*X[F3_comps_01,])+
                           colSums(Aformer*X[F2_comps_00,]) +colSums(Aformer*X[F2_comps_01,])+
                           colSums(Aformer*X[F1_comps_00,]) +colSums(Aformer*X[F1_comps_01,])+
                           colSums(Aformer*X[F0_comps_00,]) +colSums(Aformer*X[F0_comps_01,])+ 
                           colSums(Aformer*X[S4_comps_00,]) +colSums(Aformer*X[S4_comps_01,])+
                           colSums(Aformer*X[S3_comps_00,]) +colSums(Aformer*X[S3_comps_01,])+
                           colSums(Aformer*X[S2_comps_00,]) +colSums(Aformer*X[S2_comps_01,])+
                           colSums(Aformer*X[S1_comps_00,]) +colSums(Aformer*X[S1_comps_01,])+
                           colSums(Aformer*X[S0_comps_00,]) +colSums(Aformer*X[S0_comps_01,])+ 
                           colSums(Aformer*X[T4_comps_00,]) +colSums(Aformer*X[T4_comps_01,])+
                           colSums(Aformer*X[T3_comps_00,]) +colSums(Aformer*X[T3_comps_01,])+
                           colSums(Aformer*X[T2_comps_00,]) +colSums(Aformer*X[T2_comps_01,])+
                           colSums(Aformer*X[T1_comps_00,]) +colSums(Aformer*X[T1_comps_01,])+
                           colSums(Aformer*X[T0_comps_00,]) +colSums(Aformer*X[T0_comps_01,])+  
                           colSums(Aformer*X[A_comps_00,]) +colSums(Aformer*X[A_comps_01,])                         
    )   
    
    Tcurrent=age_scale10*(colSums(Apwid*X[LT2_comps_10,])+colSums(Apwid*X[LT2_comps_11,])+
                            colSums(Apwid*X[LT1_comps_10,])+colSums(Apwid*X[LT1_comps_11,])+
                            colSums(Apwid*X[DC_comps_10,]) +colSums(Apwid*X[DC_comps_11,]) + 
                            colSums(Apwid*X[HCC_comps_10,])+colSums(Apwid*X[HCC_comps_11,])+
                            colSums(Apwid*X[F4_comps_10,]) +colSums(Apwid*X[F4_comps_11,])+
                            colSums(Apwid*X[F3_comps_10,]) +colSums(Apwid*X[F3_comps_11,])+
                            colSums(Apwid*X[F2_comps_10,]) +colSums(Apwid*X[F2_comps_11,])+
                            colSums(Apwid*X[F1_comps_10,]) +colSums(Apwid*X[F1_comps_11,])+
                            colSums(Apwid*X[F0_comps_10,]) +colSums(Apwid*X[F0_comps_11,])+ 
                            colSums(Apwid*X[S4_comps_10,]) +colSums(Apwid*X[S4_comps_11,])+
                            colSums(Apwid*X[S3_comps_10,]) +colSums(Apwid*X[S3_comps_11,])+
                            colSums(Apwid*X[S2_comps_10,]) +colSums(Apwid*X[S2_comps_11,])+
                            colSums(Apwid*X[S1_comps_10,]) +colSums(Apwid*X[S1_comps_11,])+
                            colSums(Apwid*X[S0_comps_10,]) +colSums(Apwid*X[S0_comps_11,])+ 
                            colSums(Apwid*X[T4_comps_10,]) +colSums(Apwid*X[T4_comps_11,])+
                            colSums(Apwid*X[T3_comps_10,]) +colSums(Apwid*X[T3_comps_11,])+
                            colSums(Apwid*X[T2_comps_10,]) +colSums(Apwid*X[T2_comps_11,])+
                            colSums(Apwid*X[T1_comps_10,]) +colSums(Apwid*X[T1_comps_11,])+
                            colSums(Apwid*X[T0_comps_10,]) +colSums(Apwid*X[T0_comps_11,])+  
                            colSums(Apwid*X[A_comps_10,]) +colSums(Apwid*X[A_comps_11,])                         
    )   
    
    
    
    return(list(lt=lt,dc=dc,hcc=hcc,s0=s0,s1=s1,s2=s2,s3=s3,s4=s4,a0=a0,f0=f0,f1=f1,f2=f2,f3=f3,f4=f4,
                t0=t0,t1=t1,t2=t2,t3=t3,t4=t4,Tformer=Tformer,Tcurrent=Tcurrent))
  })
  
}

plot_allcomps=function(times,X,chronic_nums10){
  
  y=100*colSums(X[chronic_nums10,])/colSums(X[361:540,])
  df=data.frame(times,y)
  p1=ggplot(data=df, aes(x=times+1950, y=y)) + geom_line(color="red",size=2) +
    xlab("Year") + ylab("Percentage") + ggtitle("% of current PWID with HCV from 1950 to 2030")+
    scale_y_continuous(breaks = c(seq(0,70,10)),labels=function(x) format(x, big.mark = ",", scientific = FALSE),expand=c(0,0),lim=c(0,70))+
    scale_x_continuous(breaks = c(seq(1950,2031,5)),expand=c(0,0))+
    theme(plot.background = element_rect(fill = "white"),axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #scale_x_discrete(labels=c("2015","2017","2019","2021","2023","2025","2027","2029","2031"))
  return(p1)  
}##BFD5E3

plot_areacomps=function(times,C){
  y=c(C$f0,C$f1,C$f2,C$f3,C$f4,C$dc,C$hcc,C$lt)
  Stage=c(rep("F0",81),rep("F1",81),rep("F2",81),rep("F3",81),rep("F4",81),rep("DC",81),rep("HCC",81),rep("LT",81))
  t=rep(0:80,8)
  df=data.frame(y,t,Stage)
  df$Stage=factor(df$Stage,levels=c("F0","F1","F2","F3","F4","DC","HCC","LT"))
  p1=ggplot(data=df, aes(x=t+1950, y=y,fill=Stage)) + geom_area() +
    xlab("Year") + ylab("Population") + ggtitle("Base model - Number in chronic compartments")+
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE), lim = c(0, 2500000),expand = c(0,0),breaks=seq(0,2500000,250000))+
    scale_x_continuous(breaks = c(seq(1950,2030,5)),expand = c(0,0))+
    theme(plot.background = element_rect(fill = "white"),axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #scale_x_discrete(labels=c("2015","2017","2019","2021","2023","2025","2027","2029","2031"))
  return(p1)  
}

plot_twolines=function(times,fout1,fout2){
  
  y1=fout1
  y2=fout2
  df=data.frame(times,y1,y2)
  p=ggplot(data=subset(df,times>=71), aes(x=times+1950)) + geom_line(aes(y=y1,colour="baseline"),size=2)+
    geom_line(aes(y=y2,colour="intervention"),size=2)+
    scale_colour_manual("", 
                        breaks = c("baseline", "intervention"),
                        values = c("blue", "red")) +
    xlab("Year") + ylab("Number of people") + ggtitle("Number of current PWID with HCV from 2021 to 2030")+
    
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE),expand=c(0,0),lim=c(250000,750000))+
    scale_x_continuous(breaks = c(seq(2021,2030,1)),expand=c(0,0))+
    theme(plot.background = element_rect(fill = "white"),axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #scale_x_discrete(labels=c("2015","2017","2019","2021","2023","2025","2027","2029","2031"))
  return(p)  
}

