runmodel_intervention=function(cess_red,relapse_red){
  library(gridExtra)
  library(cowplot)
  library(ggplot2)
  if (cess_red==-1 & relapse_red==-1){
    return(NULL)
  }else{
  # This runs the model
  #  with a percenatge reduction in cessation and relapse rate
  L=coeffsof_v5()
  fitnum = 2
  if (fitnum==1){
    rin=readRDS("www/five_rawUS.rds") # from five parameter fit
  }else{
    rin=readRDS("www/three_rawUS.rds") # form three parameter fit
  }
  N=rep(0,720)
  N[1]=560;N[361]=396;N[367]=44
  param_vals=param_setup_intervention(rin,0,0)
  # the last value in param vals is an intervention indicator
  param_vals=c(param_vals,0) 
  dum_baseline=asode_v2(param_vals,N,L)
  old_vals=c(param_vals[802],param_vals[803])
  paste0("vals are",cess_red)
  param_vals=param_setup_intervention(rin,cess_red,relapse_red)
  # the last value in param vals is an intervention indicator, additionally the original vals of cess and relapse
  
  param_vals=c(param_vals,1,old_vals)
  dum_intervention=asode_v2(param_vals,N,L)
  p=plot_twolines(seq(0, 80, by = 1),dum_baseline$s1,dum_intervention$s1)
  plot(p)
  }
}

asode_v2=function(param_vals,N,L){
  library(deSolve)
  
  
  init       = N
  parameters = param_vals
  times      = seq(0, 80, by = 1)
  lt=length(times)
  X = ode(y=init, times=times, func=MODEL2, parms=parameters,method="ode45")
  X=t(X[,2:721])
  
  scale00=1380000/sum(X[L$chronic00,66]);
  scale10=540000/sum(X[L$chronic10,66]); 
  lt2=scale00*colSums(X[L$LT2_comps_00,])+scale10*colSums(X[L$LT2_comps_10,]);
  lt1=scale00*colSums(X[L$LT1_comps_00,])+scale10*colSums(X[L$LT1_comps_10,]);
  lt=lt1+lt2;
  dc=scale00*colSums(X[L$DC_comps_00,])+scale10*colSums(X[L$DC_comps_10,]);
  hcc=scale00*colSums(X[L$HCC_comps_00,])+scale10*colSums(X[L$HCC_comps_10,]);
  f4=scale00*colSums(X[L$F4_comps_00,])+scale10*colSums(X[L$F4_comps_10,]);
  f3=scale00*colSums(X[L$F3_comps_00,])+scale10*colSums(X[L$F3_comps_10,]);
  f2=scale00*colSums(X[L$F2_comps_00,])+scale10*colSums(X[L$F2_comps_10,]);
  f1=scale00*colSums(X[L$F1_comps_00,])+scale10*colSums(X[L$F1_comps_10,]);
  f0=scale00*colSums(X[L$F0_comps_00,])+scale10*colSums(X[L$F0_comps_10,]);
  

  fout2= (f0[67]+f1[67])/(f0[67]+f1[67]+f2[67]+f3[67]+f4[67]+hcc[67]+dc[67]+lt[67]);
  # fout1 = 100*colSums(X[L$chronic10,])/colSums(X[361:540,])
  fout1 = scale10*colSums(X[L$chronic10,])
  C=getallpop(X,L,rep(1,9),rep(1,9),length(times),scale00,scale10)
  p1=plot_areacomps(times,C)
  p2=plot_allcomps(times,X,L$chronic10)
  
  return(list(s1=fout1,s2=fout2,p1=p1,p2=p2))
}


MODEL2 <- function(time, state, parameters) {
  
  Mvec=parameters[1:400]
  M=t(matrix(Mvec,nrow=20,ncol=20))
  Mdashvec=parameters[401:800]
  Mdash=t(matrix(Mdashvec,nrow=20,ncol=20))
  extra_parms_vals=parameters[801:803]
  phiminusvals=parameters[804:808]
  phiplussvals=parameters[809:813]
  Magevec=parameters[814:894]
  age_matrix=t(matrix(Magevec,nrow=9,ncol=9))
  curr_mort_pwid=parameters[895:903]
  curr_mort_former=parameters[904:912]
  death_rate_dc=parameters[913]
  death_rate_hcc=parameters[914]
  death_rate_lt1=parameters[915]
  death_rate_lt2=parameters[916]
  
  if (parameters[917]==1){# intervention
     if (time<72){
       # 2020
       extra_parms_vals[2]=parameters[918]
       extra_parms_vals[3]=parameters[919]
     }
  }
  
  mort_current = t(matrix(rep(curr_mort_pwid,20),9,20))
  mort_former = t(matrix(rep(curr_mort_former,20),9,20))
  
  Irows = c(6,7,8,9,10,11,12,13,14,15) # rows of infectious PWID
  death_vec= c(-log(1-death_rate_dc), -log(1-death_rate_hcc), -log(1-death_rate_lt1), -log(1- death_rate_lt2)); #DC,HCC,LT1,LT2
  t=time
  N=state
  # percentage PWID infected
  dum=matrix(rep(Irows,each=9),nrow=9)+matrix(seq(360,520,20),nrow=9,ncol=10)
  top10_index = matrix(t(dum),nrow=90,ncol=1)
  dum=matrix(rep(Irows,each=9),nrow=9)+matrix(seq(540,700,20),nrow=9,ncol=10)
  top11_index =  matrix(t(dum),nrow=90,ncol=1)
  top = sum(N[top10_index])+sum(N[top11_index]) # 40 is j = 0, 60 is j = 1 and i = 1 for both
  bot = sum(N[361:length(N)]) # 40 is j = 0, 60 is j = 1 and i = 1 for both
  I = top/bot;
  
  # incidence function , note put in =51
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
  endy=length(phiminusvals)
  phi = scaleI*extra_parms_vals[1]*I
  dum = Mdash[5,5]
  Mdash[phiminusvals[1:(endy-1)]]=-phi
  Mdash[phiminusvals[endy]]=-phi+dum; #dum = -r_svr4DC-rsvr4HCC
  Mdash[phiplussvals]=phi;
  
  X00 = matrix(N[1:(20*9)],20,9)
  X01 = matrix(N[181:(180+20*9)],20,9)
  X10 = matrix(N[361:(360+20*9)],20,9)
  X11 = matrix(N[541:(540+20*9)],20,9)
  
  # with deaths set as zero - this should be closed
  # age movement also set as zero
  # there are treatment so no failed treatment
  # so only D10 and D00
  # only for the youngest age group since no transistion
  d00=M%*%X00-extra_parms_vals[2]*X00+extra_parms_vals[3]*X10-mort_former*X00+t(age_matrix%*%t(X00))
  d01=M%*%X01-extra_parms_vals[2]*X01+extra_parms_vals[3]*X11-mort_former*X01+t(age_matrix%*%t(X01))
  d10=Mdash%*%X10+extra_parms_vals[2]*X00-extra_parms_vals[3]*X10-mort_current*X10+t(age_matrix%*%t(X10))
  d11=Mdash%*%X11+extra_parms_vals[2]*X01-extra_parms_vals[3]*X11-mort_current*X11+t(age_matrix%*%t(X11))
  d10[1,1] =d10[1,1] + sum(mort_former*X00) + sum(mort_former*X01) + sum(mort_current*X10) + sum(mort_current*X11) +
    death_vec[1]*sum(X00[12,])  + death_vec[2]*sum(X00[13,])  + death_vec[3]*sum(X00[14,])  + death_vec[4]*sum(X00[15,]) +
    death_vec[1]*sum(X01[12,])  + death_vec[2]*sum(X01[13,])  + death_vec[3]*sum(X01[14,])  + death_vec[4]*sum(X01[15,]) +
    death_vec[1]*sum(X10[12,])  + death_vec[2]*sum(X10[13,])  + death_vec[3]*sum(X10[14,])  + death_vec[4]*sum(X10[15,]) +
    death_vec[1]*sum(X11[12,])  + death_vec[2]*sum(X11[13,])  + death_vec[3]*sum(X11[14,])  + death_vec[4]*sum(X11[15,])  
  return(list(c(as.vector(d00), as.vector(d01), as.vector(d10), as.vector(d11))))
  
}

param_setup_intervention=function(rin,cess_red,relapse_red){
  if (length(rin) == 5){
    r6=rin[1];
    r5=rin[2];
    r4=rin[3];
    r2=rin[4];
    r3=rin[5];
  }else{
    r6=rin[1];
    r2=rin[2];
    r3=rin[3];
    r5=0.027*(1-relapse_red/100) # default values
    r4=17*(1-cess_red/100)
  }
  
  # This creates two matrices for the static parameters
  # M is for former i=0
  # Mdash is for current i=1, i.e. PWID
  parm_names = c('delta','r_AF0','w','r_SVR4DC','r_SVR4HCC','r_F0F1','r_F1F2','r_F2F3','r_F3F4','r_F4DC','r_F4HCC',
                 'r_DCHCC','r_DCLT','r_DCdeath','r_HCCLT','r_HCCdeath',
                 'r_LT1death','r_LT2death')
  
  parm_current=rep(0,18)
  parm_former=rep(0,18)
  parm_current[1]=0.26; parm_former[1] = parm_current[1]; # spontaneous clearance
  parm_current[2]=52/12; parm_former[2] = parm_current[2];# acute duration is 12 weeks
  parm_current[3]=0; parm_former[3] = parm_current[3];# duration to return to susceptible after treated - no treats for burden model
  parm_current[4]=0; parm_former[4] = parm_current[4];# SVR4 to DC
  parm_current[5]=0; parm_former[5] = parm_current[5];# SVR4 to HCC
  parm_current[6] = -log(1- r2 ); parm_former[6] = -log(1-r3); # F0 to F1 current then former
  parm_current[7] = -log(1- 0.085 ); parm_former[7] = -log(1- 0.074); # F1 to F2 current then former
  parm_current[8] = -log(1-0.085); parm_former[8] = -log(1-0.106); # F2 to F3 current then former
  parm_current[9] = -log(1-0.130); parm_former[9] = -log(1-0.105); # F3 to F4 current then former
  parm_current[10] = -log(1-0.037); parm_former[10] = parm_current[10]; # F4 to DC current then former
  parm_current[11] = -log(1-0.01); parm_former[11] = parm_current[11]; # F4 to HCC current then former
  parm_current[12] = -log(1-0.068); parm_former[12] = parm_current[12]; # DC to HCC
  parm_current[13] = -log(1-0.033); parm_former[13] = parm_current[13]; # DC to LT
  parm_current[15] = -log(1-0.1); parm_former[15] = parm_current[15]; # HCC to LT
  parm_current[14] = -log(1-0.138) ; parm_former[14] = parm_current[14]; # DC to death 
  parm_current[16] = -log(1-0.605); parm_former[16] = parm_current[16]; # HCC death
  parm_current[17] = -log(1-0.169); parm_former[17] = parm_current[17]; # LT to death year 1
  parm_current[18] = -log(1-0.034); parm_former[18] = parm_current[18]; # LT to death year 2
  
  phi = 0;# place holder
  
  Mdash=matrix(data = 0,nrow=20,ncol=20);
  phiminusvals1=1;phiminusvals2=22;phiminusvals3=43;phiminusvals4=64;phiminusvals5=85
  phiplussvals1=6;phiplussvals2=28;phiplussvals3=49;phiplussvals4=70;phiplussvals5=91
  Mdash[1,1]=-phi;Mdash[1,6]=parm_current[1]*parm_current[2];Mdash[1,16]=parm_current[3];
  Mdash[2,2]=-phi;Mdash[2,17]=parm_current[3];
  Mdash[3,3]=-phi;Mdash[3,18]=parm_current[3];
  Mdash[4,4]=-phi;Mdash[4,19]=parm_current[3];
  Mdash[5,5]=-phi-parm_current[4]-parm_current[5];Mdash[5,20]=parm_current[3];
  Mdash[6,1]=phi;Mdash[6,6]=-parm_current[2];
  Mdash[7,6]=(1-parm_current[1])*parm_current[2];Mdash[7,7]=-parm_current[6];
  Mdash[8,2]=phi;Mdash[8,7]=parm_current[6];Mdash[8,8]=-parm_current[7];
  Mdash[9,3]=phi;Mdash[9,8]=parm_current[7];Mdash[9,9]=-parm_current[8];
  Mdash[10,4]=phi;Mdash[10,9]=parm_current[8];Mdash[10,10]=-parm_current[9];
  Mdash[11,5]=phi;Mdash[11,10]=parm_current[9];Mdash[11,11]=-parm_current[10]-parm_current[11];
  Mdash[12,11]=parm_current[10]+parm_current[4];Mdash[12,12] = -parm_current[12]-parm_current[13]-parm_current[14];
  Mdash[13,11]=parm_current[11]+parm_current[5];Mdash[13,12]=parm_current[12];Mdash[13,13]=-parm_current[15]-parm_current[16];
  Mdash[14,12]=parm_current[13];Mdash[14,13]=parm_current[15];Mdash[14,14]=-1-parm_current[17];
  Mdash[15,14]=1;Mdash[15,15]=-parm_current[18];
  Mdash[16,16]=-parm_current[3];
  Mdash[17,17]=-parm_current[3];
  Mdash[18,18]=-parm_current[3];
  Mdash[19,19]=-parm_current[3];
  Mdash[20,20]=-parm_current[3];
  
  
  M=matrix(data = 0,nrow=20,ncol=20);
  M[1,1]=-phi;M[1,6]=parm_former[1]*parm_former[2];M[1,16]=parm_former[3];
  M[2,2]=-phi;M[2,17]=parm_former[3];
  M[3,3]=-phi;M[3,18]=parm_former[3];
  M[4,4]=-phi;M[4,19]=parm_former[3];
  M[5,5]=-phi-parm_former[4]-parm_former[5];M[5,20]=parm_former[3];
  M[6,1]=phi;M[6,6]=-parm_former[2];
  M[7,6]=(1-parm_former[1])*parm_former[2];M[7,7]=-parm_former[6];
  M[8,2]=phi;M[8,7]=parm_former[6];M[8,8]=-parm_former[7];
  M[9,3]=phi;M[9,8]=parm_former[7];M[9,9]=-parm_former[8];
  M[10,4]=phi;M[10,9]=parm_former[8];M[10,10]=-parm_former[9];
  M[11,5]=phi;M[11,10]=parm_former[9];M[11,11]=-parm_former[10]-parm_former[11];
  M[12,11]=parm_former[10]+parm_former[4];M[12,12] = -parm_former[12]-parm_former[13]-parm_former[14];
  M[13,11]=parm_former[11]+parm_former[5];M[13,12]=parm_former[12];M[13,13]=-parm_former[15]-parm_former[16];
  M[14,12]=parm_former[13];M[14,13]=parm_former[15];M[14,14]=-1-parm_former[17];
  M[15,14]=1;M[15,15]=-parm_former[18];
  M[16,16]=-parm_former[3];
  M[17,17]=-parm_former[3];
  M[18,18]=-parm_former[3];
  M[19,19]=-parm_former[3];
  M[20,20]=-parm_former[3];
  
  age_matrix=matrix(data = 0, nrow=9,ncol=9)
  age_matrix[1,1] = -1/5;
  age_matrix[2,1] = 1/5;age_matrix[2,2] = -1/5;
  age_matrix[3,2] = 1/5;age_matrix[3,3] = -1/5;
  age_matrix[4,3] = 1/5;age_matrix[4,4] = -1/10;
  age_matrix[5,4] = 1/10;age_matrix[5,5] = -1/10;
  age_matrix[6,5] = 1/10;age_matrix[6,6] = -1/10;
  age_matrix[7,6] = 1/10;age_matrix[7,7] = -1/10;
  age_matrix[8,7] = 1/10;age_matrix[8,8] = -1/10;
  age_matrix[9,8] = 1/10;
  
  curr_mort_pwid=c(0.96 ,0.96 ,1.12 ,0.18 ,0.22 ,0.53 ,1.38 ,4.28 ,14.96)/1000; # mortality per year for PWID
  curr_mort_former=c(0.044 ,0.051 ,0.062 ,0.1 ,0.222 ,0.534 ,1.376 ,4.282 ,14.956 )/1000; # mortality per year for PWID
  
  
  mort_current = t(matrix(rep(curr_mort_pwid,20),ncol=20))
  mort_former = t(matrix(rep(curr_mort_former,20),ncol=20))
  extra_parms_nams=c('piv','relapse','nu'); # infection rate (piv instead of pi), relapse to IDU, 1/duration of injecting span
  
  extra_parms_vals=c(r6,-log(1-r5),1/r4) # was 5.6%, nu was 1/17
  rout=c(r6,r5,r4,r2,r3)
  
  # starting values t = 1950:2030
  # S001 = 560 # formwer PWID compartment 161
  # S101 = 400 # current PWID compartment 361
  # F101 = 40 # acutely infected with HCV  366 (A is comp 6 in X list)
  N0=matrix(nrow=20*9*4,ncol=1);
  N0[1]=560;
  N0[361]=396;
  N0[367]=44;
  Mvec=as.vector(t(M))
  Mdashvec=as.vector(t(Mdash))
  Magevec=as.vector(t(age_matrix))
  
  death_rate_dc=  -log(1-0.138) 
  death_rate_hcc= -log(1-0.605) 
  death_rate_lt1= -log(1-0.169) 
  death_rate_lt2= -log(1-0.034) 
  
  param_vals=c(Mvec,Mdashvec,extra_parms_vals,
               phiminusvals1,phiminusvals2,phiminusvals3,phiminusvals4,phiminusvals5,
               phiplussvals1,phiplussvals2,phiplussvals3,phiplussvals4,phiplussvals5,
               Magevec,
               curr_mort_pwid,curr_mort_former,
               death_rate_dc,death_rate_hcc,death_rate_lt1,death_rate_lt2)
  
  return(param_vals)
}

