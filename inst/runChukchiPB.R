#### analyze data from polar bears and their tracks

data(Chukchi_PB_data)

Area_table = Area_table[-757,]  #removing 1 record with little area surveyed and no grid ID

##extracting coordinates from geom 
grid<-cbind(Grid_list[[1]], st_coordinates(st_centroid(Grid_list[[1]])))
colnames(grid)[colnames(grid)=="X"] <- "easting"
colnames(grid)[colnames(grid)=="Y"] <- "northing"

#format russian data
Which_obs = which(pb_dist$distance<600)
Distances = pb_dist$distance[Which_obs]

#Convert Chukotka time and date to Alaskan time/dat

d1<-strptime(pb_ru$Time_date, format="%m/%d/%Y %H:%M", tz = "Asia/Anadyr")
pb_ru$correct_dt <- as.POSIXct(d1, tz="Asia/Anadyr")  
attributes(pb_ru$correct_dt)$tzone <- "US/Alaska" 

I_us = I_regehr = I_PBSG = rep(0,1354)
I_us[I_US]=1
Regehr_study_area = st_transform(Regehr_study_area,crs=st_crs(grid))
PBSG = st_transform(PBSG,area,crs=st_crs(grid))
I_regehr = 1*(is.na(as.numeric(st_within(st_centroid(grid),Regehr_study_area)))==FALSE)
I_PBSG = 1*(is.na(as.numeric(st_within(st_centroid(grid),PBSG)))==FALSE)

#function to drop the geometry of an st object (easier to work w/ the data.frame)
st_drop_geometry <- function(x) {
  if(inherits(x,"sf")) {
    x <- st_set_geometry(x, NULL)
    class(x) <- 'data.frame'
  }
  return(x)
}

# #calculate number of photographs, number with tracks for each unique time & grid cell surveyed
tracks_df = st_drop_geometry(tracks)
jday = floor(julian(tracks_df[,"correct_dt"], origin = as.POSIXct("2016-01-01", tz = "US/Alaska")))
tracks$jday = jday
tracks$day = as.numeric(jday - min(jday))+1
tracks$cell_time = paste(Cell_lookup_table[tracks_df$GridID],tracks$day)
tracks$I_track = as.numeric(tracks_df[,"detect_bear_track"]=="Y")
tracks_df = st_drop_geometry(tracks)
n_cells = nrow(Grid_list[[1]])
t_steps = max(tracks$day)
n_st = n_cells*t_steps

Unique_ct = paste(Cell_lookup_table[Area_table[,"Grid_ID"]],Area_table[,"Day"])
Unique_tracks = unique(tracks_df[,"cell_time"])
n_i = length(Unique_ct)
N_i = T_i = S_i = C_i = U_i = rep(0,n_i)
for(i in 1:n_i){
  Which = which(tracks_df[,"cell_time"]==Unique_ct[i])
  N_i[i]=length(Which)
  if(N_i[i]>0)T_i[i]=sum(tracks_df[Which,"I_track"])  ### number of pics with tracks on them per grid cell (?)
  S_i[i]=(Area_table[i,"Day"]-1)*n_cells+Cell_lookup_table[Area_table[i,"Grid_ID"]]
}

Grid_df = st_drop_geometry(grid)
##match up Rus data with big Chukchi grid
rus_cell_lookup_table<-c()

for (i in 1:nrow(pb_ru)){
  Which_id = which(Grid_df[,"objectid"]==pb_ru[i,"objectd"])
  rus_cell_lookup_table[i]<-Which_id
}


#calculate area surveyed, and number of tracks per km transect for each unique time & grid cell surveyed in RUSSIA
#tracks_df = st_drop_geometry(tracks)
jday_ru = floor(julian(pb_ru[,"correct_dt"], origin = as.POSIXct("2016-01-01", tz = "US/Alaska")))
pb_ru$jday = jday_ru
pb_ru$day = as.numeric(jday_ru - min(jday))+1
pb_ru$cell_time = paste(rus_cell_lookup_table,pb_ru$day)
I_may = (pb_ru$day>24)  #indicator for which russian survey days were in May


#for visually detected tracks
Unique_ct_ru = paste(rus_cell_lookup_table,pb_ru[,"day"])
Unique_tracks_ru = unique(pb_ru[,"cell_time"])
n_i_ru = length(Unique_ct_ru)
T_i_ru = S_i_ru = C_i_ru = U_i_ru = KM_i_ru = rep(0,n_i_ru) ##N_i_ru =
for(i in 1:n_i_ru){
  Which_ru = which(pb_ru[,"cell_time"]==Unique_ct_ru[i])
  KM_i_ru[i]=sum(pb_ru[Which_ru, "LENGTH"])/1000  ##transect length
  if(KM_i_ru[i]>0)T_i_ru[i]=sum(pb_ru[Which_ru,"N_tracks"]) ####estimating number of tracks per grid cell
  S_i_ru[i]=(pb_ru[i,"day"]-1)*n_cells+rus_cell_lookup_table[i]   ##site-time index for each sample
}

#indicator for sampled/not sampled
Y_s = matrix(0,n_cells,t_steps)
Y_s[S_i]=1


# #bear counts US
CellDay_on = c("457 35","1023 40","415 7")
Which = rep(0,3)
Which = c()
for(i in 1:3)Which[i]=which(Unique_ct==CellDay_on[i])
C_i[Which]=1


CellDay_off = c("1229 15","1309 9","1341 9","1342 9")

Which = rep(0,4)
for(i in 1:4)Which[i]=which(Unique_ct==CellDay_off[i])
U_i[Which]=1

#bear counts for RUSSIA (photographed ADULTS)
  #for pb_ru[,"Adults"]>0 use pb_ru[,"cell_time"]
C_i_ru = pb_ru$Photo
C_i_ru[is.na(C_i_ru)]=0
Tmp = pb_ru$Visual_all
Tmp[is.na(Tmp)]=0
C_i_ru = C_i_ru+Tmp
#remove duplicates per Irina 10/16/2019
C_i_ru[c(160,250,286,327,439,510)]=C_i_ru[c(160,250,286,327,439,510)]-1


###group size vector

##us portion

G_i<-c(1,1,1,1,1,1,1,3)
#rus portion  - all IR, photo and visual detections  ##not sure if we want all detections or instrument based only?

G_i_ru<-c(3,	2,	1,	1,	1,	1,	3,	1,	1,	1,	1,	1,	2,	3,	1,	2,	1,	1,	1,	1,	3,	1,	1,	1,	1,	2,	3,	1,	1,	1,	1,	1,	1,	1,	1,	2,	1,	1,	1,	1,	2,	2,	1,	2,	1,	2,	1,	1,	1)
#rus instuments only (IR and photo)
#G_i_ru<-c(3,	2,	1,	1,	1,	1,	3,	1,	1,	1,	1,	1,	2,	3,	1,	2,	1,	1,	1,	1,	3,	1,	1,	1,	1,	2,	3,	1,	1,	1)

G<-c(G_i, G_i_ru)


# # assemble design matrix for spatial covariates
Grid_df = st_drop_geometry(grid)
RSF = rep(Grid_df[,"rsf"],t_steps)
dist_land = rep(Grid_df[,"dist_land"],t_steps)
dist_land = dist_land/mean(dist_land)
dist_land2 = dist_land^2
Sea_ice = I_water = rep(0,n_st)
for(it in 1:t_steps){
  Sea_ice[n_cells*(it-1)+c(1:n_cells)]=st_drop_geometry(Grid_list[[it]])[,"sea_ice"]
  I_water[n_cells*(it-1)+c(1:n_cells)]=(Sea_ice[n_cells*(it-1)+c(1:n_cells)]<0.01)
}
Coords = st_coordinates(st_centroid(grid))
Coords = t(t(Coords)/colMeans(Coords))
easting = rep(Coords[,1],t_steps)
northing = rep(Coords[,2],t_steps)

Xt_s = data.frame("Intercept"=rep(1,n_st),"iwater"=I_water)
Xt_s[,"dist_land"]=dist_land
Xt_s[,"dist_land2"]=dist_land2
Xt_s[,"ice"]=Sea_ice
Xt_s[,"ice2"]=Sea_ice^2
Xt_s[,"RSF"]=RSF
Xt_s = as.matrix(Xt_s)


GAM_data = data.frame(Dummy=rep(1,n_cells*t_steps),dist_land=dist_land,
                      sea_ice = Sea_ice,I_water=I_water,
                      easting=easting,
                      northing=northing,
                      RSF=RSF)
gam_setup = mgcv::gam(Dummy ~ s(dist_land, bs = "cs",k=6) + s(sea_ice, bs = "cs",k=6) +
                  s(easting, bs = "cs",k=6) + s(northing,bs="cs",k=6) + s(RSF,bs="cs",k=6),
                data = GAM_data,fit=FALSE)

S_dist_land = gam_setup$smooth[[1]]$S[[1]]
S_sea_ice = gam_setup$smooth[[2]]$S[[1]]
S_easting = gam_setup$smooth[[3]]$S[[1]]
S_northing = gam_setup$smooth[[4]]$S[[1]]
S_RSF = gam_setup$smooth[[5]]$S[[1]]

S_list = list(S_dist_land,S_sea_ice,S_easting,S_northing,S_RSF) 
S_combined = Matrix::bdiag(S_list)         # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S

#For report, used for constructing plots----
dist_land=seq(min(GAM_data$dist_land),max(GAM_data$dist_land),by = 0.1)
sea_ice = seq(min(GAM_data$sea_ice),max(GAM_data$sea_ice),by = 0.05)
easting = seq(min(GAM_data$easting),max(GAM_data$easting),by = 0.2)
northing = seq(min(GAM_data$northing),max(GAM_data$northing),by = 0.02)
RSF = seq(min(GAM_data$RSF),max(GAM_data$northing),by = 0.04)

landReport = mgcv::PredictMat(gam_setup$smooth[[1]],data = data.frame(dist_land))
seaiceReport = mgcv::PredictMat(gam_setup$smooth[[2]],data = data.frame(sea_ice))
eastingReport = mgcv::PredictMat(gam_setup$smooth[[3]],data = data.frame(easting))
northingReport = mgcv::PredictMat(gam_setup$smooth[[4]],data = data.frame(northing))
RSFReport = mgcv::PredictMat(gam_setup$smooth[[5]],data = data.frame(RSF))

designMatrixForReport = list(landReport,seaiceReport,eastingReport,northingReport,RSFReport)




# #compile TMB
TmbFile = "./inst/ChukchiPB_TMB"
compile(paste0(TmbFile,".cpp"),"-O1 -g",DLLFLAGS="")

grid_df = st_drop_geometry(grid)
loc_s = data.frame(x=grid_df[,"easting"],y=grid_df[,"northing"])
P_i = as.numeric(Area_table[,"Area_m2"]/st_area(grid[1,]))
A_s = 1-grid_df[,"land_cover"]
# 

strip_width=600
P_i_ru = as.numeric(1.2*KM_i_ru/st_area(grid[1,]))*1000000

a_photo_us = 0.012  

Output = matrix(0,3,14)  
Lambda_s = matrix(0,n_st,3)
Z_s = matrix(0,n_st,3)
colnames(Output)=c("iN_RSF","iN_tracks","it_land","it_ice","it_RSF","it_RE","converge","N","N.SE","Nall","Nall.SE","logL","K","AIC")
Output = data.frame(Output)
counter = 1

Options_vec = c("SE"=0,"Tracks"=0,"RE"=0)
Xc_s = as.matrix(Matrix::bdiag(gam_setup$X[,-1])) 

# Data combined US and RUS
Data = list( "Options_vec"=Options_vec, "N_i"=N_i,"T_i"=T_i,"C_i"=C_i,"U_i"=U_i,"P_i"=P_i,"A_s"=A_s,"S_i"=S_i-1,"Y_s"=Y_s, "G_min1"=G-1,"T_i_ru"=T_i_ru,"C_i_ru"=C_i_ru,"U_i_ru"=U_i_ru,"P_i_ru"=P_i_ru,"S_i_ru"=S_i_ru-1,"KM_i_ru"=KM_i_ru, "Xt_s"=Xt_s,"Xc_s"=Xc_s,"ho_n"=23,"ho_c"=16,"a_photo_us"=a_photo_us)
Data$S = S_combined
Data$Sdims = Sdims
Data$designMatrixForReport = Matrix::bdiag(designMatrixForReport)
Data$I_water = I_water
Data$I_us = I_us
Data$I_regehr = I_regehr
Data$I_PBSG = I_PBSG 


Data$p_trials=12
Data$p_detects=8
Data$Distances=Distances
Data$strip_width = 600
Data$g0 = 1.0
Beta_init = rep(0,ncol(Data$Xt_s))
which_water_tracks = which(names(Xt_s[1,])=="iwater")
Beta_init[1] = 0.0
Beta_init[which_water_tracks]=-50  
Alpha_init = rep(0,ncol(Data$Xc_s))

logN = log(1000)
Params = list("Beta"=Beta_init,"Alpha"=Alpha_init,"logN"=logN,"log_adjust"=0,"log_group_size_min1"=log(mean(G)-1),"log_w"=0,"tracks_eff"=0)
Params$log_sigma = 5.543
Params$logit_p = .5
Params$log_lambda = rep(rep(0,length(Sdims)))
# Random
Random= c("Alpha","log_lambda","logit_p","log_sigma")

# Fix parameters
Map = list()
Tmp = c(1,NA)
if(length(Beta_init)>2)Tmp=c(Tmp,(2:(length(Beta_init)-1)))
Map[["Beta"]] = factor(Tmp)  #for water effect on tracks; water effect on bears directly in .cpp

Output = matrix(0,3,23)  
Lambda_s = matrix(0,n_st,3)
Z_s = matrix(0,n_st,3)
colnames(Output)=c("iN_RSF","iN_tracks","it_land","it_ice","it_RSF","it_RE","converge","N","N.SE","Nall","Nall.SE","logL","K","AIC","logL_count","PBSG","Regehr","US","Russia","SE_PBSG","SE_Regehr","SE_US","SE_Russia")
Output = data.frame(Output)
counter = 1
G0 = c(.6,.8,1.0)

Reports=SDs=vector("list",3)

for(ig in 1:3){
  Data$g0 = G0[ig]
  dyn.load( dynlib(TmbFile) )
  Start_time = Sys.time()
  #setwd( "C:/Users/paul.conn/git/OkhotskST/OkhotskSeal/src/")
  Obj = MakeADFun( data=Data, parameters=Params, map=Map, silent=FALSE) #random=Random, 
  Obj$fn( Obj$par )
  
  # Run
  #Lower = -Inf
  #Upper = Inf
  start.time = Sys.time()
  Lower = -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
  Upper = 50
  #Opt = Optimize(obj=Obj,newtonsteps=2,lower=-50,upper=50,loopnum=1,bias.correct=TRUE)
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))         #
  #Opt[["diagnostics"]] = data.frame( "Param"=names(Obj$par), "Lower"=-Inf, "Est"=Opt$par, "Upper"=Inf, "gradient"=Obj$gr(Opt$par) )
  Report = Obj$report()
  cat(Sys.time() - start.time)
  
  SD = sdreport( Obj, par.fixed=Opt$par,bias.correct=FALSE )
  
  Output$converge[counter]=Opt$convergence
  Output$N[counter] = Report$N
  Output$Nall[counter] = Report$N_all
  Output$Nall.SE[counter] = SD$sd[2]
  Output$K[counter] = length(Opt$par)
  Output$AIC[counter] =  2*Report$jnll+2*length(Opt$par)
  Output$logL[counter] = -Report$jnll
  Output$logL_count[counter] = -sum(Report$jnll_comp[6:8])
  
  Lambda_s[,counter]=Report$Lambda_s
  Z_s[,counter] = Report$Z_s
  
  Survey_day = c(1:55)
  if(ig==2){
    US=Report$N_US_t
    Russia=Report$N_rus_t
    PBSG=Report$N_PBSG_t
    Regehr=Report$N_Regehr_t
  }
  Output$US[counter]=SD$value["N_US"]
  Output$Russia[counter]=SD$value["N_rus"]
  Output$PBSG[counter]=SD$value["N_PBSG"]
  Output$Regehr[counter]=SD$value["N_Regehr"]
  
  names(SD$sd)=names(SD$value)
  
  Output$SE_US[counter]=SD$sd["N_US"]
  Output$SE_Russia[counter]=SD$sd["N_rus"]
  Output$SE_PBSG[counter]=SD$sd["N_PBSG"]
  Output$SE_Regehr[counter]=SD$sd["N_Regehr"]
  
  Reports[[ig]]=Report
  SDs[[ig]]=SD
  
  counter=counter+1
  
}


Plot_df = data.frame(Abundance=c(US,Russia,PBSG,Regehr),
                     Area=c(rep("US",55),rep("Russia",55),rep("PBSG",55),rep("Regehr",55)),
                     Day = rep(1:55,4))
comp_plot = ggplot(Plot_df)+geom_line(aes(x=Day,y=Abundance,col=Area),size=2)+
  theme(text=element_text(size=14))

pdf('PB_area_composition_spatial.pdf')
  comp_plot
dev.off()



Report=Reports[[2]]  #base plots on g(0) = 0.8
#version of plotting function for sf
plot_N_map_sf<-function(N,Grid,leg.title="Abundance",myPalette=NULL,reverse=0,limits=NULL,title=NULL,grid_lines=FALSE){
  require(ggplot2)
  require(sf)
  require(RColorBrewer)
  if(is.null(myPalette))myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
  Grid$N=N
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),plot.title = element_text(hjust = 0.5))
  if(grid_lines==TRUE)p1=ggplot()+geom_sf(data=Grid,aes(fill=N)) + tmp.theme
  else p1=ggplot()+geom_sf(data=Grid,aes(fill=N),colour=NA) + tmp.theme
  if(!reverse)p1=p1+scale_fill_gradientn(colours=myPalette(100),name=leg.title,limits=limits)
  else p1=p1+scale_fill_gradientn(colours=rev(myPalette(100)),name=leg.title,limits=limits)
  if(!is.null(title))p1=p1+ggtitle(title)
  p1
}

Z_s = Report$Z_s
### plot ice & track estimates
library(grid)
library(gridExtra)
IceTracks=vector("list",6)
plot_t = 1
Z= Z_s[c(1:n_cells)+n_cells*(plot_t-1)] 
IceTracks[[1]]=plot_N_map_sf(N=Xt_s[,"ice"][c(1:n_cells)+n_cells*(plot_t-1)],Grid=Grid_list[[1]],leg.title="Ice",myPalette=colorRampPalette(brewer.pal(6, "Blues")),reverse=1,limits=c(0,1),title="7 Apr")
IceTracks[[4]]=plot_N_map_sf(N=Z,Grid=Grid_list[[1]],leg.title="Tracks")
plot_t = 28
Z= Z_s[c(1:n_cells)+n_cells*(plot_t-1)] 
IceTracks[[2]]=plot_N_map_sf(N=Xt_s[,"ice"][c(1:n_cells)+n_cells*(plot_t-1)],Grid=Grid_list[[1]],leg.title="Ice",myPalette=colorRampPalette(brewer.pal(6, "Blues")),reverse=1,limits=c(0,1),title="4 May")
IceTracks[[5]]=plot_N_map_sf(N=Z,Grid=Grid_list[[1]],leg.title="Tracks")

plot_t = 55
Ice = Xt_s[,"ice"]
Ice[which(Ice<0)]=0
Z= Z_s[c(1:n_cells)+n_cells*(plot_t-1)] 
IceTracks[[3]]=plot_N_map_sf(N=Ice[c(1:n_cells)+n_cells*(plot_t-1)],Grid=Grid_list[[1]],leg.title="Ice",myPalette=colorRampPalette(brewer.pal(6, "Blues")),reverse=1,limits=c(0,1),title="31 May")
IceTracks[[6]]=plot_N_map_sf(N=Z,Grid=Grid_list[[1]],leg.title="Tracks")



Abund = vector("list",3)
plot_t = 1
Lambda_s=Report$Lambda_s*Report$group_size
N = Lambda_s[c(1:n_cells)+n_cells*(plot_t-1)] 
Abund[[1]]=plot_N_map_sf(N=N,Grid=Grid_list[[1]],leg.title=expression(hat(N)["s,t"]),limits=c(0,10))
Abund[[1]]=Abund[[1]]+scale_fill_viridis_c(limits = c(0,20), breaks = c(0, 5, 10, 15,20),values=c(0,.05,.1,.4,1))
plot_t = 28
N = Lambda_s[c(1:n_cells)+n_cells*(plot_t-1)] 
Abund[[2]]=plot_N_map_sf(N=N,Grid=Grid_list[[1]],leg.title=expression(hat(N)["s,t"]),limits=c(0,10))
Abund[[2]]=Abund[[2]]+scale_fill_viridis_c(limits = c(0,20), breaks = c(0, 5, 10, 15,20),values=c(0,.05,.1,.4,1))

plot_t = 55
N = Lambda_s[c(1:n_cells)+n_cells*(plot_t-1)] 
Abund[[3]]=plot_N_map_sf(N=N,Grid=Grid_list[[1]],leg.title=expression(hat(N)["s,t"]),limits=c(0,10))
Abund[[3]]=Abund[[3]]+scale_fill_viridis_c(limits = c(0,20), breaks = c(0, 5, 10, 15,20),values=c(0,.05,.1,.4,1))

pdf("IceTracksAbund_spatial.pdf")
grid.arrange(arrangeGrob(IceTracks[[1]],IceTracks[[2]],IceTracks[[3]],IceTracks[[4]],IceTracks[[5]],IceTracks[[6]],Abund[[1]],Abund[[2]],Abund[[3]],nrow=3))
dev.off()


#smooth effects
SD=SDs[[2]]
muSpline = SD$value[names(SD$value)=="splineForReport"]
sdSpline<-SD$sd[names(SD$value)=="splineForReport"]

n_entries = length(muSpline)
Plot.df = data.frame(Covariate=rep(c(rep("dist_land",length(dist_land)),rep("sea_ice",length(sea_ice)),
                                     rep("easting",length(easting)),
                                     rep("northing",length(northing)), rep("RSF",length(RSF))),2),
                     x = rep(c(dist_land,sea_ice,easting,rev(northing),RSF),2),
                     Smooth = muSpline, Lower = muSpline-1.96*sdSpline, Upper = muSpline+1.96*sdSpline)

bplot = ggplot()+geom_line(data=Plot.df,aes(x=x,y=Smooth))+
  geom_ribbon(data=Plot.df,aes(x=x,ymin=Lower,ymax=Upper),alpha=0.4)+
  facet_wrap(~Covariate,scales="free")+xlab('')
pdf('pb_spatial_smooths.pdf')
bplot
dev.off()



#use 'jittered' goodness-of-fit instead
set.seed(11111)
Pval = rep(0,6)
Residuals = rep(0,n_i*3+n_i_ru*2)


Resids = pbinom(Data$T_i-1,Data$N_i,Report$Ptrack_i)+runif(n_i)*dbinom(Data$T_i,Data$N_i,Report$Ptrack_i)
Resids[Resids>1]=0.999
hist(Resids) 
Resid_binned = hist(Resids)$counts
Expected = n_i / 10
Xsq = sum((Resid_binned-Expected)^2/Expected)
Pval[1] = 1-pchisq(Xsq,df=9)   #chi square p-value for each bin
Residuals[1:n_i]=Resids

Resids = ppois(Data$T_i_ru-1,Report$Lambda_ru_i)+runif(n_i_ru)*dpois(Data$T_i_ru,Report$Lambda_ru_i)
Resids[Resids>1]=0.999
hist(Resids)  # yuck!  Def zero inflation *and* pos outliers.  Well, we just want broadscale trends in track index anyway so not sure this matters
Resid_binned = hist(Resids)$counts
Expected = n_i_ru / 10
Xsq = sum((Resid_binned-Expected)^2/Expected)
Pval[2] = 1-pchisq(Xsq,df=9)   #chi square p-value for each bin
Residuals[(n_i+1):(n_i+n_i_ru)]=Resids

#look at tracks another way
hist(T_i_ru/KM_i_ru,breaks=100)  #one extreme outlier but not many KM for that one
hist(T_i/N_i)   
#we've got to have inflated sample sizes due to bears leaving more than 1 track and them staying around for awhile
#hard to interpret; but good thing we dont' care about variance of the track index

#residuals for bear counts
Resids = ppois(Data$C_i-1,Report$E_count_on)+runif(n_i)*dpois(Data$C_i,Report$E_count_on)
Resids[Resids>1]=0.999
hist(Resids) 
Resid_binned = hist(Resids)$counts
Expected = n_i / 10
Xsq = sum((Resid_binned-Expected)^2/Expected)
Pval[3] = 1-pchisq(Xsq,df=9)   #chi square p-value for each bin
Residuals[(n_i+n_i_ru+1):(2*n_i+n_i_ru)]=Resids

Resids = ppois(Data$U_i-1,Report$E_count_off)+runif(n_i)*dpois(Data$U_i,Report$E_count_off)
Resids[Resids>1]=0.999
hist(Resids) 
Resid_binned = hist(Resids)$counts
Expected = n_i / 10
Xsq = sum((Resid_binned-Expected)^2/Expected)
Pval[4] = 1-pchisq(Xsq,df=9)   #chi square p-value for each bin
Residuals[(2*n_i+n_i_ru+1):(3*n_i+n_i_ru)]=Resids

Resids = ppois(Data$C_i_ru-1,Report$E_count_on_ru)+runif(n_i_ru)*dpois(Data$C_i_ru,Report$E_count_on_ru)
Resids[Resids>1]=0.999
hist(Resids) 
Resid_binned = hist(Resids)$counts
Expected = n_i_ru / 10
Xsq = sum((Resid_binned-Expected)^2/Expected)
Pval[5] = 1-pchisq(Xsq,df=9)   #chi square p-value for each bin
Residuals[(3*n_i+n_i_ru+1):(3*n_i+2*n_i_ru)]=Resids




Resids.df = data.frame(Residual = Residuals)
Resids.df$Type = c(rep("Tracks",n_i+n_i_ru),rep("Bears-on",n_i),rep("Bears-off",n_i),
                   rep("Bears-on",n_i_ru))
Resids.df$Country = c(rep("U.S.",n_i),rep("Russia",n_i_ru),rep("U.S.",n_i*2),rep("Russia",n_i_ru*1))
Resids.df$Type = factor(Resids.df$Type)
Resids.df$Country = factor(Resids.df$Country)

Labels = rep('',5)
for(i in 1:5){
  Labels[i] = paste0('p=',format(Pval[i],nsmall=2,digits=2))
}
Labels[1]="p=0.00"

Labels.df = data.frame(Label = Labels,Type=as.factor(c("Tracks","Tracks","Bears-on","Bears-off","Bears-on")),
                       Country=as.factor(c("U.S.","Russia","U.S.","U.S.","Russia")),
                       x=0.9,y=200)

myplot=ggplot()+geom_histogram(data=Resids.df,aes(x=Residual),breaks=c(0:10)/10)+facet_grid(Type~Country)+ylab('Frequency')+xlab('Jittered CDF')+
  geom_text(data=Labels.df,aes(label=Labels,x=x,y=y))
pdf('GOF_bears_tracks.pdf')
myplot
dev.off()

#average # dependents / female
tabulate(Data$G_min1)  #nine had a single dependent, 6 had 2 dependents - average = 1.4


#95% log-based credible intervals
#total abundance
for(ig in 1:3){
  cat(paste("total: ",Reports[[ig]]$N_all,'\n'))
  C = exp(1.96 * sqrt(log(1+(SDs[[ig]]$sd[2]/Reports[[ig]]$N_all)^2)))
  CI_lo<- Reports[[ig]]$N_all/C
  CI_hi<- Reports[[ig]]$N_all*C
  cat(paste("total: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))
}

#Regehr area
which_sd = which(names(SD$value)=="N_Regehr")
for(ig in 1:3){
  cat(paste("Regehr: ",SDs[[ig]]$value[which_sd],'\n'))
  C = exp(1.96 * sqrt(log(1+(SDs[[ig]]$sd[which_sd]/SDs[[ig]]$value[which_sd])^2)))
  CI_lo<- SDs[[ig]]$value[which_sd]/C
  CI_hi<- SDs[[ig]]$value[which_sd]*C
  cat(paste("Regehr: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))
}

#PBSG
which_sd = which(names(SD$value)=="N_PBSG")
for(ig in 1:3){
  cat(paste("PBSG: ",SDs[[ig]]$value[which_sd],'\n'))
  C = exp(1.96 * sqrt(log(1+(SDs[[ig]]$sd[which_sd]/SDs[[ig]]$value[which_sd])^2)))
  CI_lo<- SDs[[ig]]$value[which_sd]/C
  CI_hi<- SDs[[ig]]$value[which_sd]*C
  cat(paste("PBSG: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))
}

#Russia
which_sd = which(names(SD$value)=="N_rus")
for(ig in 1:3){
  cat(paste("Russia: ",SDs[[ig]]$value[which_sd],'\n'))
  C = exp(1.96 * sqrt(log(1+(SDs[[ig]]$sd[which_sd]/SDs[[ig]]$value[which_sd])^2)))
  CI_lo<- SDs[[ig]]$value[which_sd]/C
  CI_hi<- SDs[[ig]]$value[which_sd]*C
  cat(paste("Russia: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))
}


#US
which_sd = which(names(SD$value)=="N_US")
for(ig in 1:3){
  cat(paste("US: ",SDs[[ig]]$value[which_sd],'\n'))
  C = exp(1.96 * sqrt(log(1+(SDs[[ig]]$sd[which_sd]/SDs[[ig]]$value[which_sd])^2)))
  CI_lo<- SDs[[ig]]$value[which_sd]/C
  CI_hi<- SDs[[ig]]$value[which_sd]*C
  cat(paste("US: (", round(CI_lo,0), "-", round(CI_hi,0), ")",  sep=""))
}