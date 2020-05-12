############################
#####distance sampling#####

data(Chukchi_PB_data)
###plot data  
par(mfrow=c(3,1))
hist(pb_dist$distance[pb_dist$photo==1], breaks=c(seq(0,1500,by=50)), xlab="Distance, m", ylab="Frequency", main="Photo") # photo (including photo+vis)
hist(pb_dist$distance[pb_dist$vis==1], breaks=c(seq(0,1500,by=50)),  xlab="Distance, m", ylab="Frequency", main="Visual") # visual (including vis+photo)
hist(pb_dist$distance, breaks=c(seq(0,1500,by=50)), xlab="Distance, m", ylab="Frequency", main="All data") #all_data  (no duplicates)

pb.fixed <- Distance:::checkdata(pb_dist)
pb.fixed$platform <-as.factor(pb.fixed$platform)

result_hn <- ddf(dsmodel=~mcds(key="hn", formula=~1), data=pb.fixed$data, method="ds",     #AIC=  583.7962
                 meta.data=list(width=600, point=FALSE))
result_hr <- ddf(dsmodel=~mcds(key="hr", formula=~1), data=pb.fixed$data, method="ds",     #AIC=  586.0528 
                 meta.data=list(width=600, point=FALSE))
result_unif <- ddf(dsmodel=~mcds(key="unif", formula=~1), data=pb.fixed$data, method="ds",     #AIC=  601.3114 
                   meta.data=list(width=600, point=FALSE))
result_gamm <- ddf(dsmodel=~mcds(key="gamma", formula=~1), data=pb.fixed$data, method="ds",     #AIC=  643.808 
                   meta.data=list(width=600, point=FALSE))

par(mfrow = c(2, 2))
plot.title <- "Combined distance data HN"
plot(result_hn, lwd = 2, lty = 1, pch = 1, cex = 1, nc = 10, main = plot.title)
#ddf.gof.results <- ddf.gof(result_hn, lwd = 2, lty = 1, pch = ".", cex = 0.5, col = c(1,2))
pdf('distance_hn.pdf')
plot.title <- "Combined distance data HN"
plot(result_hn, lwd = 2, lty = 1, pch = 1, cex = 1, nc = 10, main = plot.title)
dev.off()

plot.title <- "Combined distance data HR"
plot(result_hr, lwd = 2, lty = 1, pch = 1, cex = 1, nc = 10, main = plot.title)
#ddf.gof.results <- ddf.gof(result_hr, lwd = 2, lty = 1, pch = ".", cex = 0.5, col = c(1,2))

plot.title <- "Combined distance data Uniform"
plot(result_unif, lwd = 2, lty = 1, pch = 1, cex = 1, nc = 10, main = plot.title)
#ddf.gof.results <- ddf.gof(result_unif, lwd = 2, lty = 1, pch = ".", cex = 0.5, col = c(1,2))

plot.title <- "Combined distance data GAMMA"
plot(result_gamm, lwd = 2, lty = 1, pch = 1, cex = 1, nc = 10, main = plot.title)
#ddf.gof.results <- ddf.gof(result_gamm, lwd = 2, lty = 1, pch = ".", cex = 0.5, col = c(1,2))


## see if we need adjustment (tried for HN and HR only)
result_hn_cos <- ds(pb_dist, key = c("hn"), adjustment = c("cos"),
                    convert.units = 0.001, truncation=list(left=0,right=600), dht.group = FALSE)
result_hn_herm <- ds(pb_dist, key = c("hn"), adjustment = c("herm"),
                     convert.units = 0.001, truncation=list(left=0,right=600), dht.group = FALSE)
result_hn_poly <- ds(pb_dist, key = c("hn"), adjustment = c("poly"),
                     convert.units = 0.001, truncation=list(left=0,right=600), dht.group = FALSE)


result_hr_cos <- ds(pb_dist, key = c("hr"), adjustment = c("cos"),
                    convert.units = 0.001, truncation=list(left=0,right=600), dht.group = FALSE)
result_hr_herm <- ds(pb_dist, key = c("hr"), adjustment = c("herm"),
                     convert.units = 0.001, truncation=list(left=0,right=600), dht.group = FALSE)
result_hr_poly <- ds(pb_dist, key = c("hr"), adjustment = c("poly"),
                     convert.units = 0.001, truncation=list(left=0,right=600), dht.group = FALSE)


##see if we can treat photo/vis as different platforms (but they are dependant, so probably not, just for the sake of plotting)
result_hn_platf <- ddf(dsmodel=~mcds(key="hn", formula=~platform), data=pb.fixed$data, method="ds",     #AIC=  584.6212   ##added platform (vis or photo as a factor)
                       meta.data=list(width=600, point=FALSE))  ##AIC 584.6212 
par(mfrow=c(1,1))
plot.title <- "Combined distance data HN with platform as a factor"
plot(result_hn_platf, lwd = 2, lty = 1, pch = 1, cex = 1, nc = 10, main = plot.title)
#ddf.gof.results <- ddf.gof(result_hn_platf, lwd = 2, lty = 1, pch = ".", cex = 0.5, col = c(1,2))


##plot data (photoand visual and visual only) side by side

par(mfrow=c(1,1))
l <- list(pb_dist$distance[pb_dist$photo==1],pb_dist$distance[pb_dist$photo==0])
multhist(l, breaks=c(seq(0,1500,by=50)), xlab="Distance, m", ylab="Frequency")
legend("topright", legend=c("Photo/Photo and Visual", "Visual only"), col=c("darkgrey", "lightgrey"), cex = 1.2, pch=15, bty="n")

pdf('distance_both.pdf')
multhist(l, breaks=c(seq(0,1500,by=50)), xlab="Distance, m", ylab="Frequency")
legend("topright", legend=c("Photo/Photo and Visual", "Visual only"), col=c("darkgrey", "lightgrey"), cex = 1.2, pch=15, bty="n")
dev.off()

###fit detection function separately for visual only and for photo+vis (no duplicates)

pb_photo<- subset(pb.fixed$data, pb.fixed$data$photo==1)       ## detected by photo only or by photo+visually 
pb_vis<- subset(pb.fixed$data, pb.fixed$data$photo==0)      ##detected only visually

result_hn_photo <- ddf(dsmodel=~mcds(key="hn", formula=~1), data=pb_photo, method="ds",   
                       meta.data=list(width=600, point=FALSE))  ##AIC 368.899 
result_hn_vis <- ddf(dsmodel=~mcds(key="hn", formula=~1), data=pb_vis, method="ds",     
                     meta.data=list(width=600, point=FALSE))  #AIC   215.7222 




plot.title <- "Distance data HN for bears detected only visually and on photo"
plot(result_hn_vis, lwd = 2, lty = 1, pch = 1, cex = 1,nc=10, main = plot.title)


#plot.title <- "Distance data HN for photographed bears (photo/photo+vis)"
plot(result_hn_photo, lwd = 2, lty = 1, pch = 1, cex = 1, nc=10, main = plot.title, add=T)

## if add=T, plots one function on top of another on the same scale...