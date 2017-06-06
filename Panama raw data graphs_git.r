##Panama graphs
R
library(rethinking)

#below is codw for curly braces for figure 1
#http://stackoverflow.com/questions/23869526/plot-curly-bracket-for-axes-area
CurlyBraces <- function(x0, x1, y0, y1, pos = 1, direction = 1, depth = 1) {

    a=c(1,2,3,48,50)    # set flexion point for spline
    b=c(0,.2,.28,.7,.8) # set depth for spline flexion point

    curve = spline(a, b, n = 50, method = "natural")$y * depth

    curve = c(curve,rev(curve))

    if (pos == 1){
        a_sequence = seq(x0,x1,length=100)
        b_sequence = seq(y0,y1,length=100)  
    }
    if (pos == 2){
        b_sequence = seq(x0,x1,length=100)
        a_sequence = seq(y0,y1,length=100)      
    }

    # direction
    if(direction==1)
        a_sequence = a_sequence+curve
    if(direction==2)
        a_sequence = a_sequence-curve

    # pos
    if(pos==1)
        lines(a_sequence,b_sequence, lwd=1.5,   xpd=NA) # vertical
    if(pos==2)
        lines(b_sequence,a_sequence, lwd=1.5, xpd=NA) # horizontal

}


mono_index <- read.csv("~/Dropbox/Panama Data/mono_indexing.csv") #without HF and WZ

mono_index$mono <- mono_index$monos
ages <- subset(mono_index, select=c(yob,mono) )

setwd("~/Dropbox/Panama Data/Panama Shared/github stuff")
d <- read.csv("~/Dropbox/Panama Data/Panama Shared/github stuff/panama_data_14days.csv" , header=TRUE)

d$date_index <- as.integer(as.factor(d$date_index))

d <- merge(d,ages,by="mono")
d$age <- 2015-d$yob 
d$adult <- ifelse(d$age<7 , 0 , 1)
#datalists
d$nobs500.1 <- ifelse(d$nobs500==0 , 1 , d$nobs500)
d$season <- ifelse(d$date_index > 48 , 1 ,0)
mono_index$age <- 2015-mono_index$yob


tech_index <- as.vector(sort(unique(d$tech)))
tech_index2 <- tech_index
tech_index2[8] <- "OPEN"
tech_index2[9] <- "FAIL"

yobs <- mono_index[,3:5]
yobs$mono_index <- yobs$X1.21
d <- merge(yobs,d,by="mono_index")
d$yob_orderhi <- d$yob_order + 0.3
#d$monoindexlo <- d$mono_index - 0.1
col_index= c("darkgreen","blue","red","gold","grey","violet","orange")
col_index2= c("darkgreen","blue","red","gold","grey","violet","orange","black","black")
col_index2bg= c("darkgreen","blue","red","gold","grey","violet","orange","black","black")


#next get graps for individuals foraging and individuals observing and tidy up larger one

##figure s4 code below
cairo_pdf("Heat_map_techs_wide.pdf", height=7,width=20)

par(mar=c(3.5,3.25,0.5,3) , oma=c(0,0,0,0) )
plot(yob_order ~ timestep , data=subset(d,fruit_index<1500 & tech_index==1 & open==1) , pch=""  , ylim= c(-1,23) , xlim=c(-7,nrow(d) + 7) , xlab=NA , yaxt="n", xaxt="n"  , ylab=NA , xaxs="i")

for (i in 1:max(d$tech_index)){
	points(yob_order~timestep, data=subset(d,fruit_index<1500 & tech_index==i & open==1) , pch=21 , cex=1.3 , col=col_index[i] , bg=col.alpha(col_index[i] , 0.25) )
		points(yob_order~timestep, data=subset(d,fruit_index<1500 & tech_index==i & open==0) , pch=4 , cex=1.3 , col=col_index[i] )
}
axis(4, at=mono_index$yob_order , labels=mono_index$monos, las=2 , cex.axis=1 , tick=FALSE , lwd.ticks=0 , hadj=0.6)
axis(2, at=mono_index$yob_order , labels=mono_index$yob, las=2 , cex.axis=1 , tick=FALSE , lwd.ticks=0 , hadj=0.5)
axis(1, at=seq(from=0,to=1440,by=100) , labels=TRUE , las=1 , tick=TRUE , lwd.ticks=1  , tck=-0.02)
axis(1, at=seq(from=0,to=1440,by=10) , labels=FALSE, cex.axis=0.5 , tick=TRUE , lwd.ticks=1 , tck=-0.01)
mtext("Individual ID" , side=4 , line=1.9, outer=FALSE , cex=1.8)
mtext("Year of Birth" , side=2 , line=1.9, outer=FALSE , cex=1.8)
mtext("Panama fruit processing event time (N=1441 fruits)" , side=1 , line=2.3, outer=FALSE , cex=1.8)

legend("bottom", inset=c(0.1,0), pch=c(15,15,15,15,15,15,15,21,4) , col=col_index2, tech_index2, horiz=TRUE , cex=1.4 , bty = "n" ,  bg=col.alpha(col_index2[i] , 0.25) )

dev.off()


df <- table(d$tech,d$date_index,d$mono_index)

sort(unique(d$date))
cairo_pdf("fig1_techs_id_days.pdf", height=4,width=10 )
par(mar=c(0,0,0,0) , oma=c(2.75,2.75,1.25,2.75))
plot(-2,-2, xlim=c(1,75) , ylim=c(0,23) , xlab="", ylab="" ,yaxt='n' ,xaxt='n')
pointlist=c(15,16,17,18,1,4,20)
for(i in 1:7){
	for(j in 1:75){
		for(k in 1:23){
			if(df[i,j,k] > 0){
				points(j,max(d$yob_order[d$mono_index==k]) , col=col_index[which(df[,j,k] == max(df[,j,k]), arr.ind = TRUE)] , pch=pointlist[which(df[,j,k] == max(df[,j,k]), arr.ind = TRUE)] , cex=1)
				#points(j,max(d$yob_order[d$mono_index==k]) , col=col_index[which(df[,j,k] == max(df[,j,k]), arr.ind = TRUE)] , pch=pointlist[which(df[,j,k] == max(df[,j,k]), arr.ind = TRUE)] , cex=0.9)
			}
		}
	}
}

axis(4, at=mono_index$yob_order , labels=mono_index$monos, las=2 , cex.axis=0.75 , tick=FALSE , lwd.ticks=0 , hadj=0.9)
axis(2, at=mono_index$yob_order , labels=mono_index$yob, las=2 , cex.axis=0.75 , tick=FALSE , lwd.ticks=0 , hadj=0.3)
axis(1, at=seq(from=0,to=70,by=10) , labels=TRUE, las=1 , cex.axis=1 , tick=TRUE , lwd.ticks=1 , padj=-1,tck=-0.04)
axis(1, at=seq(from=0,to=75,by=5) , labels=FALSE , las=1 , tick=TRUE , lwd.ticks=1  , tck=-0.02)
axis(1, at=seq(from=0,to=75,by=1) , labels=FALSE, cex.axis=0.25 , tick=TRUE , lwd.ticks=1 , tck=-0.01)
#axis(3, at=c(1,4,11,33,49,54,67) , labels=c("01/13","12/13","01/14","02/14","12/14","01/15","02/15"), cex.axis=0.55 , tick=TRUE , lwd.ticks=1 , tck=-0.02 , padj=2 )
axis(3, at=c(mean(c(1,3)), mean(c(4,48)) , mean(c(49,75)) ) , labels=c("Jan 2013" , "Dec 2013-Feb 2014" , "Dec 2014-Feb 2015"), cex.axis=0.7 , tick=FALSE , lwd.ticks=1 , tck=-0.02 , padj=1.2 )
CurlyBraces(x0=1,  x1=3,  y0=24, y1=24, pos = 2, direction = 1, depth=0.75)
CurlyBraces(x0=4,  x1=48,  y0=24, y1=24, pos = 2, direction = 1, depth=0.75)
CurlyBraces(x0=49,  x1=75,  y0=24, y1=24, pos = 2, direction = 1, depth=0.75)

mtext("Experimental Days (N=75)" , side=1 , line=1.5, outer=FALSE , cex=1.5)
mtext("Individual ID" , side=4 , line=1.5, outer=FALSE , cex=1.3)
mtext("Year of Birth" , side=2 , line=1.5, outer=FALSE , cex=1.3)
legend("bottom", pch=pointlist , col=col_index, tech_index, horiz=TRUE , cex=.85 , bty = "n" ,  , inset=-.025)

dev.off()

#df[fruit,day,monkey]

df[,,16]
max(df[,,16])
points(1:75, df[3,,16] )

#cairo_pdf("RCI_shade.pdf", height=7,width=21)

#fig s5 code below
cairo_pdf("RCI_star.pdf", height=6,width=18)
par(mar=c(0,0,0,0) , oma=c(3.5,3.25,0.5,3))
plot(yob_order ~ timestep , data=subset(d,fruit_index<1500 & tech_index==1 & open==1) ,  pch=""  , ylim= c(-2,23) , xlim=c(-7,nrow(d) + 7) , xlab=NA , yaxt="n" , ylab=NA , xaxs="i" )
for (i in 1:nrow(d)){
	for (j in 1:23){
		if ( grepl(mono_index$monos[j],d[i,"RCI"])==TRUE  )
			#points( mono_index$yob_order[j] ~ d$timestep[i]  , pch=15 , cex=2 , col=col.alpha("black" , .2) )
			points( mono_index$yob_order[j] ~ d$timestep[i]  , pch="*" , cex=1.2 , col=1 )
		}
}
for (i in 1:nrow(d)){
		if ( grepl("WZ",d[i,"RCI"])==TRUE  )
			#points( 0 ~ d$timestep[i]  , pch=15 , cex=2 , col=col.alpha("black" , .25) )
			points( 0 ~ d$timestep[i]  , pch="*" , cex=1.2 , col=1 )

}
for (i in 1:nrow(d)){
		if ( grepl("HF",d[i,"RCI"])==TRUE  )
			#points( -1 ~ d$timestep[i]  , pch=15 , cex=2 , col=col.alpha("black" , .25) )
			points( -1 ~ d$timestep[i]  , pch="*" , cex=1.2 , col=1 )

}

axis(4, at=mono_index$yob_order , labels=mono_index$monos, las=2 , cex.axis=1 , tick=FALSE , lwd.ticks=0 , hadj=0.6)
axis(2, at=mono_index$yob_order , labels=mono_index$yob, las=2 , cex.axis=1 , tick=FALSE , lwd.ticks=0 , hadj=0.5)
axis(4, at=c(0,-1) , labels=c("WZ","HF"), las=2 , cex.axis=1 , tick=FALSE , lwd.ticks=0 , hadj=0.6)
axis(2, at=c(0,-1) , labels=c(2014,2014), las=2 , cex.axis=1 , tick=FALSE , lwd.ticks=0 , hadj=0.5)
axis(1, at=seq(from=0,to=1440,by=100) , labels=TRUE , las=1 , tick=TRUE , lwd.ticks=1  , tck=-0.02)
axis(1, at=seq(from=0,to=1440,by=10) , labels=FALSE, cex.axis=0.5 , tick=TRUE , lwd.ticks=1 , tck=-0.01)
mtext("Individual ID" , side=4 , line=1.9, outer=FALSE , cex=1.8)
mtext("Year of Birth" , side=2 , line=1.9, outer=FALSE , cex=1.8)
mtext("Panama fruit processing event time (N=1441 fruits)" , side=1 , line=2.3, outer=FALSE , cex=1.8)

dev.off()

d[order(as.Date(d$V3, format="%d/%m/%Y")),]



########plot raw data old graphgs used in conference talks that are not in publication
d$date_index <- as.integer(as.factor(d$date_index))
d$num_open <- 0
d$num_processed <- 0
for (i in 1:75){
    #plot( sum(d$open[d$date_index==i]) / length(unique(d$fruit_index[d$date_index==i]))  ~  i, data=d)
    #points( (d$num_open/d$num_processed) ~ d$date_index, data=d)
    d$num_open <- ifelse( d$date_index==i, sum(d$open[d$date_index==i]),d$num_open )
    d$num_processed <- ifelse( d$date_index==i, length(unique(d$fruit_index[d$date_index==i])),d$num_processed )
}
plot( ( d$num_open/d$num_processed) ~ d$date_index, col="black", xlab="Experimental Days (N=75)", ylab="Proportion of Fruits Opened" , pch=19 , cex.lab=1.5)
plot( d$s1~d$date_index , col="white")

mat <- matrix(, nrow = 75, ncol = 8)
for (i in 1:75){
    for( j in 1:7) {
        mat[i,j] <- length(unique(d$fruit_index[d$date_index==i & d$tech_index==j]))/length(unique(d$fruit_index[d$date_index==i]))
        
    }
}

cairo_pdf("raw_data_splines.pdf", height=7,width=10)

par(oma=c(1,1,0,0))
matlo <- mat
for (i in 1:7){
	matlo[,i] <- ifelse(matlo[,i]==0, 0.001, matlo[,i])
	matlo[,i] <- ifelse(matlo[,i]==1, 0.999, matlo[,i])

}
tech_id <- as.vector(sort(unique(d$tech)))
mat[,8] <- c(1:75)
#col_ind <- c("black","red","gold","orange", "cyan" , "blue" , "magenta" )
plot(mat[,i]~mat[,8], ylim=c(0,1) , col="white" , xlab="Experimental Days (N=75)" , ylab="daily proportion of observed techniques" , cex.lab=1.5)
#legend("top", inset=0.01, c(tech_id) , fill=col_index, border=col_index, horiz=TRUE,cex=0.75,bty = "n")

for (i in 1:7){
    points(matlo[,i]~matlo[,8], col=col.alpha(col_index[i], alpha=0.3) , pch=19 )
    ss <- smooth.spline( x=matlo[,8],y=logit(matlo[,i]) , spar=0.95 ) 
    lines( ss$x, logistic(ss$y), col=col_index[i] , lw=4)

}

dev.off()
library(rms)
library(Hmisc)
rcspline.plot(mat[,8],mat[,3] , model="logistic", ylim=c(-4,4), lty=1)

