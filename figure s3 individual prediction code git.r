#####BRENDAN GRAPHS OF PREDICTIONS##########

####Richard multiplicative model plots for individuals
library(rethinking)
#setwd("~/Dropbox/Panama Data/Panama Shared")
load("Rstan_mods14dayCHisBA.RData")

#lassuming global model is fit

post <- extract(fit_global_age)
d <- read.csv("~/panama_data_14days.csv" , header=TRUE) #data file in github

Softmax <- function(x){
exp(x)/sum(exp(x))
} #softmax function to simplify code

##use same data list that was used when fitting stan model
d_global <- list(
    N = nrow(d),                                                                        #length of dataset
    J = length( unique(d$mono_index) ),                                                 #number of individuals
    K=max(d$tech_index),                                                                #number of processing techniques
    tech = d$tech_index,                                                                #technique index
    y = cbind( d$y1 , d$y2 , d$y3 , d$y4 , d$y5 , d$y6 , d$y7 ),                        #individual processing times for all techniques at each bout N (individual payoff)
    s = cbind(d$s1 , d$s2 , d$s3 , d$s4 ,d$s5 ,d$s6 , d$s7 ),                           #observed counts of all K techniques to individual J (frequency-dependence)
    ps = cbind(d$ps1 , d$ps2 , d$ps3 , d$ps4 ,d$ps5 ,d$ps6 , d$ps7 ),                   #observed mean payoffs of all K techniques to individual J (payoff bias)
    ks = cbind(d$ks1 , d$ks2 , d$ks3 , d$ks4 ,d$ks5 ,d$ks6 , d$ks7 ),                   #observed matrilineal kin cues of all K techniques to individual J (matrilineal kin-bias)
    cohos = cbind(d$cos1 , d$cos2 , d$cos3 , d$cos4 ,d$cos5 ,d$cos6 , d$cos7 ),         #observed cohort cues of all K techniques to individual J (age-cohort/similarity bias)
    yobs = cbind(d$Yobs1 , d$Yobs2 , d$Yobs3 , d$Yobs4 ,d$Yobs5 ,d$Yobs6 , d$Yobs7 ),   #observed age cues of all K techniques to individual J (age-bias)
    press = cbind(d$prs1 , d$prs2 , d$prs3 , d$prs4 ,d$prs5 ,d$prs6 , d$prs7 ),         #observed rank cues of all K techniques to individual J (age-bias)
    bout = d$forg_bout,#bout is forg index here                                         #processing bout unique to individual J
    id = as.integer(as.factor(d$mono_index)),                                           #individual ID
    N_effects=8,                                                                        #number of parameters to estimates
    age=d$age.c                                                                         #centered age of individuals
)

#scale payoffs/cues by dividing by max value
d1 <- d_global
d1$yobs <- d1$yobs / max(d1$yobs)
d1$cohos <- d1$cohos / max(d1$cohos)
d1$ks <- d1$ks / max(d1$ks)
d1$ps <- d1$ps / max(d1$ps)
d1$press <- d1$press / max(d1$press)

d1$fruit_index <- d$fruit_index
d1$date_index <- d$date_index
Preds = array(0,dim=c(nrow(d),7,23)) #predictions for all individuals, all techniques, across all timesteps
Preds2 = array(0,dim=c(nrow(d),7)) ##predictions for all individuals, all techniques, at times when they foraged

lambda = mean(post$lambda)
AC=PrS=PrA=lin_mod=s_temp=rep(0,7) #stroage slots for calculating predictions

for ( i in 1:max(d1$N) ) {
    #if ( d1$bout[i]==1 ) {
        #// calculate new individual's parameter values
        phi = mean(logistic(post$mu[,1] + post$a_id[,d1$id[i],1] + post$b_age[,1]*ages$age.c[d1$id[i]] ))
        fconf = mean(exp( post$mu[,3] + post$a_id[,d1$id[i],3] ))
        fpay = mean( post$mu[,4] + post$a_id[,d1$id[i],4] )
        gamma = mean(logistic( post$mu[,2] + post$a_id[,d1$id[i],2] + post$b_age[,2]*ages$age.c[d1$id[i]]))
        ##PICK UP HERE
        fkin = mean(( post$mu[,5] + post$a_id[,d1$id[i],5] ))
        fpres = mean(( post$mu[,6] + post$a_id[,d1$id[i],6] ))
        fcoho = mean(( post$mu[,7] + post$a_id[,d1$id[i],7] ))
        fyob = mean(( post$mu[,8] + post$a_id[,d1$id[i],8] ))
        #}
        #//update attractions
        for ( j in 1:max(d1$tech) ) {
            if ( d1$bout[i] > 1 ) {
                AC[j] = (1-phi)*AC[j] + phi*d1$y[i-1,j]
            } else {
                AC[j] = 0;
            }
        }#//j

        for (j in 1:max(d1$tech)){PrA[j] = exp(lambda*AC[j])/sum(exp(lambda*AC))}
        
        #//conformity aspect below
        if ( d1$bout[i] > 1 ) {
            if (sum( d1$s[i,] ) > 0 ) {

                #// compute non-frequency cues as log-linear model
                for ( j in 2:max(d1$tech) ) {
                    lin_mod[j] = exp( fpay*d1$ps[i,j] + fkin*d1$ks[i,j] + fpres*d1$press[i,j] + fcoho*d1$cohos[i,j] + fyob*d1$yobs[i,j])
                }
                lin_mod[1] = 1; #// aliased outcome

                #// compute frequency cue
                for ( j in 1:max(d1$tech) ){ s_temp[j] = d1$s[i,j]^fconf }
                for ( j in 1:max(d1$tech) ){ lin_mod[j] = lin_mod[j] * s_temp[j] }

                for (j in 1:7){PrS[j] = lin_mod[j]/sum(lin_mod)}
                
                for(j in 1:7){ Preds[d1$fruit_index[i],j,d1$id[i]] = (1-gamma)*PrA[j] + gamma*PrS[j] 
                               Preds2[i,j] = (1-gamma)*PrA[j] + gamma*PrS[j]
                }

            } else {
                for(j in 1:7){ Preds[d1$fruit_index[i],j,d1$id[i]]= PrA[j] 
                               Preds2[i,j] = PrA[j]}
            }
        } else {
            for(j in 1:7){ Preds[d1$fruit_index[i],j,d1$id[i]]= PrA[j]
                           Preds2[i,j] = PrA[j] }
         }
     }#//i  





col_index= c("darkgreen","blue","red","gold","grey","violet","orange")

d$date_index <- as.integer(as.factor(d$date_index))
matF <- matrix(, nrow = 75, ncol = 8)
matP <- matrix(, nrow = 75, ncol = 8)
mat <- matrix(, nrow = 75, ncol = 8)

for (i in 1:75){
    for( j in 1:7) {
        matF[i,j] <- length(unique(d$fruit_index[d$date_index==i & d$tech_index==j]))/length(unique(d$fruit_index[d$date_index==i]))
        matP[i,1] <- mean(d$ps1[d$date_index==i])
        matP[i,2] <- mean(d$ps2[d$date_index==i])
        matP[i,3] <- mean(d$ps3[d$date_index==i])
        matP[i,4] <- mean(d$ps4[d$date_index==i])
        matP[i,5] <- mean(d$ps5[d$date_index==i])
        matP[i,6] <- mean(d$ps6[d$date_index==i])
        matP[i,7] <- mean(d$ps7[d$date_index==i])
    }
}

tech_id <- as.vector(sort(unique(d$tech)))
mat[,8] <- c(1:75)

    for(k in 1:23){
        for (i in 1:1440){
        for (j in 1:7){
            Preds[i+1,j,k] <- ifelse(Preds[i+1,j,k]==0 & Preds[i,j,k]!=0,Preds[i,j,k],Preds[i+1,j,k] )
        }
    }
}


##plots where probability is constant between fruits
#Preds[d1$fruit_index[i],j,d1$id[i]]
mono_index_all <- read.csv("~/Dropbox/Panama Data/mono_indexing_all.csv") #with HF and WZ
mono_index_all$mono <- mono_index_all$monos
mono_index <- mono_index_all

##below plot was not in paper but plots time series of behavioral change against all foraging behaviors, useful for analyst

pdf("individual_panama_predictions_fruit_pred_all.pdf",width=11,height=8.5) 

par(mfrow=c(4, 1) , oma = c(5,3,1,0) , mar=c(3,3,1.5,0))
par(cex = 0.6)
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

for(k in 1:23){
    plot(d1$date_index~d1$fruit_index, ylim=c(0,1.2) , col="white" , xlab="" , ylab="" , cex.lab=.8 ,main=mono_index$mono[k], axes=FALSE, cex.lab=1.1)
           # legend("top", inset=0.01, c(tech_id) , fill=col_index, horiz=TRUE,cex=0.6,bty = "n")
    axis(1, at = seq(from=0 , to=1450, by = 100) , tck=-0.03)
    axis(2, at = seq(from=0 , to=1 , by = 0.2) , tck=-0.03 )
    axis(1, at = seq(from=0 , to=1450, by = 25) , tck=-0.015 , labels=FALSE )
    axis(2, at = seq(from=0 , to=1 , by = 0.1) , tck=-0.015 , labels=FALSE)
    for (j in 1:7){    
        for (i in 1:1441){
           points(Preds[i,j,k]~i, col=col.alpha(col_index[j], alpha=0.99) , pch=20 , cex=0.5 )
        }
        lls <- d$timestep[d$fruit_index<1500 & d$tech_index==j & d$open==1 & d$mono_index==k]
        llf <- d$timestep[d$fruit_index<1500 & d$tech_index==j & d$open==0 & d$mono_index==k]
        llo <-  d$timestep[grepl(mono_index$monos[k],d[,"RCI"])==TRUE & d$tech_index==j]

            #points( mono_index$yob_order[j] ~ d$timestep[i]  , pch=15 , cex=2 , col=col.alpha("black" , .2) )
        points(llo,rep(1.15,length(llo))  , pch="*" , cex=1.2 , col=col_index[j] )
        points(lls,rep(1.05,length(lls)), pch=23 , cex=0.75 , col=col_index[j] , bg=col.alpha(col_index[j] , 0.5) )
        points(llf,rep(1.1,length(llf)), pch=4 , cex=0.75 , col=col_index[j] , bg=col.alpha(col_index[j] , 0.5) )    
    
    }
}

mtext("fruit processing event time", side = 1, outer = TRUE, cex = 1.5, line = .1)
mtext("probability of choosing technique", side = 2, outer = TRUE, cex = 1.5, line = .1)

legend(x=200 , y=-.5, inset= -.1, tech_id , fill=col_index, border=col_index ,horiz=TRUE,cex=1, bty = "n", xpd=NA)
legend(x=500 , y=-.65, inset= -.1, c(" tech succesful", "tech failure", "observed tech") , col=1, horiz=TRUE, cex=1, pch=c(23,4,8), xpd=NA , bty="n" , bg=col.alpha("black" , 0.5) )

dev.off()


PredsPop = array(0,dim=c(75,7))
PredsPop2 = array(0,dim=c(75,7))
PredsMono = array(0,dim=c(75,7,23))

for(i in 1:75){
    for(k in 1:23){
        for (j in 1:7){
            PredsMono[i,j,k] <- mean(Preds[d1$fruit_index[d1$date_index==i & d1$tech==j],j,k])
        }
    }
}

for(i in 1:75){
        for (j in 1:7){
            PredsPop[i,j] <- mean(Preds[d1$fruit_index[d1$date_index==i & d1$tech==j],j,])
        }
    }


pdf("individual_panama_predictions2.pdf",width=8.5,height=11) 
par( mfrow=c(12, 2) , mar=c(1,1,1,1) , oma=c(4,4,.5,.5) )
par(cex = 0.5)
par(tcl = -0.2)
par(mgp = c(2, 0.6, 0))

for(k in 1:23){
        plot(mat[,1]~mat[,8], ylim=c(0,1) , col="white" , xlab="" , ylab="" ,main=""  , cex.axes=0.4)
        #legend("top", inset=0.01, c(tech_id) , fill=col_index,border=col_index, horiz=TRUE,cex=0.75,bty = "n")
#
    for(i in 1:75){
        for (j in 1:7){
            points(PredsMono[i,j,k]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
            lines (PredsMono[i,j,k]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
 
        }
mtext(mono_index$mono[k], side = 3, line = -2, adj = 0.01, cex = .8) 


    }
}
 plot(mat[,1]~mat[,8], ylim=c(0,1) , col="white" , xlab="" , ylab="" ,main="" , xaxt="n" , yaxt="n" , axes=FALSE)
 legend("top", inset=0.1, c(tech_id[1:3]) , fill=col_index[1:3],border=col_index[1:3], horiz=TRUE,cex=1.2,bty = "n")
 legend("bottom", inset=0.1, c(tech_id[4:7]) , fill=col_index[4:7],border=col_index[4:7], horiz=TRUE,cex=1.2,bty = "n")
mtext("Experimental Days (N=75)", side = 1, line = 1.2, cex = 1, outer=TRUE) 
mtext("daily average probability of choosing technique", side = 2, line = 1.2, cex = 1.2 , outer=TRUE) 

#axis(1, at = seq(from=min(mat[,8]) , to=max(mat[,8]) , by = 5 ), labels=F  , tck=-0.01)
#axis(2, at = seq(from=0 , to=1, by = 0.25) ,tck=-0.01 , labels=TRUE )
dev.off()


        plot(mat[,1]~mat[,8], ylim=c(0,1.1) , col="white" , xlab="Experimental Days (N=75)" , ylab="probability of choosing technique" , cex.lab=1.5 )
    for(i in 1:75){
        for (j in 1:7){
            points(PredsPop[i,j]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
            lines( smooth.spline(date_rep, y=PredsPop[,j] , spar=.9) , col=col_index[j] , lw=1)
        }
    }


###old
for(i in 1:75){
        for (j in 1:7){
            PredsPop2[i,j] <- mean(Preds2[d1$fruit_index[d1$date_index==i & d1$tech==j],j])
        }
    }
###trial
d1$row.seq <- seq(1:d1$N)
d$row.seq <- seq(1:nrow(d))

for(i in 1:75){
        for (j in 1:7){
            PredsPop2[i,j] <- mean(Preds2[d$row.seq[d$date_index==i & d$tech_index==j],j])
        }
    }

d[d$date_index==1 & d$tech_index==1,]


PredsPop2[is.nan(PredsPop2)] = 0
date_rep <- seq(1:75)


#raw data
mat[,8] <- c(1:75)

matlo <- mat
for (i in 1:7){ 
    matlo[,i] <- ifelse(matlo[,i]==0, 0.001, matlo[,i])
    matlo[,i] <- ifelse(matlo[,i]==1, 0.999, matlo[,i])
}


###code for figure s3
cairo_pdf("individual_panama_predictions_with_pop.pdf",width=8.5,height=11)
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,25,13,14,15,16,17,18,19,20,21,22,23,24,25), nrow = 13, ncol = 2)
##set up the plot
layout(m)
par( mar=c(1,1,0.6,0.6) , oma=c(2,4,.1,.1) )
par(cex = 0.5)
par(tcl = -0.2)
par(mgp = c(2, 0.6, 0))
        plot(mat[,1]~mat[,8], ylim=c(0,1) , col="white" , xlab="" , ylab="" ,main=""  , cex.axes=0.4)
    for(i in 1:75){
        for (j in 1:7){
            points(PredsPop2[i,j]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 ,  )
            #ss <- smooth.spline(x=date_rep, y=logit(PredsPop2[,j]) , spar=0.94  ) 
            #lines( ss$x, logistic(ss$y), col=col_index[j] , lw=4)
        }
    }

mtext("Population Mean", side = 3, line = -2, adj = 0.01, cex = .8) 

for(k in 1:23){
        plot(mat[,1]~mat[,8], ylim=c(0,1) , col="white" , xlab="" , ylab="" ,main=""  , cex.axes=0.4)
        #legend("top", inset=0.01, c(tech_id) , fill=col_index,border=col_index, horiz=TRUE,cex=0.75,bty = "n")
#
    for(i in 1:75){
        for (j in 1:7){
            points(PredsMono[i,j,k]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
            lines (PredsMono[i,j,k]~i ,col=col.alpha(col_index[j], alpha=0.99) , pch=18 )
 
        }
mtext(mono_index$mono[k], side = 3, line = -2, adj = 0.01, cex = .8) 
mtext(mono_index$yob[k], side = 3, line = -3, adj = 0.01, cex = .5) 


    }
}

 plot(mat[,1]~mat[,8], ylim=c(0,0.8) , col="white" , xlab="" , ylab="" ,main="" , xaxt="n" , yaxt="n" , axes=FALSE)
 legend("top", inset=0.1, c(tech_id[1:3]) , fill=col_index[1:3],border=col_index[1:3], horiz=TRUE,cex=1.3,bty = "n")
 legend("bottom", inset=0.1, c(tech_id[4:7]) , fill=col_index[4:7],border=col_index[4:7], horiz=TRUE,cex=1.3,bty = "n")
mtext("Experimental Days (N=75)", side = 1, line = 0.01, cex = 1.2, outer=TRUE) 
mtext("daily average probability of choosing technique", side = 2, line = 1.2, cex = 1.2 , outer=TRUE) 

#axis(1, at = seq(from=min(mat[,8]) , to=max(mat[,8]) , by = 5 ), labels=F  , tck=-0.01)
#axis(2, at = seq(from=0 , to=1, by = 0.25) ,tck=-0.01 , labels=TRUE )
dev.off()

