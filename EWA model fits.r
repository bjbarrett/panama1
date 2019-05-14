R
require(rethinking)
require(rstan)

#load data
setwd("~/Dropbox/Panama Data/Panama Shared/github stuff") #set your local diirectory to wherever these files are
d <- read.csv("panama_data_14days.csv" , header=TRUE)
mono_index <- read.csv("mono_indexing.csv")

mono_index$mono <- mono_index$monos
ages <- subset(mono_index, select=c(yob,mono) )
ages$age <- 2015-ages$yob 
ages$age.c <-ages$age - mean(ages$age)

##prepare data list for stan model
d_global <- list(
N = nrow(d),																		#length of dataset
J = length( unique(d$mono_index) ),													#number of individuals
K=max(d$tech_index),																#number of processing techniques
tech = d$tech_index,																#technique index
y = cbind( d$y1 , d$y2 , d$y3 , d$y4 , d$y5 , d$y6 , d$y7 ),						#individual processing times for all techniques at each bout N (individual payoff)
s = cbind(d$s1 , d$s2 , d$s3 , d$s4 ,d$s5 ,d$s6 , d$s7 ),							#observed counts of all K techniques to individual J (frequency-dependence)
ps = cbind(d$ps1 , d$ps2 , d$ps3 , d$ps4 ,d$ps5 ,d$ps6 , d$ps7 ),					#observed mean payoffs of all K techniques to individual J (payoff bias)
ks = cbind(d$ks1 , d$ks2 , d$ks3 , d$ks4 ,d$ks5 ,d$ks6 , d$ks7 ),					#observed matrilineal kin cues of all K techniques to individual J (matrilineal kin-bias)
cohos = cbind(d$cos1 , d$cos2 , d$cos3 , d$cos4 ,d$cos5 ,d$cos6 , d$cos7 ),			#observed cohort cues of all K techniques to individual J (age-cohort/similarity bias)
yobs = cbind(d$Yobs1 , d$Yobs2 , d$Yobs3 , d$Yobs4 ,d$Yobs5 ,d$Yobs6 , d$Yobs7 ),   #observed age cues of all K techniques to individual J (age-bias)
press = cbind(d$prs1 , d$prs2 , d$prs3 , d$prs4 ,d$prs5 ,d$prs6 , d$prs7 ),			#observed rank cues of all K techniques to individual J (age-bias)
bout = d$forg_bout,#bout is forg index here											#processing bout unique to individual J
id = as.integer(as.factor(d$mono_index)),											#individual ID
N_effects=8,																		#number of parameters to estimates
age=d$age.c 																		#centered age of individuals
)

#scale payoffs/cues by dividing by max value
d1 <- d_global
d1$yobs <- d1$yobs / max(d1$yobs)
d1$cohos <- d1$cohos / max(d1$cohos)
d1$ks <- d1$ks / max(d1$ks)
d1$ps <- d1$ps / max(d1$ps)
d1$press <- d1$press / max(d1$press)

#parameter list to save in posterior extractions
parlistglobalage=c("lambda" ,"a_id" , "mu" ,"Bpay","Bkin","Bpres","Bcoho","Byob", "fconf", "dev" , "log_lik" , "b_age" , "Rho" , "sigma" )
parlistpay=c("lambda" ,"a_id" , "mu" ,"Bpay", "dev" , "log_lik", "Rho" , "sigma")
parlistfreq=c("lambda" ,"a_id" , "mu" ,"fconf", "dev" , "log_lik", "Rho" , "sigma")
parlistind=c("lambda" ,"a_id" , "mu" , "dev" , "log_lik", "Rho" , "sigma")



#global model with age assuming file is in working directory
fit_global_age <- stan( file = 'PN_social_global_age.stan', data = d1 , 
	iter = 1000, warmup=500 , chains=3, cores=3,
	control=list( adapt_delta=0.98 ) ,pars=parlistglobalage )

###simpler models not including age

#individual learning
fit_i <- stan( file = 'PN_indiv.stan', data = d3 , 
	iter = 3000, warmup=1500, chains=1, cores=1, pars=parlistind,
	control=list( adapt_delta=0.98 ) )
#frequency dependent-sea
fit_freq <- stan( file = 'PN_social_freq.stan', data = d4 , 
	iter = 3000, warmup=1500 , chains=3, cores=3, pars=parlistfreq,
	control=list( adapt_delta=0.98 ) )
#payoff/cue bias 
fit_pay <- stan( file = 'PN_social_pay.stan', data = d4 , 
	iter = 3000, warmup=1500 , chains=1, cores=1, pars=parlistpay,	
	control=list( adapt_delta=0.98 ) )

#frequency-dependent no age effects

#payoff/cur biased no age effects

#individual learning no age effects

###traceplots
traceplot(fit_global_age, pars=c("mu" , "lambda"))
traceplot(fit_global_age, pars=c("b_age"))
traceplot(fit_global_age, pars=c("sigma"))

post <- extract(fit_global_age) #extract posterior samples from model

###########main effects posterior graphs aka fig s1##################
cairo_pdf("Main_effects_postgraphs.pdf",height=8,width=8)
par(mfrow=c(3,3), oma = c(0, 0, 0, 0) , mar=c(3,.5,2,.5))
par(cex = 0.6)
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

dens(logistic(post$mu[,1]) , main=expression(a.~ weight~of~new~experience~(paste(phi)) ) , xlim=c(0,1), ylab='' , xlab= '' , col="white", yaxt='n')##phi
abline(v=median(logistic(post$mu[,1])) , col="cornflowerblue" ) 
shade( density(logistic(post$mu[,1])) , lim= as.vector(HPDI(logistic(post$mu[,1]), prob=0.9999)) , col = col.alpha("cornflowerblue", 0.25))
#curve( logistic(dnorm( x , 0 , 1)) , lty=2 , add=TRUE , col="cornflowerblue" )
dens( logistic(rnorm( 10000 , 0 , 1)) , lty=2  , col="cornflowerblue" , adj=1 , add=TRUE)

dens(logistic(post$mu[,2]) , main=expression(b.~weight~of~social~information~(paste(gamma)))  ,  xlim=c(0,1) , ylab='', xlab= '' , col="white", yaxt='n') #gamma
abline(v=median(logistic(post$mu[,2])) , col="cornflowerblue") 
shade( density(logistic(post$mu[,2])) , lim= as.vector(HPDI(logistic(post$mu[,2]), prob=0.9999)) , col = col.alpha("cornflowerblue", 0.25))
#curve( logistic(dnorm( x , 0 , 1)) , lty=2 , add=TRUE , col="cornflowerblue" )
dens( logistic(rnorm( 10000 , 0 , 1)) , lty=2  , col="cornflowerblue" , adj=1 , add=TRUE)

dens(exp(post$mu[,3]) , main=expression(c.~strength~of~frequency~dependence~(paste( "\u0192"[c]))),  xlim=c(0,4) , xlab="c. strength of frequency dependence" , col="white", ylab='', yaxt='n')##fconf
abline(v=median(exp(post$mu[,3])) , col="red" ) #fconf
shade( density(exp(post$mu[,3])) , lim= as.vector(HPDI(exp(post$mu[,3]), prob=0.9999)) , col = col.alpha("red", 0.25))
#curve( exp(dnorm( x , 0 , 1)) , lty=2 , add=TRUE , col="cornflowerblue" )
dens( exp(rnorm( 10000 , 0 , 1)) , lty=2  , col="red" , adj=1 , add=TRUE)

dens(post$mu[,4]  ,main=expression(d.~strength~of~payoff~bias~(paste(beta)[pay])) ,  xlim=c(-4,4), xlab='' , ylab='',col="white", yaxt='n')##fpay
abline( v=median(post$mu[,4]) , col="orange"  )
shade( density(post$mu[,4]) , lim= as.vector(HPDI((post$mu[,4]), prob=0.9999)) , col = col.alpha("orange", 0.25))
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )
#dens( rnorm( 10000 , 0 , 1) , lty=2  , col="orange" , adj=1 , add=TRUE)

dens(post$mu[,5]  ,main=expression(e.~strength~of~kin~bias~(paste(beta)[kin])) ,  xlim=c(-4,4), ylab='',xlab='', col="white", yaxt='n')##fpay
abline( v=median(post$mu[,5]) , col="orange" ) 
shade( density(post$mu[,5]) , lim= as.vector(HPDI((post$mu[,5]), prob=0.9999)) , col = col.alpha("orange", 0.25))
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(post$mu[,6] ,main=expression(f.~strength~of~rank~bias~(paste(beta)[rank])) ,  xlim=c(-4,4), ylab='', col="white", xlab='', yaxt='n')##fpay
abline( v=median(post$mu[,6]) , col="orange" )
shade( density(post$mu[,6]) , lim= as.vector(HPDI((post$mu[,6]), prob=0.9999)) , col = col.alpha("orange", 0.25))
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(post$mu[,7]  ,main=expression(g.~strength~of~age-cohort~bias~(paste(beta)[coho])) ,  xlim=c(-4,4), ylab='', xlab='' , col="white", yaxt='n')##fpay
abline( v=median(post$mu[,7]) , col="orange" ) 
shade( density(post$mu[,7]) , lim= as.vector(HPDI((post$mu[,7]), prob=0.9999)) , col = col.alpha("orange", 0.25))
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(post$mu[,8]  ,main=expression(h.~strength~of~age~bias~(paste(beta)[age])) , col="white", xlim=c(-4,4), ylab='',xlab='', yaxt='n')
abline( v=median(post$mu[,8]) , col="orange" ) 
shade( density(post$mu[,8]) , lim= as.vector(HPDI((post$mu[,8]), prob=0.9999)) , col = col.alpha("orange", 0.25))
#curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )
curve( dnorm( x , 0 , 1) , lty=2 , add=TRUE , col="orange" )

dens(as.vector(post$lambda) , main=expression(i.~sensitivity~to~individual~payoff~(paste(lambda))) , ylab='', xlab='' , col="white" , yaxt='n', xlim=c(0,25) )
abline( v=median(post$lambda) , col="black" ) 
shade( density(post$lambda) , lim= HPDI(as.vector(post$lambda), prob=0.9999) , col = col.alpha("black", 0.25) )
curve( dexp( x , rate=1) , lty=2 , add=TRUE , col="black" )
#dens( rexp (10000 , rate=1/5) , lty=2  , col="black" , adj=1 , add=TRUE)
#dens( rexp (10000 , rate=1/5) , lty=2  , col="black" , adj=1 , add=TRUE)

dev.off()

##########estimates of sigma graph aka fig s2########################
cairo_pdf("varef_sigma.pdf",height=4,width=7.5)

col.sig <- c("cornflowerblue","cornflowerblue","red", "orange" , "orange", "orange" , "orange" , "orange" )
par(mfrow=c(2,4), oma = c(0, 0, 0, 0) , mar=c(3,.5,2,.5))
par(cex = 0.6)
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
lab.name <- c("\u03A6","\u03B3","\u0192c" , "\u03B2pay" , "\u03B2kin" , "\u03B2rank" , "\u03B2coho" , "\u03B2age")
for (i in 1:8){
	dens(post$sigma[,i] , xlim=c(0,4) ,  col="white", yaxt='n'  , main= bquote( sigma [.(lab.name[i])])  )
	curve( dexp( x , rate=3) , lty=2 , add=TRUE , col=col.sig[i] )
	shade( density(post$sigma[,i]) , lim= HPDI(as.vector(post$sigma[,i]), prob=0.999) , col = col.alpha(col.sig[i], 0.25) )
}
dev.off()


##########################phi age effects graph fig 3a################
post <- extract(fit_global_age)
vfphi <- matrix(0,nrow=length(post$lambda),ncol=23)
varefphi <- rep(0,23)
for(i in 1:23){vfphi[,i] <- logistic(post$mu[,1] + post$a_id[,i,1] + post$b_age[,1]*ages$age.c[i] ) }
for(i in 1:23){varefphi[i] <- mean(vfphi[,i])}

rando.samps <- sample(1:nrow(post$lambda), size=100, replace = FALSE, prob = NULL)
     
     sample.int(n, size = n, replace = FALSE, prob = NULL)
cairo_pdf("phi_age_varef.pdf", width=8 , height=8)
par(mar=c(5,5,0.5,0.5))
plot(varefphi~ ages$age.c , ylab="attraction toward new experience (\u0278)" , xlab="age (years)" , pch=19 , col="orange" ,xlim=c( (min(ages$age.c)-1) , (max(ages$age.c) + 2)), ylim=c(0,.6) , cex=1.5 , cex.lab=2.4  , xaxt="n" , yaxt="n" , axes=FALSE)

axis(1, at = seq(from=(min(ages$age.c)-1) , to=(max(ages$age.c) + 2), by = 5) , labels=seq(from=0 , to=25 , by=5), tck=-0.02 , cex.axis=1.5)
axis(2, at = seq(from=0 , to=0.6, by = 0.2) , tck=-0.02 , cex.axis=1.5)
axis(1, at = seq(from=(min(ages$age.c)-1) , to=(max(ages$age.c) + 2), by = 1 ), labels=F  , tck=-0.01, cex.axis=1.5)
axis(2, at = seq(from=0 , to=0.6, by = 0.1) ,tck=-0.01 , labels=F , cex.axis=1.5)

age.seq <- seq(from=min(ages$age.c) , to=max(ages$age.c) , length=25 )
pred.mean <- sapply(age.seq , function(z)
	 mean(logistic(post$mu[,1] + post$b_age[,1]*z) ))

pred.ci<- sapply(age.seq , function(z)
	 HPDI(logistic(post$mu[,1] + post$b_age[,1]*z) ))

pred.sims <- sapply(age.seq , function(z)
	 logistic(post$mu[,1] + post$b_age[,1]*z) )

for (i in rando.samps){lines(age.seq , pred.sims[i,] , lw=2 , col=col.alpha( "orange" , alpha = 0.1 ))}
text(ages$age.c, varefphi, mono_index$mono, cex=.4, col="black")
lines(age.seq , pred.mean , lw=2)
dev.off()


##############gamma age effects graph fig 3b###################

vfgamma <- matrix(0,nrow=length(post$lambda),ncol=23)
varefgamma <- rep(0,23)
for(i in 1:23){vfgamma[,i] <- logistic(post$mu[,2] + post$a_id[,i,2] + post$b_age[,2]*ages$age.c[i] ) }
for(i in 1:23){varefgamma[i] <- mean(vfgamma[,i])}
lines(age.seq , pred.mean , lw=2)

cairo_pdf("gamma_age_varef.pdf", width=8 , height=8)
par(mar=c(5,5,0.5,0.5))
plot(varefgamma~ ages$age.c , ylab="weight given to social information (\u03B3)" , xlab="age (years)" , pch=19 , col="cornflowerblue" ,xlim=c( (min(ages$age.c)-1) , (max(ages$age.c) + 2)), ylim=c(0,.6) , cex=1.7 , cex.lab=2.4  , xaxt="n" , yaxt="n" , axes=FALSE )
axis(1, at = seq(from=(min(ages$age.c)-1) , to=(max(ages$age.c) + 2), by = 5) , labels=seq(from=0 , to=25 , by=5), tck=-0.02, cex.axis=1.5)
axis(2, at = seq(from=0 , to=0.6, by = 0.2) , tck=-0.02 , cex.axis=1.5)
axis(1, at = seq(from=(min(ages$age.c)-1) , to=(max(ages$age.c) + 2), by = 1 ), labels=F  , tck=-0.01, cex.axis=1.5)
axis(2, at = seq(from=0 , to=0.6, by = 0.1) ,tck=-0.01 , labels=F , cex.axis=1.5)
age.seq <- seq(from=min(ages$age.c) , to=max(ages$age.c) , length=25 )
pred.mean <- sapply(age.seq , function(z)
	 mean(logistic(post$mu[,2] + post$b_age[,2]*z) ))
pred.ci<- sapply(age.seq , function(z)
	 HPDI(logistic(post$mu[,2] + post$b_age[,2]*z) ))

pred.sims <- sapply(age.seq , function(z)
	 logistic(post$mu[,2] + post$b_age[,2]*z) )
for (i in rando.samps){lines(age.seq , pred.sims[i,] , lw=2 , col=col.alpha( "cornflowerblue" , alpha = 0.1 ))}
text(ages$age.c, varefgamma, mono_index$mono, cex=.4, col="black")
lines(age.seq , pred.mean , lw=2)

dev.off()
