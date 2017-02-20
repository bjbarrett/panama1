require(rethinking)
require(truncnorm)
setwd("~/github stuff")

Softmax <- function(x){
exp(x)/sum(exp(x))
} #softmax function to simplify code
do <- read.csv("~/Dropbox/Panama Data/Panama Shared/github stuff/panama_data_14days.csv" , header=TRUE)

#data sims
n <- 25 # number of individuals
nbouts <- 70 #rounds foraged
techmeans <- c( mean(do$y1[do$tech_index==1 & do$open==1]),
            mean(do$y2[do$tech_index==2 & do$open==1]),
            mean(do$y3[do$tech_index==3 & do$open==1]),
            mean(do$y4[do$tech_index==4 & do$open==1]) ) # list of tech means for successes from data

techvar <- c( sd(do$y1[do$tech_index==1 & do$open==1]),
            sd(do$y2[do$tech_index==2 & do$open==1]),
            sd(do$y3[do$tech_index==3 & do$open==1]),
            sd(do$y4[do$tech_index==4 & do$open==1]) ) #list of tech variances for successes from data

techprsucceed <- c(0.511,0.378,0.885,0.665) #probability technique leads to success


#parameter sims
fc.sim <- log(.4)				#frequency dependecy parameter on log scale
phi.sim <- logit(.15) ## attraction updating parameter on log-odds scale
gamma.sim <- logit(.14) ## influence of social info parameter on log-odds scale
k.lambda <- 20           ##sensitivity to individual payoffs
beta.p <- 1.1          ##contribution of payoff bias cues

#varying effects offsets for individuals
gamma.sim_i <- rnorm( n , mean=0 , sd=.69 ) 
phi.sim_i <- rnorm( n , mean=0 , sd=.66 ) 
fc.sim_i <- rnorm( n , mean=0 , sd=.15 ) 
beta.p_i <- rnorm( n , mean=0 , sd=.2 )

logistic(gamma.sim+gamma.sim_i) ##simulated gammas for all n individuals
logistic(phi.sim+phi.sim_i) ##simulated phis for all n individuals
exp(fc.sim + fc.sim_i) ##simulated fc for all individuals
beta.p + beta.p_i ##simulated fc for all individuals

#simulate data
dsim_s <- data.frame( i=0 , bout=0 , tech=0 , y1=0 , y2=0, y3=0 , y4=0 , s1=0 , s2=0 , s3=0 , s4=0, ps1=0 , ps2=0 , ps3=0 , ps4=0 , A1=0 , A2=0 , A3=0 , A4=0)
therow <- 1
AC <- matrix(1,ncol=4,nrow=n) #attraction scores for each tech
AC[,3] <- 0 #all have equal prob of other behaviors, except option 3
AC[1,] <- c(.5,.5,.6,.5) #individual 1 does option 3 most of the time, it is best option
S1 <- S2 <- S3 <- S4 <- rep(0,n+1) # num of individuals choosing each tech in previous bout
PS1 <- PS2 <- PS3 <- PS4 <- rep(0,nbouts+1) # empty vector for mean observed in previous rounds
lin_mod <- s_temp <-  rep(0,4)
for ( r in 1:nbouts ) {
    for ( i in 1:n ) {  

		prtech_i <-  Softmax(k.lambda*AC[i,]) 
        my.bp <- beta.p + beta.p_i[i]  #payoff weight for individual i
        my.gam <- logistic( gamma.sim + gamma.sim_i[i] ) #social info weight for individual i
        my.phi <- logistic( phi.sim + phi.sim_i[i] ) #social info weight for individual i
        my.fconf <- exp( fc.sim + fc.sim_i[i]) 
        prtech_sp <- c(PS1[r],PS2[r],PS3[r],PS4[r]) #how will this fit into 0 payoffs
        prtech_su <- c(S1[r],S2[r],S3[r],S4[r])

        #//conformity aspect below
        if ( r > 1 ) {
            if (sum( prtech_su ) > 0 ) {

                #// compute non-frequency cues as log-linear model
                for ( j in 2:4 ) {
                    lin_mod[j] <- exp( my.bp*prtech_sp[j])
                }
                lin_mod[1] <- 1 #// aliased outcome

                #// compute frequency cue
                for ( j in 1:4 ){ s_temp[j] <- prtech_su[j]^my.fconf}
                for ( j in 1:4 ){lin_mod[j] <- lin_mod[j] * s_temp[j]}

                prtech_s <- lin_mod/sum(lin_mod)
                prtech <- (1-my.gam)*prtech_i + my.gam*prtech_s

            } else {
                prtech <- prtech_i
            }
        } else {
            prtech <- prtech_i
         }
# choose tech
        tech <- sample( 1:4 , size=1 , prob=prtech)
        yield <- rtruncnorm( 1 , a=0 , b=Inf, mean=techmeans[tech] , sd=techvar[tech] )
# update attractions
        yields <- rep(0,4)
        succeed <- sample(rbinom(1,1,techprsucceed[tech]))
        succeed_i <- sample(rbinom(1,1,techprsucceed[tech]))

        yields[tech] <- succeed*(yield) #makes payoff 1/time if succeed, 0 if not
        #yields[tech] <- succeed#makes payoff 1 if succeed, 0 if not
        #yields[tech] <- yield#makes payoff yield
        for (k in 1:4){
        	AC[i,k] <- (1-my.phi)*AC[i,k] + my.phi*yields[k]
        }

        dsim_s[therow,] <- c( i , r , tech , yields[1] , yields[2] , yields[3] , yields[4] , S1[r] , S2[r] , S3[r] , S4[r], PS1[r] , PS2[r] , PS3[r] , PS4[r] , AC[i,1] , AC[i,2] , AC[i,3] , AC[i,4])
        therow <- therow + 1
    } #i
    PS1[r+1] <- ifelse(is.nan(mean( dsim_s$y1[dsim_s$bout==r & dsim_s$tech==1] )), 0, mean( dsim_s$y1[dsim_s$bout==r & dsim_s$tech==1] ) )
    PS2[r+1] <- ifelse(is.nan(mean( dsim_s$y2[dsim_s$bout==r & dsim_s$tech==2] )), 0, mean( dsim_s$y2[dsim_s$bout==r & dsim_s$tech==2] ) )
    PS3[r+1] <- ifelse(is.nan(mean( dsim_s$y3[dsim_s$bout==r & dsim_s$tech==3] )), 0, mean( dsim_s$y3[dsim_s$bout==r & dsim_s$tech==3] ) )
    PS4[r+1] <- ifelse(is.nan(mean( dsim_s$y4[dsim_s$bout==r & dsim_s$tech==4] )), 0, mean( dsim_s$y4[dsim_s$bout==r & dsim_s$tech==4] ) )
    S1[r+1] <- length( dsim_s$tech[dsim_s$tech==1 & dsim_s$bout==r] )
    S2[r+1] <- length( dsim_s$tech[dsim_s$tech==2 & dsim_s$bout==r] )
    S3[r+1] <- length( dsim_s$tech[dsim_s$tech==3 & dsim_s$bout==r] )
    S4[r+1] <- length( dsim_s$tech[dsim_s$tech==4 & dsim_s$bout==r] )

 }

o <- order( dsim_s$i )
dsim2 <- dsim_s[o,]

#plot simulated data
plot(s1/n ~ bout, data=dsim2, col="red" , ylim=c(0,1) , pch=19 , xlab="Time (Foraging Bouts)" , ylab="Proportion of Individuals Choosing Option" , main="4 options, green higest payoff", xlim=c(2,70) )
points(s2/n ~ bout, data=dsim2 , col="gold", pch=19)
points(s3/n ~ bout, data=dsim2 , col="green", pch=19)
points(s4/n ~ bout, data=dsim2 , col="blue", pch=19)

#turn simulated data to a list for r-stan
ds <- list(
N = nrow(dsim2),
J = length( unique(dsim2$i)),
K=max(dsim2$tech),
tech = dsim2$tech,
y = cbind( dsim2$y1 , dsim2$y2 , dsim2$y3 , dsim2$y4 ),
s = cbind(dsim2$s1 , dsim2$s2 , dsim2$s3 , dsim2$s4 ),
ps = cbind(dsim2$ps1 , dsim2$ps2 , dsim2$ps3 , dsim2$ps4)  ,
id = dsim2$i,
bout = dsim2$bout ,
N_effects=4
)
parlistcombo=c("lambda" ,"a_id" , "mu" ,"Bpay", "fconf", "dev" , "log_lik")

#fit model in stan
fit_combo <- stan( file = 'PN_social_combo.stan', data = ds , 
    iter = 2000, warmup=1000, chains=2, cores=3, pars=parlistcombo, 
    control=list( adapt_delta=0.98 ) )

post <- extract(fit_combo) # extract posterior

#plot to recover main effects
par(mfrow=c(2, 2), oma = c(0, 0, 0, 0) , mar=c(3,.5,2,.5))
par(cex = 0.6)
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

dens(logistic(post$mu[,1]) , main=expression(paste(phi)) , xlim=c(0,1), ylab='' , xlab= "a. weight of new experience" , col="white", yaxt='n' , cex.lab=1.5)##phi
abline( v=logistic(phi.sim) , col="cornflowerblue" ) 
shade( density(logistic(post$mu[,1])) , lim= as.vector(HPDI(logistic(post$mu[,1]), prob=0.9999)) , col = col.alpha("cornflowerblue", 0.25))

dens(logistic(post$mu[,2]) , main=expression(paste(gamma))  ,  xlim=c(0,1) , ylab='', xlab= "b. weight of social information" , col="white", yaxt='n', , cex.lab=1.5) #gamma
abline(v=logistic(gamma.sim) , col="cornflowerblue") 
shade( density(logistic(post$mu[,2])) , lim= as.vector(HPDI(logistic(post$mu[,2]), prob=0.9999)) , col = col.alpha("cornflowerblue", 0.25))

dens(exp(post$mu[,3]) , main=expression(paste( "\u0192"[c])),  xlim=c(0,4) , xlab="c. strength of frequency dependence" , col="white", ylab='', yaxt='n',  cex.lab=1.5)##fconf
abline(v=exp(fc.sim) , col="red" ) #fconf
shade( density(exp(post$mu[,3])) , lim= as.vector(HPDI(exp(post$mu[,3]), prob=0.9999)) , col = col.alpha("red", 0.25))

dens(post$mu[,4]  ,main=expression(paste(beta)[pay]) ,  xlim=c(-4,4), xlab="d. strength of payoff bias" , ylab='',col="white", yaxt='n', cex.lab=1.5)##fpay
abline( v=beta.p , col="orange"  ) 
shade( density(post$mu[,4]) , lim= as.vector(HPDI((post$mu[,4]), prob=0.9999)) , col = col.alpha("orange", 0.25))

dens(as.vector(post$lambda) , main=expression(paste(lambda)) , ylab='', xlab="i. sensitivity to individual payoff" , col="white" , yaxt='n', cex.lab=1.5)
abline( v=k.lambda , col="black" ) 
shade( density(post$lambda) , lim= HPDI(as.vector(post$lambda), prob=0.9999) , col = col.alpha("black", 0.25) )

dens(post$lambda)

#plot varying effects (recommend exploring parameter space to see under what conditions these can be recovered. 
#if task is too easy to learn individually, varying efects for gamma are hard to recover
#if behavior changes quickly at population level and the past is not different from recent behavior varying effects of phi are hard to recover

gam.k <- logistic(gamma.sim + gamma.sim_i)
gam.pred <- rep(0,n)
for(i in 1:n){gam.pred[i] <- mean(logistic(post$mu[,2] + post$a_id[,i,2]))}
plot(gam.k,gam.pred , pch=19 , col="orange" , xlim=c(0,0.7) , ylim=c(0,0.7) )
abline(a = 0, b = 1)


phi.k<- logistic(phi.sim + phi.sim_i)
phi.pred <- rep(0,n)
for(i in 1:n){phi.pred[i] <- median(logistic(post$mu[,1] + post$a_id[,i,1]))}
plot(phi.k,phi.pred , pch=19 , col="slateblue" , xlim=c(0,0.6) , ylim=c(0,0.6) )
abline(a = 0, b = 1)

