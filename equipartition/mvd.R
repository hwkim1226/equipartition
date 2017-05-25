# 1. find velocity dispersion as a function of mass, and plot it. (log-log plots)
# 2. fit power-law and Bianchini model, and write the table of coefficients.
library(data.table)

# round keeping the last 3 digits
round1000 <- function(a)
{
	return(0.001*round(1000*a))
}

# for last 4 digits
round10000 <- function(a)
{
	return(0.0001*round(10000*a))
}

# In the following, we will calculate the velocity dispersion. This function calculates the error of the velocity dispersion by bootstrap sampling.
# Ntry is the number of resampling runs. (the higher, the more accurate)
Ntry <- 1000
count <- 0
esd <- function(x) #estimate of the uncertainty of standard error estimate using resampling
{
        sd_i <- numeric(0)
        for (i in 1:Ntry)
        {
                sd_i <- c(sd_i, sd(sample(x, replace = T)))
        }
        return(sd(sd_i)) #do not divide by sqrt(Ntry) because we want the uncertainty on a single measurement (we measure sd once on the true dataset, we do not take the mean of the sds of the resampled datasets)
}

countn <- function(x)
{
	count <<- count + 1
}

# number of bins in mass
step <- 8

x <- matrix(ncol=step)
y <- matrix(ncol=step)
err <- matrix(ncol=step)
logx <- matrix(ncol=step)
logy <- matrix(ncol=step)
logerr <- matrix(ncol=step)
plotmvd <- function(datfile) # function to make velocity dispersion data for a given file c_????.dat
{
	
	gcdat <- read.table(datfile, header = T) #read the file
	
	countn()

	minmass <- min(sort(unique(gcdat$V7))) #find a minimum of mass
	maxmass <- max(sort(unique(gcdat$V7))) #find a maximum of mass		

	qgrid <- seq(from = 0.0, to = 1.0, by = 1.0/step)

	grid <- quantile(gcdat$V7, qgrid) # each bins has the equal number of stars.
	m <- grid[1:step]
	grid[step+1] <- 1.01*grid[step+1] # expanding a little bit to include all stars

	lgrid <- quantile(log10(gcdat$V7), qgrid) # same thing with log...	
	logm <- lgrid[1:step]
	lgrid[step+1] <- 1.01*lgrid[step+1]

# centering and calculating velocities
	vx <- gcdat$V4 - median(gcdat$V4)
	vy <- gcdat$V5 - median(gcdat$V5)
	vz <- gcdat$V6 - median(gcdat$V6)
	vel <- sqrt( (vx)**2 + (vy)**2 + (vz)**2 )

# lvd will hold the log of velocity dispersion. lerrvd will hold its error. lbin is used to select stars inside the bin.
	lbin <- numeric(0)
	lvd <- numeric(0)
	lerrvd <- numeric(0)
	bin <- numeric(0)
	vd <- numeric(0)
	errvd <- numeric(0)	

### for linear plot
	#for (i in 1:step)
	#{
	#	bin <- gcdat$V7 >= grid[i] & gcdat$V7 < grid[i+1]
	#	vd <- c(vd, sd(vel[bin]) )
	#	errvd <- c(errvd, esd(vel[bin]) )
	#}

### for log plot
	for (i in 1:step)
	{
		#lbin <- log10(gcdat$V7) >= lgrid[i] & log10(gcdat$V7) < lgrid[i+1]
		#lvd <- c(lvd, log10( sd(vel[lbin]) ) )
		#lerrvd <- c(lerrvd, esd(vel[lbin])/(log(10.0)*sd(vel[lbin])) )
		bin <- gcdat$V7 >= grid[i] & gcdat$V7 < grid[i+1] # selecting the stars inside the radial bin
		lvd <- c(lvd, log10( sd(vel[bin]) ) ) # calculating the velocity dispersion
		eba <- esd(vel[bin])/(log(10.0)*sd(vel[bin])) # error bars
		lerrvd <- c(lerrvd, eba )
	}

# putting them together
	x <<- rbind(x, m)
	y <<- rbind(y, vd)
	err <<- rbind(err, errvd)

	logx <<- rbind(logx, logm)
	logy <<- rbind(logy, lvd)
	logerr <<- rbind(logerr, lerrvd)
}
# in the original files, the function was defined twice.

# choose the alphas for which to run the plots and the fit
alpha <- seq(0.5, 3.0, 0.1)
folder <- numeric(0)

# read folders corresponding to alphas
for (i in 1:length(alpha))
{

	folder[i] <- paste("../10times/",as.character(10*alpha[i]),"/all", sep="")	

	if( alpha[i] < 1.0 || (alpha[i] > 1.0 && alpha[i] < 1.4) || alpha[i] > 2.5 )
	{
		folder[i] <- paste("../isolated/32k",as.character(10*alpha[i]),"/all", sep="")
	}
}


for (j in 1:length(alpha)) # plot and fit for each alphas
{
	datfiles <- system(paste("ls ",folder[j],"/*.dat", sep=""), intern = T) # this command executes shell commands, intern = T means it returns the output of the command (if you remove intern = T datfiles would be empty)
	sapply(datfiles, plotmvd) # applies a function to the elements of an array one by one, in this case applies plotvel to every dat file name

	colors <- rainbow(6)

	x <<- x[-1,]
	y <<- y[-1,]
	err <<- err[-1,]
	logx <<- logx[-1,]
	logy <<- logy[-1,] 
	logerr <<- logerr[-1,]

	fn_logx <- paste("logx_",alpha[j] ,sep="")
	fn_logy <- paste("logy_",alpha[j] ,sep="")
	fn_logerr <- paste("logerr_",alpha[j] ,sep="")
# write tables (data for plot)
	write.table(logx, fn_logx)
	write.table(logy, fn_logy)
	write.table(logerr, fn_logerr)


# log-log

	mx <- min(logx) - 0.1 * ( max(logx) - min(logx) )
	Mx <- max(logx) + 0.2 * ( max(logx) - min(logx) )
	my <- min(logy) - 0.1 * ( max(logy) - min(logy) )
	My <- max(logy) + 0.2 * ( max(logy) - min(logy) )

	pdf(file=paste(as.character(10*alpha[j]),"_mvd.pdf",sep=""))
	plot.new()
	par(mar=c(10,10,8,8))
	plot.window( xlim=c(-6, -4), ylim=c(-0.6, -0.4 ))

	axis(1); axis(2, las=1); box()

	title(main="", sub="",
		xlab='log m', ylab=expression(paste("log ", sigma, sep="")),
		mgp=c(5,1,0),
		cex.lab=2)

	k <- 0

	for (i in 1:count)
	{
		par(new=T)

		fit <- lm(logy[i,] ~ logx[i,], weights=(1/(logerr[i,])**2) )
		fit2 <- lm(logy[i,] ~ x[i,], weights=(1/(logerr[i,])**2))
		a2=fit2$coefficients[[1]]
		b2=fit2$coefficients[[2]]
		xs2 <- seq(from = 1E-7, to = 1E-3, length.out = 500000)
	    ys2 <- a2 + b2*xs2
	
#		B = -1/2meq
		coefftable <- data.frame(c("power-law","exponential"), c(round1000(fit$coefficients[[2]]), round10000(-0.5/fit2$coefficients[[2]])), c(round1000(fit$coefficients[[1]]), round1000(fit2$coefficients[[1]])), c(datfiles[i], datfiles[i]))

		cota <- paste("coefftable_", alpha[j], sep = "")
		write.table(coefftable, cota, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# only one in four is actually plotted to reduce confusion in the plot.
		if(i %% 4 == 1)
		{
			k <- k + 1
#			lines(logx[i,], logy[i,], col=colors[k])
			points(logx[i,], logy[i,], pch=20, col=colors[k])
			arrows(logx[i,], logy[i,], logx[i,], logy[i,] + logerr[i,], length = 0.04, angle = 90, col=colors[k])
			arrows(logx[i,], logy[i,], logx[i,], logy[i,] - logerr[i,], length = 0.04, angle = 90, col=colors[k])

			abline(a=fit$coefficients[[1]], b=fit$coefficients[[2]], col=colors[k])
			lines(log10(xs2), ys2, col = colors[k], lty = 2, lwd = 2) 
		}
	}



name <- numeric(0)

dev.off()

count <- 0
x <<- matrix(ncol=step)
y <<- matrix(ncol=step)
err <<- matrix(ncol=step)
logx <<- matrix(ncol=step)
logy <<- matrix(ncol=step)
logerr <<- matrix(ncol=step)

}

#warnings()

