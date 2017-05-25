### find velocity dispersion as a function of mass, and plot it. (log-log plots)
library(data.table)

round1000 <- function(a)
{
	return(0.001*round(1000*a))
}

round10000 <- function(a)
{
	return(0.0001*round(10000*a))
}

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

step <- 8

x <- matrix(ncol=step)
y <- matrix(ncol=step)
err <- matrix(ncol=step)

logx <- matrix(ncol=step)
logy <- matrix(ncol=step)
logerr <- matrix(ncol=step)


plotmvd <- function(datfile) #function to make a plot for a given file c_????.dat
{
	
	gcdat <- read.table(datfile, header = F) #read the file
	countn()

	minmass <- min(sort(unique(gcdat$V7))) #find a minimum of mass
	maxmass <- max(sort(unique(gcdat$V7))) #find a maximum of mass		

	qgrid <- seq(from = 0.0, to = 1.0, by = 1.0/step)

	grid <- quantile(gcdat$V7, qgrid)
	m <- grid[1:step]
	grid[step+1] <- 1.01*grid[step+1]

	lgrid <- quantile(log10(gcdat$V7), qgrid) #make a grid		
	logm <- lgrid[1:step]
	lgrid[step+1] <- 1.01*lgrid[step+1]


	vx <- gcdat$V4 - median(gcdat$V4)
	vy <- gcdat$V5 - median(gcdat$V5)
	vz <- gcdat$V6 - median(gcdat$V6)
	vel <- sqrt( (vx)**2 + (vy)**2 + (vz)**2 ) #calculate velocities


	lbin <- numeric(0)
	lvd <- numeric(0)
	lerrvd <- numeric(0)

	bin <- numeric(0)
	vd <- numeric(0)
	errvd <- numeric(0)	

	for (i in 1:step)
	{
### for linear plot
		bin <- gcdat$V7 >= grid[i] & gcdat$V7 < grid[i+1]
		vd <- c(vd, sd(vel[bin]) )
		errvd <- c(errvd, esd(vel[bin]) )
	}

		for (i in 1:step)
	{
### for log plot
		lbin <- log10(gcdat$V7) >= lgrid[i] & log10(gcdat$V7) < lgrid[i+1]
		lvd <- c(lvd, log10( sd(vel[lbin]) ) )
		lerrvd <- c(lerrvd, esd(vel[lbin])/(log(10.0)*sd(vel[lbin])) )
	}
       

	x <<- rbind(x, m)
	y <<- rbind(y, vd)
	err <<- rbind(err, errvd)

	logx <<- rbind(logx, logm)
	logy <<- rbind(logy, lvd)
	logerr <<- rbind(logerr, lerrvd)
}

datfiles <- system("ls c_????.dat", intern = T) #this command executes shell commands, intern = T means it returns the output of the command (if you remove intern = T datfiles would be empty)
sapply(datfiles, plotmvd) #applies a function to the elements of an array one by one, in this case applies plotvel to every dat file name
count

colors <- rainbow(20)
colors <- c(colors[1], colors[5], colors[10], colors[15], colors[20])

x <<- x[-1,]
y <<- y[-1,]
err <<- err[-1,]

logx <<- logx[-1,]
logy <<- logy[-1,] 
logerr <<- logerr[-1,]


# plot 1: log-log
mx <- min(logx) - 0.1 * ( max(logx) - min(logx) )
Mx <- max(logx) + 0.2 * ( max(logx) - min(logx) )
my <- min(logy) - 0.1 * ( max(logy) - min(logy) )
My <- max(logy) + 0.2 * ( max(logy) - min(logy) )



numb <- system("pwd | tr '/' '\n' | tail -1", intern=T)
pdf(file=paste(numb,"_mvd.pdf",sep=""))
plot.new()
par(mar=c(10,10,8,8))
#par(cex.axis=1.5)
plot.window(xlim=c(mx,Mx), ylim=c(my,My))

axis(1); axis(2, las=1); box()

title(main="", sub="",
	  xlab='log m', ylab=expression(paste("log", sigma, sep="")),
	  mgp=c(5,1,0))


for (i in 1:count)
	{
	par(new=T)

# fit
	points(logx[i,], logy[i,], pch=16, col=colors[i])
	fit <- lm(logy[i,] ~ logx[i,], weights=(1/(logerr[i,])**2) )
	chi <- sum((fit$residuals/logerr[i,])**2)
	fit2 <- lm(logy[i,] ~ x[i,], weights=(1/(logerr[i,])**2))
	chi2 <- sum((fit2$residuals/logerr[i,])**2)
	abline(a=fit$coefficients[[1]], b=fit$coefficients[[2]], col=colors[i])
	a2=fit2$coefficients[[1]]
	b2=fit2$coefficients[[2]]
	xs2 <- seq(from = 1E-7, to = 1E-3, length.out = 500000)
    ys2 <- a2 + b2*xs2
	lines(log10(xs2), ys2, col = colors[i], lty = 2, lwd = 2) 
	
#errorbar
	arrows(logx[i,], logy[i,], logx[i,], logy[i,] + logerr[i,], length = 0.1, angle = 90, col=colors[i])
	arrows(logx[i,], logy[i,], logx[i,], logy[i,] - logerr[i,], length = 0.1, angle = 90, col=colors[i])
	}

plottime <- read.table("plottime.dat", header=F)$V1

name <- numeric(0)
for (i in 1:count)
{
	name[i] <- as.expression(bquote(.(paste(plottime[i])) ~ x ~ 't'[cc]))
}

legend(Mx-0.2*(Mx-mx), My, legend=name, pch=19, col=colors, bty="n", cex=0.8)

dev.off()
