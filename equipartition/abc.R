# This program calculates time in nbody units of last snapshot, 1% and 99% quantile in mass, the number of particles in last snapshot, and the time of core collapse.
# This should be run from the folders with a folder with this structure:
#
# ls ~/Desktop/10times/
#10  12  14  16  18  20  22  24  5  7  9      abc.Rout  histm.Rout  merge.Rout
#11  13  15  17  19  21  23  25  6  8  abc.R  histm.R   merge.R     tcc.dat
# and each folder inside is like this:
#1  10  2  3  4  5  6  7  8  9  all  scripts  test_t.sh  this
#this script is used only on  simulations with 10 different runs (nbot on those with only one run)

library(data.table)

# % radius in mass
rpercent <- function(r, m, p)
{
	M <- sum(m)
	rm <- data.table(radius=r, mass=m)
	rsort <- rm[order(radius),]
	ma <- cumsum(rsort[,mass])
	ra <- rsort[,radius]
	rad <- ra[which.min(abs(ma - p*M))]
	return(rad)
}

r5 <- numeric(0)
# core collapse time: the time when r5 is minimum
findtcc <- function(datfile) 
{
	gcdat <- read.table(datfile, header = F) #read the file
	
	m <- gcdat$V7
	r <- sqrt( (gcdat$V1-median(gcdat$V1))**2 + (gcdat$V2-median(gcdat$V2))**2 + (gcdat$V3-median(gcdat$V3))**2 )
#	r2D <- sqrt( (gcdat$V2-median(gcdat$V2))**2 + (gcdat$V3-median(gcdat$V3))**2 )

	rad_5 <- rpercent(r, m, 0.05)
	r5 <- log10(rad_5)

	return(r5)
}

folders <- system("ls | grep -v abc", intern = T)
folders

c0000s <- character(0)
clast <- character(0)
nlast <- numeric(0)
tcc <- numeric(0)


for (f in c(10, 14:25))
{
	for (i in 1:10)
	{
		command <- paste("ls ", f, "/", i, "/c_0000.dat", sep = "")
		command2 <- paste("ls ", f, "/", i, "/c_*.dat | tail -n 1", sep = "")
		command3 <- paste("head -n 1 `", command2, "` | cut -b 18-", sep = "")
		command4 <- paste("wc `", command2, "` | cut -b -7", sep = "")
		command5 <- paste("ls ", f, "/", i, "/c_*.dat", sep = "")
		command6 <- paste("ls ", f, "/", i, "/time.dat", sep="")
		
		c0000s <- c(c0000s, system(command, intern = T))
		clast <- c(clast, system(command3, intern = T))
		nlast <- c(nlast, as.numeric(system(command4, intern = T)) - 2)

		cfiles <- system(command5, intern=T)
		time <- system(command6, intern=T)
		
		ttable <- read.table(time, header=F)
		r5 <- sapply(cfiles, findtcc)

		tcc <- c(tcc, ttable[which.min(r5),1])
	
	}
}

tcc

masses1 <- numeric(0)
masses99 <- numeric(0)

for(c0000 in c0000s)
{
	m <- read.table(c0000, header = F)$V7
	masses1 <- c(masses1, quantile(m, 0.01))
	masses99 <- c(masses99, quantile(m, 0.99))
}


data.frame(clast, masses1, masses99, nlast, tcc)
