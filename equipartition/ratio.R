# Calculate the ratio between radii which contain 50% and 99% mass of total cluster

library(data.table)

# calculate the number of stars in each snapshots
command <- paste("wc c_*.dat | awk '{print $1}'")
np <- as.numeric(system(command, intern=T))
np <- np[-length(np)]-2
tt <- read.table("time.dat")
# find time (this is also possible)
# command3 <- paste("head -n 1 `", command2, "` | cut -b 18-", sep = "")
table <- cbind(np, tt)
colnames(table) <- c("num", "t")

# calculate the radius which contains p% mass of total cluster
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
r50 <- numeric(0)
r99 <- numeric(0)

# find r5
findr5 <- function(datfile) 
{
	gcdat <- read.table(datfile, header = F) #read the file
	
	m <- gcdat$V7
	r <- sqrt( (gcdat$V1-median(gcdat$V1))**2 + (gcdat$V2-median(gcdat$V2))**2 + (gcdat$V3-median(gcdat$V3))**2 )
	r2D <- sqrt( (gcdat$V2-median(gcdat$V2))**2 + (gcdat$V3-median(gcdat$V3))**2 )

	rad_5 <- rpercent(r, m, 0.05)

	r5 <<- c(r5, rad_5)
}

# find r99/r50
findratio <- function(datfile) 
{
	gcdat <- read.table(datfile, header = F) #read the file
	
	m <- gcdat$V7
	r <- sqrt( (gcdat$V1-median(gcdat$V1))**2 + (gcdat$V2-median(gcdat$V2))**2 + (gcdat$V3-median(gcdat$V3))**2 )
	r2D <- sqrt( (gcdat$V2-median(gcdat$V2))**2 + (gcdat$V3-median(gcdat$V3))**2 )

	rad_50 <- rpercent(r, m, 0.5)
	rad_99  <- rpercent(r, m, 0.99)

	r50 <- rad_50
	r99 <- rad_99

	return(r99/r50)
}

#find r5 for each snapshots
command2 <- paste("ls c_*.dat", sep = "")
cfiles <- system(command2, intern=T)
sapply(cfiles, findr5)
table <- cbind(table, r5)

tcc <- table$t[which.min(table$r5)] 

findratio(cfiles[table$t==tcc])
findratio(cfiles[1])
findratio(cfiles[length(cfiles)])



q()


command4 <- paste("wc `", command2, "` | cut -b -7", sep = "")
command6 <- paste("ls ", f, "/", i, "/time.dat", sep="")

		
c0000s <- c(c0000s, system(command, intern = T))
clast <- c(clast, system(command3, intern = T))
nlast <- c(nlast, as.numeric(system(command4, intern = T)) - 2)

#find the core collapse time
time <- system(command6, intern=T)
		
ttable <- read.table(time, header=F)
radtable <- data.frame(time=ttable, r5=r5) # leftside: column name

tcc <- c(tcc, radtable$time[which.min(radtable$r5)])


masses1 <- numeric(0)
masses99 <- numeric(0)

for(c0000 in c0000s)
{
	m <- read.table(c0000, header = F)$V7
	masses1 <- c(masses1, quantile(m, 0.01))
	masses99 <- c(masses99, quantile(m, 0.99))
}


data.frame(clast, masses1, masses99, nlast, tcc)
