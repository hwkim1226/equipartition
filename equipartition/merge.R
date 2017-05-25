# Merge files of simulations with 10 runs at 21 different times (0.0tcc to 2.0tcc by 0.1tcc)
# all이라는 폴더에 파일이 저장되므로 빈 폴더가 먼저 있어야 함.

# read a file (I copied the output of talpha.R to tcc.dat manually.)
tcc <- read.table("tcc.dat", header=T)

# find the mean of tcc (10 snapshots)
meantcc <- function(alpha)
{
	a <- mean(tcc[tcc$alpha==alpha,]$t)
	return(2*round(0.5*a))
}


alphas <- sort(unique(tcc$alpha))
tccs <- sapply(alphas, meantcc)
data.frame(alphas, tccs)

# the function that finds the folder for a given alpha
findf <- function(alpha)
{
	folder <- ""
	if( !(alpha < 1.0 || (alpha > 1.0 && alpha < 1.4) || alpha > 2.5) )
	{
		folder <- as.character(10*alpha)
	}
	else
	{		
		folder <- paste("../isolated/32k",as.character(10*alpha), sep="")		
	}

	return(folder)
}

# checks if a folder is multiple run or single run (time.dat is in the upper level folder in single runs)
# 0.5~0.9, 1.1~1.3, 2.6~3.0 : TRUE
# 1.0, 1.4~2.5 : FALSE
istimethere <- function(alpha)
{
	folder <- findf(alpha)
	command <- paste("ls ", folder, "/time.dat", sep = "")
	ret <- system(command, intern = T)
	if (length(ret) == 0)
	{
		 ret <- F
	}
	else
	{
		ret <- T
	}
	return(ret)
}
sapply(seq(from=0.5, to=3.0, by=0.1), istimethere)

# returns *.dat files inside the folder
listthedat <- function(folder)
{
	dats <- system(paste("ls ", folder, "/c_*.dat", sep = ""), intern = T)
	return(dats)
}

# find the nearest snapshot for given time in nbody units in a given folder
fileatt_folder <- function(folder, t)
{
		timedatpath <- paste(folder, "/time.dat", sep="")
		timedat <- read.table(timedatpath)
		dats <- listthedat(folder)
		matchingt <- which.min(abs(timedat$V1 - t))
		datfilename <- dats[(1:length(timedat$V1))[matchingt]]
	#	first2lines <- system(paste("head -n 2 ", datfilename, sep =""), intern = T)
		return(datfilename)
}

# find the nearest snapshots for given time in nbody units in a given alpha
# If the alpha corresponds to simulation with multiple runs, it finds many files. Otherwise, it just runs fileatt_folder.
filesatt <- function(alpha, t)
{
	folder <- findf(alpha)
	if(istimethere(alpha))
	{
		return(fileatt_folder(folder, t))
	}
	else
	{
		againfile <- function(i)
		{
			newf <- paste(folder, "/", i, sep = "")
			return(fileatt_folder(newf, t))
		}
		thefiles <- sapply(1:10, againfile)
		return(thefiles)
	}
}

# converts t/tcc to t (in nbody units)
tovertcc_to_t <- function(alpha, tovertcc)
{
	return(tovertcc*meantcc(alpha))
}

# find the nearest snapshots for given time in t/tcc in a given alpha
findfilesatt <- function(alpha, tovertcc)
{
	filesatt(alpha, tovertcc_to_t(alpha, tovertcc))
}

# center the snapshots in position & velocity
center <- function(tab0)
{
	tab <- tab0
	tab$V1 <- tab$V1 - median(tab$V1)
	tab$V2 <- tab$V2 - median(tab$V2)
	tab$V3 <- tab$V3 - median(tab$V3)
	tab$V4 <- tab$V4 - median(tab$V4)
	tab$V5 <- tab$V5 - median(tab$V5)
	tab$V6 <- tab$V6 - median(tab$V6)
	return(tab)
}

# put files together (for multiple runs at a given t/tcc)
mergethefilesattovertcc <- function(alpha, tovertcc)
{
	thefiles <- findfilesatt(alpha, tovertcc)
	if(length(thefiles) > 1)
	{
		tab <- read.table(thefiles[1], header = F)
		tab <- center(tab)		
		for (file in thefiles[2:length(thefiles)])
		{
			tab0 <- read.table(file, header = F)
			tab0 <- center(tab0)
			tab <- rbind(tab, tab0)
		}
	}
	else
	{
			tab <- read.table(thefiles[1], header = F)
			tab <- center(tab)
	}
	wheretowrite <- paste(findf(alpha), "/all/", tovertcc, ".dat", sep="")
	write.table(tab, wheretowrite)
}

# run the function one alpha at a time...
# there was a problem running them altogether.

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(0.6, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(0.7, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(0.8, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(0.9, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.0, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.1, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.2, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.3, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.4, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.5, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.6, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.7, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.8, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(1.9, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.1, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.2, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.3, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.4, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.5, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.6, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.7, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.8, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(2.9, time)
}

for(time in seq(from = 0.0, to = 2.0, by = 0.1))
{
	mergethefilesattovertcc(3.0, time)
}

warnings()

