# Makes a plot of tcc as a function of alpha and fits various functions to it
library(data.table)

ta <- read.table("talpha.txt", header = T)
# (I copied the output of abc.R to tcc.txt manually.)
ten <- read.table("tcc.txt", header = F)$V6

alp <- numeric(0)
alp[1:10] <- 1.0
alp[11:20] <- 1.4
alp[21:30] <- 1.5
alp[31:40] <- 1.6
alp[41:50] <- 1.7
alp[51:60] <- 1.8
alp[61:70] <- 1.9
alp[71:80] <- 2.0
alp[81:90] <- 2.1
alp[91:100] <- 2.2
alp[101:110] <- 2.3
alp[111:120] <- 2.4
alp[121:130] <- 2.5

ten <- cbind(ten, alp)
colnames(ten) <- c("ts", "alphas")

ta <- rbind(ta, ten)
rownames(ta) <- NULL

ta


#ta <- ta[-6,]

pdf("ta.pdf")

# trh=590.0
plot(ta$alphas, ta$ts/590.0, pch = 16, xlab = expression(alpha), ylab=expression(t/t[r*h]))

cols <- rainbow(5)

t <- ta$ts/590.0
a <- ta$alphas

#linear fit
fi <- lm(t ~ a)
lines(c(-1, 3), fi$coefficients[1] + fi$coefficients[2]*c(-1,3), col=cols[1])

lt <- log10(ta$ts/590.0)
la <- log10(ta$alphas)
lfi <- lm(lt ~ la)

xx <- seq(from=0, to=3, by=0.1)
yy <- 10**(lfi$coefficients[1])*xx**lfi$coefficients[2]

lines(xx,yy, col=cols[2])

# quadratic fit
pfi <- lm(t ~ a+I(a^2))
pyy <- pfi$coefficients[1]+pfi$coefficients[2]*xx+pfi$coefficients[3]*xx*xx

lines(xx,pyy, col=cols[3])



#linear time log alpha fit
lilog.fit <- lm(t ~ la + 1)

#log time linear alpha fit
logli.fit <- lm(lt ~ a + 1)


#plot the linear time log alpha fit
yy <- lilog.fit$coefficients[1] + lilog.fit$coefficients[2]*log10(xx)
lines(xx, yy, col = cols[4], lwd = 2)

#plot the log time linear alpha fit
yy <- 10.0**(logli.fit$coefficients[1] + logli.fit$coefficients[2]*xx)
lines(xx, yy, col = cols[5], lwd = 2)

legends <- c('linear', 'power-law', '2nd order', 'li time log a', 'log time li a')

legend(2.0, 3.5, pch=19, legend=legends, col=cols)

#comparing the adjusted R squared
names.fit <- c("linear fit", "lilog fit", "logli fit", "loglog fit", "para fit")
adjRs <- c(summary(fi)$adj.r.squared, summary(lilog.fit)$adj.r.squared, summary(logli.fit)$adj.r.squared, summary(lfi)$adj.r.squared, summary(pfi)$adj.r.squared)

#high adjusted R-square -> better fit
data.frame(names.fit, adjRs)

summary(pfi)


dev.off()
