### Plot to show the saprsity patterns in the coefficients

pdf("plots/Coeffsettup.pdf", width=8, height=5)
plot(gamma, col="olivedrab3", cex=3, pch=16,
     xlab="j",
     main="Coefficient value",
     ylab="")
points(b, col="darkorchid", cex=2, pch=18,
       xlab="j",
       main="Coefficient value")
abline(v=8.5, lty=6, lwd=2)
abline(v=14.5, lty=6, lwd=2)
box() 
legend(1.5,-.6,
       legend=c("gamma","mu"),
       col=c("olivedrab3","darkorchid"),
       cex=1.5,
       pch=c(16,18))
dev.off()
