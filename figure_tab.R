y1=c(952,1466,1407,1583,1157,928,3621,277)
y2=c(322.31,1105.07,1867.41,2165.04,1855.75,1272.51,2180.43,595.48)
y3=c(266.07,999.57,1877.61,2351.28,2208.34,1659.27,1038.93,989.92)
y4=c(924.82,1569.02,1608.52,1367.51,1029.91,706.57,3621.22,563.41)
y5=c(598.20,1436.88,1868.75,1931.85,1722.57,1368.68,985.39,1478.65)
y6=c(2107.24,1560.22,1155.20,855.32,633.29,468.69,3620.64,963.33)
y7=c(2398.82,1891.34,1493.71,1179.67,931.66,735.79,851.1,2054.38)

x=c(0,1,2,3,4,5,6,7)
plot(x,y1,type='o',xlab="Count",ylab="Frequency",main="Empirical and fitted frequency",pch=1, lty=1, col=1)
lines(x,y2, type="o", pch=2, lty=2, col=2)
lines(x,y3, type="o", pch=3, lty=2, col=3)
lines(x,y4, type="o", pch=4, lty=2, col=4)
lines(x,y5, type="o", pch=5, lty=2, col=5)
lines(x,y6, type="o", pch=6, lty=2, col=6)
lines(x,y7, type="o", pch=7, lty=2, col=7)
legend(0.5,3700,c(
expression(paste(plain("Observed"))),
expression(paste(plain("6IP"))),
expression(paste(plain("P"))),
expression(paste(plain("6IDW"))),
expression(paste(plain("DW"))),
expression(paste(plain("6IGEO"))),
expression(paste(plain("GEO")))
),pch=c(1,2,3,4,5,6,7),lwd=c(2,2,2,2,2,2,2),col=c(1,2,3,4,5,6,7),text.width=1.5)







