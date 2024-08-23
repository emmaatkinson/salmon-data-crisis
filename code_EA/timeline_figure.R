library(here)
library(scales)

d = read.csv(here("data_EA","Number_populationsAssessed_total.csv"))
c = read.csv(here("data_EA","NPAFC_Catch_Stat-1925-2023.csv"), skip=1)
total_catch = c[c$Country=="Canada" & c$Reporting.Area=="Whole country" & c$Species=="Total", 1:ncol(c)]
tc = as.numeric(gsub(",","",total_catch[1,7:ncol(total_catch)]))

png(here("figures","2024-june-salmon-monitoring-timeline.png"), width=8, height=4, units="in",res=750)

plot(d$Year, d$count, bty="n", xlab="Years", ylab="Number of populations monitored", ylim=c(0,3000),xlim=c(1915,2025), xaxt="n", type="l",col="white")
axis(1, at=seq(1915,2025,10), labels=TRUE)

# polygon(x=c(1866,1867,1867,1866), y=c(0,0,2500,2500), col=alpha("grey",.5), border=NA)
# polygon(x=c(1868,1869,1869,1868), y=c(0,0,2500,2500), col=alpha("grey",.5), border=NA)
polygon(x=c(1930,1935,1935,1930), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(1950,1955,1955,1950), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(1970,1980,1980,1970), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(1985,1986,1986,1985), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(1995,1998,1998,1995), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(1998,1999,1999,1998), y=c(0,0,3000,3000), col=alpha("grey",.9), border=NA)
polygon(x=c(2005,2006,2006,2005), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(2009,2010,2010,2009), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(2012,2013,2013,2012), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(2018,2019,2019,2018), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)
polygon(x=c(2021,2022,2022,2021), y=c(0,0,3000,3000), col=alpha("grey",.5), border=NA)

key.years=c(1934,1954,1975,1985,1995,1998,2005,2009,2012,2018,2021)
lines(d$Year, d$count, lwd=2, col=alpha("navy",0.7))
points(d[d$Year%in%key.years,]$Year,d[d$Year%in%key.years,]$count, pch=20,cex=1.5, col="navy")

par(new=TRUE)
plot(1925:2023, tc, type="l", lwd=2, col=alpha("darkred", 0.7), axes=FALSE, bty="n",xlab="", ylab="", lty=6)

dev.off()
