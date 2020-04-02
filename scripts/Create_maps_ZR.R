## Create_maps_ZR.R

## Code for the mapping as described in:
## Reed et al. Mapping the genetic and environmental aetiology of autistic traits in Sweden and the UK (submitted for publication and available on BioRxiv)

## Copyright Zoe Reed 2020
## Distributed under the terms of the GNU General Public License version 3

######################################
## Read in results and prepare data
######################################

load("Results.RData")

## Load in target location data as well here that were used in the previous script

## Combine results and target locations
mapMatrix<-cbind(target_locations, results)
colnames(mapMatrix)<-c("east","north", "A", "C", "E", "V", "2.5%_A", "97.5%_A", "2.5%_C", "97.5%_C", "2.5%_E", "97.5%_E")

######################################
## Mapping preparation
######################################

## The colour scale for the points on the map and histograms (19 different colours)
mapColours<-c("#20D2FF","#23BCEE","#26A7DD","#2A92CC","#2D7DBB","#3168AA","#345399", "#383E88","#3B2977","#3F1466","#541561","#69165D","#7E1859","#941955","#A91A50","#BE1C4C","#D41D48","#E91E44","#FF2040")

## Function to work out which colour each point should be. Has 19 bins for colours. Data is sorted into these bins.
calcColourIndex<-function(data,bins=numeric(19)){
  kolours <- numeric(length(data))
  bins <- bins
  dmin <- min(data)
  dmax <- max(data)
  ivl <- (dmax-dmin)/19
  louuer <- dmin-1
  if(bins[19]==0){
  	for(i in 1:length(bins)){
      bins[i] <- dmin+(ivl*(i))
    }
  }
  for(i in 1:length(bins)){
    for(row in 1:length(kolours)){
      x = data[row];
      if((x > louuer) && (x <= bins[i])){
      	kolours[row] <- i
      } else if(x > bins[19]) {
      	kolours[row]<-19      	
      } else if(x < bins[1]){
      	kolours[row]<-1
      }
    }
    louuer <- bins[i]
  }
  return(list(min=dmin,max=dmax,interval=ivl,bin=bins,index=kolours,data=data))
}


## load rgdal library (install this with install.packages("rgdal") if necessary)
library(rgdal)

## load background map image e.g. a sahpe file here
setwd("~/Maps/")
map_image<-readOGR(dsn=".",layer="map_image")

## Function to winsorise data for plotting purposes if needed i.e. outliers will make the colour distribution to extreme
winsorize <- function(x, pct=0.04) {
	low <- quantile(x, probs=pct, na.rm=T)
	high <- quantile(x, probs=1-pct, na.rm=T)
	idx <- which(x < low)
	x[idx] <- low
	idx <- which(x > high)
	x[idx] <- high
    x
}

## Create a winsorised version of results
wins<-results
wins[,1]<-winsorize(wins[,1])
wins[,2]<-winsorize(wins[,2])
wins[,3]<-winsorize(wins[,3])


######################################
## Create maps
######################################

## Output as a png file
png("map.png", width=9, height=8, units='in', res=300, pointsize=9)
m<-rbind(c(1,2,3), c(4,5,6))
layout(m, widths=c(3,3,3,3,3,3), heights=c(6,2), respect=T)

## For each map for A, C and E, first the background map is plotted, here in a grey colour
## Then the results (winsorised) are plotted at the target locations, where the colour represents the value
## Finally the axes are added, the values here will need to be changed depending on coordinate system used

## A 
plot(map_image, xaxt="n",yaxt="n",xlab="East",ylab="North", col="grey17", border="grey17", cex=2)
points(target_locations,col= mapColours[calcColourIndex(wins[,1])$index],las=1, pch=20, cex=0.5)
axis(1, at=seq(0, 1000000, by=200000), labels=c(0, 2, 4, 6, 8, 10), col="black",col.ticks="black",las=1)
axis(2,at=seq(6000000, 7600000, by=200000), labels=c(6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6), col="black",col.ticks="black",las=1)

## C
plot(map_image, xaxt="n",yaxt="n",xlab="East",ylab="North", col="grey17", border="grey17", cex=2)
points(target_locations,col= mapColours[calcColourIndex(wins[,2])$index], las=1, pch=20, cex=0.5)
axis(1, at=seq(0, 1000000, by=200000), labels=c(0, 2, 4, 6, 8, 10), col="black",col.ticks="black",las=1)
axis(2,at=seq(6000000, 7600000, by=200000), labels=c(6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6), col="black",col.ticks="black",las=1)

## E
plot(map_image, xaxt="n",yaxt="n",xlab="East",ylab="North", col="grey17", border="grey17", cex=2)
points(target_locations,col= mapColours[calcColourIndex(wins[,3])$index], las=1, pch=20, cex=0.5)
axis(1, at=seq(0, 1000000, by=200000), labels=c(0, 2, 4, 6, 8, 10), col="black",col.ticks="black",las=1)
axis(2,at=seq(6000000, 7600000, by=200000), labels=c(6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6), col="black",col.ticks="black",las=1)

## Create histograms of results data
bwidth_A<-(max(results[,1])-min(results[,1]))/19
bwidth_C<-(max(results[,2])-min(results[,2]))/19
bwidth_E<-(max(results[,3])-min(results[,3]))/19

win_width_A<-(max(wins[,1])-min(wins[,1]))/19
win_width_C<-(max(wins[,2])-min(wins[,2]))/19
win_width_E<-(max(wins[,3])-min(wins[,3]))/19

## Colours for histogram of A
res_breaks_A<-seq(min(results[,1]), max(results[,1]), by=bwidth_A) #uses results
win_breaks_A<-seq(min(wins[,1]), max(wins[,1]), by=win_width_A) #uses wins

lower_A<-res_breaks_A
upper_A<-res_breaks_A
lower_A<-lower_A[lower_A<min(win_breaks_A)]
upper_A<-upper_A[upper_A>max(win_breaks_A)]
new_breaks_A<-c(lower_A, win_breaks_A, upper_A)

hist_colours_A<-vector(length=length(new_breaks_A))

for(i in 1:length(new_breaks_A)){
	if(new_breaks_A[i]<min(win_breaks_A)){
		hist_colours_A[i]<-mapColours[1]
	}
	if(new_breaks_A[i]>max(win_breaks_A)){
		hist_colours_A[i]<-mapColours[19]
	}
}

pos_A<-min(which(hist_colours_A=="FALSE"))
hist_colours_A<-hist_colours_A[-pos_A]
first_A<-min(which(hist_colours_A=="FALSE"))
last_A<-max(which(hist_colours_A=="FALSE"))
hist_colours_A[first_A:last_A]<-mapColours

## Colours for histogram of C
res_breaks_C<-seq(min(results[,2]), max(results[,2]), by=bwidth_C) #uses results
win_breaks_C<-seq(min(wins[,2]), max(wins[,2]), by=win_width_C) #uses wins

lower_C<-res_breaks_C
upper_C<-res_breaks_C
lower_C<-lower_C[lower_C<min(win_breaks_C)]
upper_C<-upper_C[upper_C>max(win_breaks_C)]
new_breaks_C<-c(lower_C,win_breaks_C,upper_C)

hist_colours_C<-vector(length=length(new_breaks_C))

for(i in 1:length(new_breaks_C)){
	if(new_breaks_C[i]<min(win_breaks_C)){
		hist_colours_C[i]<-mapColours[1]
	}
	if(new_breaks_C[i]>max(win_breaks_C)){
		hist_colours_C[i]<-mapColours[19]
	}
}

pos_C<-min(which(hist_colours_C=="FALSE"))
hist_colours_C<-hist_colours_C[-pos_C]
first_C<-min(which(hist_colours_C=="FALSE"))
last_C<-max(which(hist_colours_C=="FALSE"))
hist_colours_C[first_C:last_C]<-mapColours

## Colours for histogram of E
res_breaks_E<-seq(min(results[,3]), max(results[,3]), by=bwidth_E) #uses results
win_breaks_E<-seq(min(wins[,3]), max(wins[,3]), by=win_width_E) #uses wins

lower_E<-res_breaks_E
upper_E<-res_breaks_E
lower_E<-lower_E[lower_E<min(win_breaks_E)]
upper_E<-upper_E[upper_E>max(win_breaks_E)]
new_breaks_E<-c(lower_E,win_breaks_E,upper_E)

hist_colours_E<-vector(length=length(new_breaks_E))

for(i in 1:length(new_breaks_E)){
	if(new_breaks_E[i]<min(win_breaks_E)){
		hist_colours_E[i]<-mapColours[1]
	}
	if(new_breaks_E[i]>max(win_breaks_E)){
		hist_colours_E[i]<-mapColours[19]
	}
}

pos_E<-min(which(hist_colours_E=="FALSE"))
hist_colours_E<-hist_colours_E[-pos_E]
first_E<-min(which(hist_colours_E=="FALSE"))
last_E<-max(which(hist_colours_E=="FALSE"))
hist_colours_E[first_E:last_E]<-mapColours


## A
hist(results[,1], main="Histogram of A", xlab="Results for A", breaks=new_breaks_A, col=hist_colours_A, cex=2, las=1)

## C
hist(results[,2], main="Histogram of C", xlab="Results for C", breaks=new_breaks_C, col=hist_colours_C, cex=2, las=1)

## E
hist(results[,3], main="Histogram of E", xlab="Results for E", breaks=new_breaks_E, col=hist_colours_E, cex=2, las=1)


dev.off()





