#Simulation model of "Dominance of the Suppressed", Farrior et al. 2016. Contact: cfarrior@gmail.com

#To run, create a folder: "c:/usr/Farrior_etal_SleeperDist/"
#enter the following in R (without the #'s): 
#setwd("c:/usr/Farrior_etal_SleeperDist/")
#source("Farrior_etal_SleeperSimulation.r")
#runV = c(filestem="PaperParameters",PA=1000,deltaT=1,Fnot=0.02,dnot=.02,G.c=6.05,G.u=.534,mu.c=0.0194/2,mu.u=0.0347-0.0194/2,mu=0.0194/2,cutT=400,landscape_m2=125000,censuscutT = 5)
#main_cohorts(runV=runV,tag=1,newplot=TRUE,plotting=TRUE)

#notes:
#When plotting = TRUE, the simulation runs much more slowly but you can watch the patch develop and go through a stand clearing disturbances.  
#During the simulation, you will seeing the size distribution of the single patch as it changes through time. The dots are the number of individuals in a small bin. Stars indicate that the bin is empty (zero impossible to plot on log/log scale). The blue line is a smoothing on those dots to show what the zeros do to the size distribution. 
#The bestfit powerlaw (Fig. 1) is also plotted for reference. 
#At the end of the simulation the size distribution at the landscape scale is plotted. This is assembled from taking snapshots of the single patch many years apart. 

#crown area allometry parameters (Individual crown area (in m2) = phi * d(in mm) ^theta)
phi = 0.03615016; theta = 1.2819275

######################################
main_cohorts = function(runV=c(filestem="default",PA=1000,deltaT=1,Fnot=0.02,dnot=0.02,G.c=6.05,G.u=0.534,mu.c=0.0097,mu.u=0.025,mu=0.0097,cutT=400,landscape_m2=125000),tag=1,newplot=FALSE,plotting=TRUE){
#runs the cohorts model
#runV: contains all of the inputs for the simulation, including:
	#filestem: the name to call the output files
#PA: size in m2 of the plot
#deltaT: timestep in years
#Fnot: individuals produced per m2 of sun-lit crown area at dnot (reproduction to seedling); also used for recruitment rate ind/m2/year
#dnot: the initial size of a seedling, mm
#G.c: the diameter growth rate, in mm/year of individuals in sun
#G.u: the diameter growth, in mm/year of individuals in the understory
#mu.c: the mortality rate in years^-1 of individuals in the sun
#mu.u: the mortality rate in years^-1 of individuals in the understory
#mu: the stand clearing disturbance rate
#cutT: the time between recording snapshots of the forest - these used to build a landscape of independent forests
#landscape_m2: the desired size of the final assembled landscape (landscape_m2 = PA*[total run time]/cutT)
#tag: an identifier for the specific run
#newplot: logical, whether to create a new plotting window
#plotting: logical, whether to plot the size distribution during the simulation, when TRUE slows the run speed considerably
######################################
	filestem = runV[names(runV)=="filestem"][[1]]
	
	#save a file log file with the parameters of the run
	if(tag==1) write.table(cbind(date(),t(runV)),paste("main_cohorts_LOG.txt",sep=""),append=TRUE,sep="\t",col.names=filestem=="start",row.names=FALSE)
	
	#pull the parameters out of runV for use.
	PA = as.numeric(runV[names(runV)=="PA"][[1]])
	deltaT = as.numeric(runV[names(runV)=="deltaT"][[1]])
	Fnot = as.numeric(runV[names(runV)=="Fnot"][[1]])
	dnot = as.numeric(runV[names(runV)=="dnot"][[1]])
	G.c = as.numeric(runV[names(runV)=="G.c"][[1]])
	G.u = as.numeric(runV[names(runV)=="G.u"][[1]])
	mu.c = as.numeric(runV[names(runV)=="mu.c"][[1]])
	mu.u = as.numeric(runV[names(runV)=="mu.u"][[1]])
	mu = as.numeric(runV[names(runV)=="mu"][[1]])
	cutT = as.numeric(runV[names(runV)=="cutT"][[1]])
	landscape_m2 = as.numeric(runV[names(runV)=="landscape_m2"][[1]])
	
	#convert desired landscape size into the length of time to run the simulation (depends on the plot area and the frequency of taking snapshots).
	maxT = landscape_m2/PA*cutT

	#used to put a limit on seedlings if needed
	maxF = PA/(phi*dnot^theta) 
	
	data = matrix(1,nrow=1,ncol=3) #main matrix with columns: (1) cohort diameter, (2) #of individuals, (3) crown class
	#initialize the matrix with initial cohort of individuals.
	data[,1] = dnot; data[,2] = min(maxF,round(PA*Fnot*deltaT)); data[,3] = 1
	
	#save initial conditions for use when stand gets wiped out. 
	ndata = data 
	
	if(newplot) {X11(width=4,height=4); par(mfrow=c(1,1),oma=c(1,1,1,1))}
	
	#main loop. Remeber data matrix has columns: (1) diameter (2) # of individual (3) crown class
	#crown class = 1 for canopy; crown class = 2 for understory
	for(t in seq(0,maxT,by=deltaT)){
		#Step 1: Mortality
		#Step 1a: Stand clearing disturbance
		if(runif(1)<mu*deltaT) data = ndata  #If stand-level disturbance, set trees back to ndata
		#Step 1b: background mortality
		for(i in seq(1,dim(data)[1])){
			#canopy individuals die yearly with probability mu.c
			if(data[i,3]==1) data[i,2] = rbinom(1,data[i,2],1-mu.c*deltaT)
			#understory individuals die at a yearly rate mu.u
			if(data[i,3]==2) data[i,2] = rbinom(1,data[i,2],1-mu.u*deltaT)
		}
		data = data[data[,2]>0,,drop=FALSE]
		#Step 2: Growth
		#canopy individuals grow at a rate G.c
		data[data[,3]==1,1] = data[data[,3]==1,1]+(G.c)*deltaT
		#understory individuals grow at a rate G.u
		data[data[,3]==2,1] = data[data[,3]==2,1]+(G.u)*deltaT
		#shouldn't be necessary but needed if one goes to variable growth. 
		data = data[data[,1]>0,,drop=FALSE]

		#Step 3: Reproduce
		#Note here, reproduction is not a function of this patch itself, but it easily could be with: babies = round(min(PA,sum(data[,1]^theta*phi*data[,2]))*Fnot*deltaT)
		babies = round(PA*Fnot*deltaT) #dispersal from far away
		if(babies>0){
			babyMatrix = matrix(1,nrow=1,ncol=3)
			CA = sum(phi*data[,1]^theta*data[,2])
			babyMatrix[,1]=dnot; babyMatrix[,2]=babies
			if(CA>PA) babyMatrix[,3]=2  #this could be wrong, but will be corrected in CCassign below.
			data = rbind(data,babyMatrix)
		}
			
		#Step 4: Assign crown class
		CA = sum(phi*data[,1]^theta*data[,2]) #total crown area of individuals in the plot
		if(CA<=PA) data[,3]=1  #if less than the ground area, no need for CCassign. Everyone is in the canopy. 
		if(CA>PA) data = CCassign(data,PA,deltaT) #if greater than the ground area, go through CCassign
		
		#Step 5: Record and plot
		if(floor(t/cutT)==t/cutT){  #records data once every cutT timesteps.
			if(dim(data)[1]>0) write.table(data,paste(filestem,tag,".txt",sep=""),sep="\t",col.names=t==cutT,row.names=FALSE,append=t!=cutT)
		}
		if(plotting) if(max(data[,1])>10) bdata = SizeDistPlot(data,sizem2=PA,xmax=5000,ymin=1e-4)
	}
	
	#plot the size distribution of the whole landscape
	tdata = read.table(paste(runV[[1]],tag,".txt",sep=""),sep="\t",header=TRUE)
	tbdata = SizeDistPlot(as.matrix(tdata),sizem2=as.numeric(runV[names(runV)=="landscape_m2"][[1]]),main=paste(runV[[1]],tag))

}# end main_cohorts

######################################
CCassign = function(data,PA,deltaT){
# main_cohorts function 
# assigns crown class based on the PPA assumption
# assumes all individuals have the same crown area and height allometries
# assumes that CAtot>PA 
######################################

	#make a vector cacaV where the ith entry is the crown area of the i'th cohort plus all cohorts with individuals of greater diameter.
	CAv = phi*data[,1]^theta
	data = data[order(CAv,decreasing=TRUE),]
	CAv = CAv[order(CAv,decreasing=TRUE)]
	cohortCAv = CAv*data[,2]
	cacaV = cumsum(cohortCAv)

	#pull out individuals that are definitely in the canopy and understory
	und = data[cacaV>PA,,drop=FALSE]
	can = data[cacaV<PA,,drop=FALSE]
	
	#split the first cohort in the understory to fill the leftover open canopy space
	canCA = max(0,sum(phi*can[,1]^theta*can[,2]))
	tosplit = und[1,,drop=FALSE]
	opencan = PA-canCA
	splitind_incan = floor(opencan/(phi*tosplit[1,1]^theta))
	und[1,2] = und[1,2] - splitind_incan
	tosplit[,2] = splitind_incan
	can = rbind(can,tosplit)
	can[,3]=1; und[,3]=2

	#piece the data back together
	data = rbind(can,und)
	#always have a tree in the canopy, even if it's bigger than the plot area 
	if(dim(can)[1]==0) data[1,3]=1 
	
	data = data[data[,2]>0,,drop=FALSE]
	return(data)
}# end CCassign


	
######################################
SizeDistPlot = function(data,logby=.14,win=3,sizem2=NaN,main="",census=5,ploton=TRUE,newplot=TRUE,col=1,xaxt="s",yaxt="s",ymin=1e-4,xmax=3500,binV=NaN,justbdata=FALSE,powerfit=list(alpha=2.1283,xmin=22,C=19.6211)){
#Minimal plotting function used in the main_cohorts simulations to plot the size distribution during the simulation run
#this function slows down the speed of the simulation quite a bit.
#data: either a vector of diameters of individuals in the community OR a matrix with the first column the diameter and the second column the number of individuals
#logby: the bin width in logspace for the diameter bins
#win: the number of bins on either side of the target bin to use to smooth over (win=3, makes a window size of 7 bins actually)
#sizem2: the size of the sampled area in m2
#main: header for the plot
#ploton: whether to plot
#newplot: whether to generate a new plotting window
#col: color to use to plot the data
#census: only needs setting if using BCI data with census <=1 where the binning of measurements was different
#xaxt and yaxt: whether to plot an xaxis and yaxis (for use with multiplot figures) 
#ymin: minimum value of y axis (individuals mm^-1 Ha^-1) for plotting
#xmax: maximum value of x (diameter) for plotting
######################################

	if(is.matrix(data)) dV = get_dV(data)
	if(is.vector(data)) dV = data

	dV = dV[order(dV)]
	dV = round(dV)
	dV = dV[dV>10]
	
	if(length(dV)==0) return(NaN)
	
	#make bin categories to group individuals by
	if(is.na(binV[1])){
		binV = exp(seq(log(min(dV)),log(max(dV))+logby*win,by=logby))
		#fix bin categories for censuses BCI censuses that were binned (census 0,1)
		if(census<=1) binV = c(seq(10,55,by=5),exp(seq(log(60),log(max(dV))+logby*win,by=logby)))
	}
	bdata = NULL
	for(i in seq(1,length(binV)-1)){
		d = dV[dV<binV[i+1]]
		d = d[d>=binV[i]]
		bdata = rbind(bdata,c(exp(mean(c(log(binV[i+1]),log(binV[i])))),length(d)/(binV[i+1]-binV[i])))
	}
	bdata[,2] = bdata[,2]/(sizem2)*10000 #generates individuals/mm/Ha
	if(justbdata) return(bdata)
	
	sbin = bdata[,1]
	if(win==0) sdata = bdata
	if(win>0){
		sdata = bdata[seq(1,win),]
		for(i in seq(win+1,length(sbin)-win)){
			n=NULL
			for(j in seq(i-win,i+win)){
				d = dV[dV<sbin[j+1]]
				d = d[d>=sbin[j]]
				d = d[!is.na(d)]
				n = c(n,length(d)/(sbin[j+1]-sbin[j]))
			}
			sdata = rbind(sdata,c(sbin[i],mean(n))) #this is not the mean in log space, because we need the zeros here
		}
	}	
	
	sdata[,2] = sdata[,2]/(sizem2)*10000 #generates individuals/mm/Ha
	
	if(ploton){
		if(newplot) plot(bdata[,1],bdata[,2],main=main,pch=20,xlab="Diameter (mm)",ylab=expression(paste("Individuals ( ",plain(mm^-1)," ", plain(Ha^-1)," )",sep="")),col=col,cex=1,log="xy",xlim=c(10,xmax),ylim=c(ymin,500),xaxt=xaxt,yaxt=yaxt)
	
		if(!newplot) points(bdata[,1],bdata[,2],pch=20,col=col,cex=1)
		zeros = bdata[bdata[,2]==0,,drop=FALSE]
		points(zeros[,1],exp(log(ymin)+log(2.5))+seq(1,dim(zeros)[1])*0,pch=8,cex=1,col=col)
		
		if(col==1) lines(sdata[,1],sdata[,2],col="blue",lwd=2)
		if(col!=1) lines(sdata[,1],sdata[,2],col=col,lwd=2)
		
		dd = seq(powerfit$xmin,max(dV))
		C = 1/(sum(dd^-powerfit$alpha))*length(dV[dV>powerfit$xmin])/sizem2*10000
		lines(dd,C*dd^-powerfit$alpha,lwd=1)
	}
	return(bdata)	
}#end SizeDistPlot


######################################
get_dV = function(data){
# function to turn data frame of cohorts into a list of diameters of individuals
# data is a matrix with columns : (1) diameter (2) number of individuals
######################################
	dV = NULL
	if(dim(data)[1]>0) for(i in seq(1,dim(data)[1])) dV = c(dV,data[i,1][[1]]+seq(1,data[i,2][[1]])*0)
	return(dV)
}#end get_dV
