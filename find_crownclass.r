#By Caroline Farrior (cfarrior@gmail.com); last date updated:180724
		
######################################
find_crownclass = function(d,dV,PA,phi=0.03615016,theta=1.2819275){ 
#a function to estimate crown class (light availability) for an individual in a forest given it's diameter and the diameter of neighbors within a neighborhood
#Note: Be sure not to include dead individuals' diameters.  
#d: the diameter of the individual of interest (mm)
#dV: a vector of the diameters of the neighbors within area PA around the focal individual (in any order) (mm)
#PA: area of the neighborhood that dV covers, in m^2
#phi, theta: a crown area allometry, crown area (in m2) = phi*d^theta, where d is in mm
#this function returns a list with both:
#layer: the number of layers of trees above the focal individual
#dstar: the diameter above which, if trees "obeyed the PPA" within the given neighborhood, trees are in the canopy 
#CAt: a "crown area index", an estimate of the amount of space filled by individuals with diameters larger than the focal individual. (just a bit more information on how sure you are that tree is 
#meaning of "obey the PPA": all of the tallest trees whose crowns fit in the canopy, are in the canopy
######################################
	dV = as.vector(dV)
	dV = dV[dV>0]
	dV = dV[order(dV,decreasing=TRUE)]
	CAv = cumsum(phi*dV^theta)
	CAt = max(CAv[dV>=d])/PA
	if(length(dV[dV>=d])==0) CAt = 0
	dstar = NaN
	if(max(CAv)<PA) dstar = 0
	if(min(CAv)>PA) dstar = dV[1]
	if(length(CAv[CAv==PA])==1) dstar = dV[CAv==PA]
	if(is.na(dstar)) dstar = mean(c(dV[CAv<PA][length(dV[CAv<PA])],dV[CAv>PA][1]))
  
  layer = NaN
  if(dstar==0) layer = 0 
  if(is.na(layer)) layer = floor(CAt)
  
	return(list(layer=layer,dstar=dstar,CAt=CAt))
}#end find_crownclass one tree centric
