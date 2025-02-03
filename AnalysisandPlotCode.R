library(ape)
library(INLA)
library(BayesX)
library(spdep)
library(tidyverse)
library(forstringr)

load("NigeriaFatalityGit.Rdata")
DDat = NigeriaFatality[["Data"]]

# Get the binary response
DDat$FatBinom <- ifelse(DDat$FATALITIES>0,1,0)



###### Calculate Moran I statistic

# Get distance matrix
dists <- as.matrix(dist(cbind(DDat$LONGITUDE, DDat$LATITUDE)))
dists.inv <- 1/(1+dists)
diag(dists.inv) <- 0

Moran.I(log(DDat$FATALITIES+1),dists.inv )

###### Fit Binary component Model Using descrete spatiotemporal model

# Get map of nigeria

map = NigeriaFatality[["MapofNigeria"]]
map_shp=bnd2sp(map)

# convert to inla graph
temp <- poly2nb(map_shp)
nb2INLA("LDN.graph", temp)
graph=inla.read.graph(filename="LDN.graph")

# Linear predictor
formulabinom = FatBinom ~ 1 + Season +factor(SubEventID)+ f(YEAR,model = "ar1") + 
  f(stateID, model = "besag", graph = graph,group = YEAR2, control.group = list(model = "ar1"))

# Estimation
Model.st_binary.f <- inla(formulabinom,
                          data =  DDat,
                          family = "binomial",
                          control.predictor = list(compute = TRUE, link = 1),
                          control.compute = list(dic = TRUE, waic=TRUE,cpo=TRUE),verbose = FALSE)

# Summary of result
summary(Model.st_binary.f)

###### Fit count component Model Using descrete spatiotemporal model

# Offset 
DDat$Popul2 =DDat$Popul/max(DDat$Popul)
DDat$LogPopul =log(DDat$Popul)

# Linear predictor
formulacount = FATALITIES ~ 1+ LogPopul + Season + factor(SubEventID)+  f(YEAR,model = "ar1") + 
  f(stateID, model = "besag", graph = graph,group = YEAR2, control.group = list(model = "ar1"))


# Estimation

Model.st_count <- inla(formulacount,
                       data =  DDat, 
                       family = "zeroinflatedpoisson1",
                       control.predictor = list(compute = TRUE, link = 1),
                       control.compute = list(dic = TRUE, waic=TRUE,cpo=TRUE),verbose = FALSE)


# Summary
summary(Model.st_count)

##### Plot probability 

state2 = NigeriaFatality[["StateAndId"]]$state2
cloc2 = NigeriaFatality[["StateAndId"]]$cloc2
shp.map = NigeriaFatality[["ShpMap"]]

# Predicted probabilities

predProb = Model.st_binary.f$summary.fitted.values$mean
DDat = cbind(DDat,predProb)
prob = DDat %>% group_by(stateID) %>% summarise(mean(predProb))

estMean=prob$`mean(predProb)`

## Plot
plt = list()
SHP.MAP=NULL
for (k in 1:(length(estMean))) {
  
  w=1:37
  #
  estMean.aux=estMean[w]
 # estlw.aux=estlw[w]
 # estup.aux=estup[w]
  
  est_data <- data.frame(region= state2,
                         est=estMean[cloc2]
  )
  
  shp.map$count=0
  
  s <- function(region,est){
    shp.map$count[which(shp.map$STATE==as.character(region))] <- est
    return(shp.map)
  }
  
  # Mean 
  
  for (i in 1:37) {
    shp.map <-s(est_data[i,1],est_data[i,2]) 
  }

  SHP.MAP = rbind(SHP.MAP,shp.map)
}

AB = st_coordinates(st_centroid(SHP.MAP$geometry[1:nrow(SHP.MAP)])) %>% as.data.frame()
SHP.MAP$LONGITUDE = AB$X
SHP.MAP$LATITUDE = AB$Y

# Plot with ggplot2

plt=SHP.MAP %>% ggplot()+
  geom_sf( aes(fill =count,group = STATE)) +
  theme_bw()+scale_fill_viridis_c(option = "D", direction = 1,name="", guide=guide_colorbar(
    barheight = unit(40, units = "mm"),
    barwidth = unit(1, units = "mm"),
  ))+labs(x="",y="",title ="")




#######################
# SPDE/ continuous spatial model
########################

# Extract coordinate

coords = cbind(DDat$LONGITUDE ,DDat$LATITUDE)

# Create mesh

mesh <- inla.mesh.2d(coords, max.edge = c(2, 2),cutoff = .1) 
plot(mesh)
 
# Year id
k <- DDat$YEAR %>% unique() %>% length()

DDat$YEAR2 <- ((DDat$YEAR)-1996)
k=DDat$YEAR2  %>% unique()%>% length()

# Random noise

DDat$iddd = 1:length(DDat$YEAR2)

# SPDE model
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            prior.range = c(1.38, 0.95), # P(range < 1.38) = 0.95
                            prior.sigma = c(10, 0.95))   # P(sigma > 10) = 0.0.95

# SPDE index
iseti <- inla.spde.make.index('i', n.spde = spde$n.spde,
                              n.group = k)

# derive projection matrix A

A <- inla.spde.make.A(mesh = mesh,
                      loc = cbind(DDat$LONGITUDE, DDat$LATITUDE), group = DDat$YEAR2) 



# # Model fitting for the Zero inlfated model for the count component for the continuous model

stk <- inla.stack(
  data = list(Y = cbind(DDat$FATALITIES)),
  A = list(A, 1),
  effects = list(iseti,
                 data.frame(Intercept.zero = 1,
                            Season = DDat$Season,
                            SubEventID = DDat$SubEventID,
                            YEARz=DDat$YEAR,
                            YEAR2z=DDat$YEAR2#,
                            #  iddz=DDat$iddd
                 )),
  tag = "Zero")


# Linear predictor
h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))

formulae <- Y ~ -1 +Season+ factor(SubEventID)+ f(YEAR2z, model="ar1")+ f(i, model = spde, group = i.group, 
    control.group = list(model = 'ar1', hyper = h.spec))

# Estimation
Model.st_count_spde  <- inla(formulae,  
                             family = c("zeroinflatedpoisson1"),
                             data = inla.stack.data(stk), 
                             control.predictor = list(compute = TRUE,
                                                      A = inla.stack.A(stk)), 
                             control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
                             ,verbose = TRUE)

######### PLot observed and oredicted density
b= Model.st_count_spde$summary.fitted.values$mean[1:15867]
data.frame(true = stk$data$data$Y,
           pred= Model.st_count_spde$summary.fitted.values$mean[1:15867],
           pred0 =b  ) %>% ggplot()+
  geom_density(aes(x=true),color = "magenta",alpha=0.2,size=1)+
  geom_density(aes(x=pred),color="orange",alpha=0.2,size=1,linetype=2)+
  theme_bw()

########

## Model fitting for the binary component for the continuous model

# Data stack

stkBin <- inla.stack(
  data = list(Y = cbind(DDat$FatBinom)),
  A = list(A, 1),
  effects = list(iseti,
                 data.frame(Intercept.zero = 1,
                            Season = DDat$Season,
                            SubEventID = DDat$SubEventID,
                            YEARz=DDat$YEAR,
                            YEAR2z=DDat$YEAR2#,
                            #iddz=DDat$iddd
                 )),
  tag = "binary")

# Linear predictor

formulae2 <- Y ~ -1 +Season+ factor(SubEventID)+f(YEAR2,model="ar1")+
  f(i, model = spde, group = i.group, 
    control.group = list(model = 'ar1', hyper = h.spec))
 
# Estimation

Model.st_binary_spde  <- inla(formulae2,  
                              family = c("binomial"),
                              data = inla.stack.data(stkBin), 
                              control.predictor = list(compute = TRUE,
                                                       A = inla.stack.A(stkBin)), 
                              control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE)
                              ,verbose = TRUE)


############## 
#SPDE plot
############## 

### Projection of the spatial effect 

# Get shape file
shp.mapState =  NigeriaFatality[["ShpMap"]]
# Get boundary
boundary=NigeriaFatality[["Boundary"]]
stepsize <- 0.08 

# Construct projection grid
nxy <- round(
  c(diff(range( boundary[, 1])), 
    diff(range( boundary[, 2]))) / stepsize)
projgrid <- inla.mesh.projector(
  mesh, xlim = range(boundary[, 1]), 
  ylim = range(boundary[, 2]), dims = nxy)

xmean <- list()

#### projection of the Binary component

Model.st_binary = Model.st_binary_spde
for (j in 1:k){
  xmean[[j]] <- inla.mesh.project(
    projgrid, Model.st_binary$summary.random$i$mean[iseti$i.group == j])
}


df <-  expand.grid(x = projgrid$x, y = projgrid$y)

datresult= NULL
datresult2 = NULL
X = 1996
for (i in 1:k) {
  
  df$mean_s <- as.vector(xmean[[i]])
  
  ind <- point.in.polygon(
    df$x, df$y,
    boundary[, 1], boundary[, 2]
  )
  
  dff <- df[which(ind == 1), ]
  dff$id=i+X
  datresult = rbind(datresult,dff)
}
datresult$idd  <- datresult$id
datresult$idcol <- ifelse(datresult$mean_s<0.0,1,1)


subshp =  st_coordinates(shp.mapState) %>%as.data.frame()

##### plot with ggplot2

p <- ggplot() +  
  geom_tile(data = datresult, aes(x=x, y=y,fill = mean_s), alpha=datresult$idcol) +
  scale_fill_viridis_c(option = "H", direction = 1,name="", guide=guide_colorbar(
    barheight = unit(60, units = "mm"),
    barwidth = unit(1.5, units = "mm")
    #title.position = "top",
    #title.vjust=0.5,
    #label.vjust = 0.5
    #label.position = "left",
    #direction = "vertical"
  ))+
  labs(title ="",
       y = "",x="") +
  theme_light()+scale_size(guide = guide_legend(direction = "vertical"))+
  geom_sf(data=shp.map,color="grey",size=0,shape = ".",fill=NA,alpha=0.05,lwd = 0.1)+
  geom_sf(data = shp.mapState,color="black",size=0,shape = ".",fill=NA,alpha=.4)

p+facet_wrap(~id,#scales = "free",
             ncol = 5)+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black'))+
  theme(legend.position = "right",
        strip.text = element_text(size = 9,
                                  margin = margin()))

ggsave("NigeriabinarySPDE.pdf", width =10.23, height = 10.23, units = "in")

