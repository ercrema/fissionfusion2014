#'@title Group Fission-Fusion Agent-Based Model
#'@description Simulates group fission-fusion dynamics
#'@param ini Initial Number of Agents. Default is 10.
#'@param P Size of the world (square root of the total number of cells). Default is 10 (i.e. a 10 x 10 world)
#'@param K Carrying capacity at each cell. Default is 200.
#'@param Kseq Time series of K. If provided K changes depending on the timestep. Default is NA.
#'@param mu Basic payoff of the agents.
#'@param c Threshold of evidence.
#'@param k Proportion of sampled neighbour agents for the model Biased Transmission
#'@param b Cooperation derived benefit in payoff.
#'@param sigma Payoff uncertainty.
#'@param timesteps Number of timsteps in the simulation.
#'@param r Reproduction Rate
#'@param omega1 Mortality Parameter 1
#'@param omega2 Mortality Parameter 2
#'@param h Fission Distance (in Chebyshev distance)
#'@param s Observation Distance (in Chebyshev distance)
#'@param z Frequency of decision making. 
#'@return A matrix of group sizes with dimensions P and timesteps + 2. The first two columns contains the coordinates of each cell, and the remaining columns contain the number of agents occupying a given cell at a particular timestep.
#'@references 
#' Crema, E.R., 2014. A simulation model of fission-fusion dynamics and long-term settlement change. Journal of Archaeological Method and Theory 21, 385–404.
#' Crema, E.R., 2015. Modelling Settlement Rank-Size Fluctuations, in: Wurzer, G., Kowarik, K., Reschreiter, H. (Eds.), Agent-Based Modeling and Simulation in Archaeology, Advances in Geographic Information Science. Springer International Publishing, pp. 161–181. https://doi.org/10.1007/978-3-319-00008-4_8
#' @examples
#' set.seed(1)
#' x = FF(timesteps=500,k=0.5)
#' gg = x[,502]
#' sort(gg[which(gg>0)],TRUE)
#' x2 = FF(timesteps=500,k=0.5,s=1,h=1)
#' gg2 = x2[,502]
#' sort(gg2[which(gg2>0)],TRUE)


FF <- function(ini = 10, P = 10, K = 200, Kseq=NA, mu = 10, c = 3,k = 1, b = 0.5, sigma = 1, timesteps = 500, r = 0.05,omega1 = 1.0, omega2 = 5, h = 100 , s = 100, z = 1)
{
  
  #Initialise model:
  
  if (is.na(Kseq)){Kseq<-rep(K,timesteps)} #If Kseq is not provided, create a flat time-series of K
  
  tmp <- initialise(ini = ini, P = P, K = Kseq[1]); #initialise model
  groupSpace <- tmp$groupSpace; #extract groupSpace
  
  Raw <- cbind(groupSpace$R, groupSpace$C);
  RawMat <- matrix(0, nrow=length(groupSpace$R), ncol=timesteps);
  Raw <- cbind(Raw, RawMat);
  
  pb <- txtProgressBar(min = 1, max = timesteps, style = 3) #initialise progress bar
  
  for (t in seq(timesteps)){
    setTxtProgressBar(pb, t) #update progressbar
    
    # STEP0 Environment Change
    groupSpace$K<-Kseq[t]
    
    # STEP1 Fitness Evaluation (computed by group):
    groupSpace <- evaluateFitness(groupSpace, mu = mu, b = b, sigma = sigma);   
    
    # STEP2 Reproduction & Death:
    groupSpace <- repDeath(groupSpace = groupSpace, mu = mu, r = r, omega1 = omega1,
                           omega2 = omega2);
    
    #Loophole in case of extinction: 
    if(sum(groupSpace$groupSize) == 0)
    {
      close(pb)
      print("extinction!");
      return(Raw);
      break();
    }
    
    #STEP 3 FissionFusion
    groupSpacePre<-groupSpace #make a carbon copy of the groupSpace
    groupSpace <- fissionfusion(groupSpace, k=k, c=c, 
                                P=P, h=h, s=s, z=z, mu=mu);
    groupSpaceAfter<-groupSpace
    
    #if (any(groupSpace$groupSize<0)){break()}
    
    
    #STEP 4 Record group size distribution
    Raw[,t+2] <- groupSpace$groupSize
    
    
    #RETURN ARGUMENTS
  }
  close(pb)
  return(Raw)
  
  
}

#'@title Initialise the model
#'@param ini Initial Number of Agents.
#'@param P Size of the world (square root of the total number of cells). 
#'@param K Carrying capacity at each cell. 
#'@return A list containing the follwing objects:
#'\itemize{
#' \item{\code{agentSet}} {A data.frame with ini rows with the following fields:  R (row coordinate); C (column coordinate); fitness (initial fitness, equal to 0); contribution (initial contribution,equal to 0); groupID (unique identifier for linking to specific groups)}.
#' \item{\code{groupSpace}} {A data.frame with P^2 rows with the following fields: R (row coordinate); C (column coordinate); occupied (whether the cell is currently occupied (1) or not (0)); preoccupied (whether the cell was currently occupied (1) or not (0) in the previous timestep); groupSize (current number of agents in the cell); pregroupSize (previous number of agents in the cell); T (total yield of the group), K (resource input size); fit (individual fitness of the agents in the cell).}
#' \item{\code{world}} {A matrix of dimensions P and P representing the world.}  
#' }

initialise <- function(P, ini, K)
{
  #Create World
  world <- matrix(0, P, P);
  #Create Agent Space
  
  agentSet <- data.frame(R=sample(1:P, ini, TRUE), C = sample(1:P, ini, TRUE),
                         fitness = numeric(length=ini), contribution = numeric(length=ini),
                         groupID = numeric(length=ini));
  
  groupSpace <- expand.grid(R=1:P,C=1:P);
  
  groupSpace <- cbind(groupSpace, occupied = rep(0, length = P^2),
                      preoccupied = rep(0, length=P^2), groupSize = numeric(length=P^2),
                      pregroupSize = numeric(length=P^2), T = numeric(length = P^2),
                      K = rep(K, length=P^2),fit = numeric(length = P^2));
  
  
  #Define Agent's group and update groupSpace and agentSet
  for (i in seq(ini))
  {
    agentSet$groupID[i] = which(groupSpace$R == agentSet$R[i]&groupSpace$C == agentSet$C[i]);
    groupSpace[agentSet$groupID[i], ]$groupSize = groupSpace[agentSet$groupID[i], ]$groupSize + 1;
    groupSpace[agentSet$groupID[i], ]$occupied = 1;
  }
  
  #Update groupSpace
  groupSpace$preoccupied = groupSpace$occupied;
  groupSpace$pregroupSize = groupSpace$groupSize;
  
  return(list(agentSet = agentSet, groupSpace = groupSpace, world = world))
}


#'@title Evaluate the agent's fitness
#'@param groupSpace data.frame with group specific info. Created with the initialise() function.
#'@param mu Basic payoff of the agents. 
#'@param b Cooperation derived benefit in payoff.
#'@param sigma Payoff uncertainty. 
#'@return Returns an updated version of groupSpace.

evaluateFitness <- function(groupSpace, mu, b, sigma)
{
  index<-which(groupSpace$occupied==1); #look for all occupied groupSpace
  
  for (i in index)
  {
    g <- groupSpace[i, ]$groupSize; #collect group size
    groupSpace[i, ]$T =  sum(rnorm(n = g, mean = mu+(g-1)^b, sd = sigma)); 
    #compute group contribution
    
    if (groupSpace[i, ]$T>groupSpace[i, ]$K){
      groupSpace[i,]$T=groupSpace[i,]$K}
    #In case of overexploitation use K instead of T
    
    groupSpace[i, ]$fit = groupSpace[i,]$T/g;
    #compute individual fitness
  }
  
  return(groupSpace)
}

#'@title Birth and Death model.
#'@param groupSpace data.frame with group specific info. Created with the initialise() function.
#'@param mu Basic payoff of the agents. 
#'@param r Reproduction Rate
#'@param omega1 Mortality Parameter 1
#'@param omega2 Mortality Parameter 2
#'@return Returns an updated version of groupSpace.

repDeath <- function(groupSpace, mu, r, omega1, omega2)
{
  index <- which(groupSpace$occupied == 1); #retrieve index of occupied patches
  
  for (i in index)
  {
    Fit <- groupSpace[i, ]$fit; #collect fitness
    G <- groupSpace[i, ]$groupSize; #collect groupSize
    births <- 0
    #births
    births <- sum(runif(G)<((Fit/mu)*r))
    #death:
    deathProb <- 1/(1+exp(1)^((omega1*Fit)-omega2)); #probability of death
    deaths <- sum(runif(G)<deathProb); #actual number of death
    
    #update group size
    G <- G+births-deaths 
    
    # IN case of extinction set everything to 0:
    if (G<=0){
      groupSpace[i,]$groupSize <- 0;
      groupSpace[i,]$occupied <- 0;
      groupSpace[i,]$T <- 0;
      groupSpace[i,]$fit <- 0;
    }
    if (G>0){
      groupSpace[i,]$groupSize <- G;
    }
  }
  
  return(groupSpace)
}


#'@title Fission Fusion Routine
#'@param groupSpace data.frame with group specific info. Created with the initialise() function.
#'@param k Proportion of sampled neighbour agents for the model Biased Transmission
#'@param c Threshold of evidence.
#'@param s Observation Distance (in Chebyshev distance)
#'@param h Fission Distance (in Chebyshev distance)
#'@param P Size of the world (square root of the total number of cells). Default is 10 (i.e. a 10 x 10 world)
#'@param z Frequency of decision making. 
#'@return Returns an updated version of groupSpace.

fissionfusion <- function(groupSpace, k, c, s, h, P, z, mu)
{
  
  #utility function for finding neighbours:
  matNeighbour <- function(D, myLoc, size)
  {
    if (size < Inf){
      sizeSeq <- (-size:size)}
    
    if (size == Inf){
      sizeSeq = 1:D}
    
    if (length(sizeSeq) < D){
      
      L <- length(sizeSeq);
      coordinates <- expand.grid(r=sizeSeq, c=sizeSeq);
      rev <- D:1;
      for (x in 1:L^2)
      {
        tmpR <- coordinates[x, 1] + myLoc[1];
        tmpC <- coordinates[x, 2] + myLoc[2];
        if (tmpR <= 0){
          tmpR <- rev[abs(tmpR)+1];
        }
        
        if (tmpC <= 0)
        {
          tmpC <- rev[abs(tmpC)+1];
        }
        if (tmpR > D){
          tmpR <- tmpR-D;
        }
        if (tmpC > D){
          tmpC <- tmpC-D
        }
        coordinates[x, 1] <- tmpR;
        coordinates[x, 2] <- tmpC;
      }
    }
    if (length(sizeSeq) >= D){
      coordinates <- expand.grid(r=1:D, c=1:D);
    }
    return(coordinates)
  }
  
  
  
  #Create an AgentSet with the following columns:
  #id       ...agents' id
  #R        ...row coordinate
  #C        ...column coordinate
  #fitness  ...fitness
  #groupSize...group Size
  #groupID  ...group ID
  #moved    ...boolean (1=decision taken; 0=decision to be taken)
  
  ########################  
  ####CREATE AGENTSET#####  
  ########################
  
  N = sum(groupSpace$groupSize);
  agentSet = data.frame(id = 1:N, R = numeric(length=N), C = numeric(length=N),
                        fitness = numeric(length=N), groupSize = numeric(length=N),
                        groupID = numeric(length=N), moved = rep(1,N));
  
  index <- which(groupSpace$occupied == 1);
  counter = 1;
  for (i in index)
  {
    G <- groupSpace[i,]$groupSize;
    inputR <- counter:(counter+G-1);
    counter <- counter+G;
    agentSet[inputR, ]$R <- groupSpace[i, ]$R;
    agentSet[inputR, ]$C <- groupSpace[i, ]$C;
    agentSet[inputR, ]$fitness <- groupSpace[i, ]$fit;
    agentSet[inputR, ]$groupSize <-G;
    agentSet[inputR, ]$groupID <-as.numeric(rownames(groupSpace[i, ]));
  }
  
  
  ################################
  #####SELECT DECISIONMAKERS######
  ################################
  
  #Randomize Order of Execution
  order <- sample(1:N)
  #Set frequency of execution
  decisionmakers <- order[runif(N)<z]
  if (length(decisionmakers)>0)
  {
    #decisionmakers still need to make their decision (moved=0) w
    #while all the other agents are treated as if they've already
    #made their choices
    agentSet[decisionmakers,]$moved=0
    
    
    
    for (x in decisionmakers)
    {
      
      #####################################    
      ##########LOOK AROUND################    
      #####################################
      
      
      #Look only at other groups
      
      #Spatially within the neighbourhood:
      myLoc <- c(agentSet[x, ]$R, agentSet[x, ]$C);
      destinations <- matNeighbour(D=P, myLoc=myLoc, size=s)
      # the following selects agents from the agentset with the coordinates of destinations 
      # but without the groupID of the focal agent
      others <- which(agentSet$R %in% destinations$r & agentSet$C %in% destinations$c &
                        agentSet$groupID != agentSet[x,]$groupID)
      
      #if other is not empty (this could happen if the group is isolated spatially)
      if (length(others) > 0){
        
        #evaluate empty patches 
        #Problem of synchronisation, as any refers to the current group space
        #however if this is referred to the <spaces> object, then it will lead to 
        #the problem of co-occurence of agents in the same location.
        #for now the emptyPatches will refer to groupSpace no to spaces
        
        #evaluate emptypatches
        tmp <- matNeighbour(D=P, myLoc=myLoc, size=h);
        tmp2 <- which(groupSpace$R%in%tmp$r & groupSpace$C%in%tmp$c);
        emptyPatches <- any(groupSpace[tmp2, ]$occupied == 0);     
        
        
        #####################################    
        ##########CHOOSE MODEL STAGE#########    
        #####################################    
        #reset the value of K (K UPPERCASE is the actual "k" used for sampling agents)
        
        K=k;  
        K<-length(others)*K;
        K<-ceiling(K);
        
        if (K > length(others)){
          K <- length(others);
        }
        
        modelIDs=sample(x=others,size=K)
        
        #Choose the best fit agent among k individuals
        modelID <- modelIDs[which(agentSet[modelIDs, ]$fitness == 
                                    max(agentSet[modelIDs, ]$fitness))[1]];
        modelF <- agentSet[modelID, ]$fitness
        modelG <- agentSet[modelID, ]$groupSize
        modelGID <- agentSet[modelID, ]$groupID
        myF <- agentSet[x, ]$fitness
        myG <- agentSet[x, ]$groupSize
        myGID <- agentSet[x, ]$groupID
        
        #####################################    
        ##########COMPARISON STAGE########   
        #####################################    
        if (agentSet[x, ]$moved != 1)
        {
          #CASE 1: G vs G#
          
          if (myG>1 & modelG>1)
          {
            if((myF>=modelF) & (myF>(mu-c)))
            {
              #STAY, DO NOTHING
              agentSet[x, ]$moved <- 1
            }
            
            if (myF <= (mu-c) & emptyPatches & ((myF>=modelF) | (modelF<=(mu-c))))
            {
              #EMERGENCY FISSION
              #Reduce former group Size
              groupSpace[myGID, ]$groupSize <- groupSpace[myGID, ]$groupSize-1;
              #set occupied to 0 if there were no more agents
              if (groupSpace[myGID, ]$groupSize == 0) {
                groupSpace[myGID,  ]$occupied <- 0;
              }
              
              #Create a new group
              
              if(length(which(groupSpace$occupied==0))>1)
              {newPlace <- sample(x=which(groupSpace$occupied == 0),size=1);}
              if(length(which(groupSpace$occupied==0))==1)
              {newPlace <- which(groupSpace$occupied == 0)}
              if (groupSpace[newPlace, ]$groupSize > 0) {print("ERROR line 159")}
              groupSpace[newPlace, ]$groupSize <- 1
              groupSpace[newPlace, ]$occupied <- 1
              agentSet[x, ]$moved <- 1
            }
            
            if ((myF <= (modelF-c) | (myF <= (mu-c))) & modelF > (mu-c)){
              #GUIDED MIGRATION
              #Change Group Sizes
              groupSpace[modelGID, ]$groupSize <- groupSpace[modelGID, ]$groupSize+1;
              agentSet[x, ]$moved <- 1
              #ensure the new group now is occupied
              groupSpace[modelGID, ]$occupied <- 1
              groupSpace[myGID, ]$groupSize <- groupSpace[myGID, ]$groupSize-1;
              if (groupSpace[myGID, ]$groupSize==0) {
                groupSpace[myGID, ]$occupied <- 0;
              }
              
            }
            
          }
          
          #CASE 2: G vs L
          
          if (myG>1 & modelG == 1){
            if(myF >= modelF){
              #STAY, DO NOTHING
              agentSet[x, ]$moved <- 1
              #  print("option4")
            }
            if ((myF < (modelF-c) | (myF<=(mu-c))) & emptyPatches)
            {
              #FISSION & EMERGENCY FISSION
              #Reduce former group Size
              groupSpace[myGID, ]$groupSize <- groupSpace[myGID, ]$groupSize-1;
              #handle local extinction
              if (groupSpace[myGID, ]$groupSize == 0) {
                groupSpace[myGID, ]$occupied <- 0}
              #Create a new group
              if(length(which(groupSpace$occupied==0))>1)
              {newPlace <- sample(x=which(groupSpace$occupied == 0),size=1);}
              if(length(which(groupSpace$occupied==0))==1)
              {newPlace <-which(groupSpace$occupied == 0)}
              groupSpace[newPlace, ]$groupSize <- 1;
              groupSpace[newPlace, ]$occupied <- 1;
              agentSet[x, ]$moved <- 1
              #  print("option5")
            }
            if ((myF > (modelF-c) | (myF>(mu-c))))
            {
              agentSet[x, ]$moved <- 1
            }
            
          }
          
          #CASE 3: L vs G
          
          if (myG == 1 & modelG > 1){
            if(myF >= modelF){
              #STAY, DO NOTHING
              agentSet[x, ]$moved <- 1
            }
            if (myF <= (modelF-c)){
              #FUSION from SINGLE
              groupSpace[modelGID, ]$groupSize <- groupSpace[modelGID, ]$groupSize+1;
              groupSpace[modelGID, ]$occupied <- 1;
              groupSpace[myGID, ]$groupSize <- groupSpace[myGID, ]$groupSize-1;
              groupSpace[myGID, ]$occupied <- 0;
              agentSet[x, ]$moved <- 1
            }
            
          }
          
          #CASE 4: L vs L
          
          if (myG==1&modelG==1)
          {
            if(myF>=mu)
            {
              #STAY, DO NOTHING
              agentSet[x, ]$moved <- 1
            }
            
            #the third condition is to ensure that 
            #the other agent did not make any decision, since this step
            # involves both agent making a decision.      
            if (myF < mu & modelF < mu & agentSet[modelID, ]$moved == 0){  
              #FUSION between SINGLES
              groupSpace[myGID, ]$groupSize <- groupSpace[myGID, ]$groupSize-1;
              groupSpace[myGID, ]$occupied <- 0;
              groupSpace[modelGID, ]$groupSize <- groupSpace[modelGID, ]$groupSize+1;
              groupSpace[modelGID, ]$occupied <- 1;
              agentSet[modelID, ]$moved <- 1;
              agentSet[x, ]$moved <- 1
            }
            
          }
        }
      }
      
      
      # CASE 5 : no other groups, emergency fission is still possible:   
      if(length(others) == 0&agentSet[x, ]$moved != 1){
        myF <- agentSet[x, ]$fitness;
        emptyPatches <- any(groupSpace$occupied == 0);
        myGID <- agentSet[x, ]$groupID;
        
        if (myF<(mu-c) & emptyPatches){
          #EMERGENCY FISSION
          #Reduce former group Size
          groupSpace[myGID, ]$groupSize <- groupSpace[myGID, ]$groupSize-1;
          if (groupSpace[myGID, ]$groupSize == 0){
            groupSpace[myGID, ]$occupied <- 0;}
          
          #Create a new group
          if(length(which(groupSpace$occupied==0))>1)
          {newPlace <- sample(x=which(groupSpace$occupied == 0),size=1);}
          if(length(which(groupSpace$occupied==0))==1)
          {newPlace <-which(groupSpace$occupied== 0)}
          groupSpace[newPlace, ]$groupSize <- 1;
          groupSpace[newPlace, ]$occupied <- 1;
          agentSet[x, ]$moved <- 1
          #print("option10")
        }
        
      }
      
      #After the decision is taken the agent is "moved"    
      agentSet[x, ]$moved <- 1;
      
    }  
  }
  return(groupSpace)
}



