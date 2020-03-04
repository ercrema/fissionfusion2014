%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       DISCUSSION  AND  CONCLUSIONS             %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
\FloatBarrier
\pagebreak
\chapter{ABM Code}
\renewcommand{\thesection}{C.\arabic{section}}

\section{Disturbance-free Model}

All versions of the agent-based models are written in R statistical computing language and can be executed once all scripts are sourced. The main function {\ttfamily FF()} will require a series of sub-functions ({\ttfamily initialise()}, {\ttfamily evaluateFitness()}, {\ttfamily repDeath()}, {\ttfamily fissionfusion()}) for its execution. 

\subsubsection{FF()}

\scriptsize{
  \begin{verbatim}
  #
  # *disturbance free model*
  # 
  #
  # MODEL PARAMETERS:
  # ini       ... Initial Number of Agents
  # P         ... Square root of the number of patches 
  # K         ... Resource Input Size
  # mu        ... Basic Payoff
  # b         ... group benefit
  # sigma     ... payoff uncertainty
  # c1        ... cost 1 for FissionFusion
  # c2        ... cost 2 for FissionFusion
  # c3        ... cost 3 for FissionFusion
  # c4        ... cost 4 for FissionFusion ///New Parameter///
  # timesteps ... number of timesteps
  # k         ... Number of sampled neighbour agents for the model Biased Transmission
  # z         ... Frequency of Decision Making (transmission rate)
  # h         ... Fission Distance (in Chebyshev distance)
  # s         ... Observation Distance (in Chebyshev distance)
  # omega1     ... Mortality Parameter 1
  # omega2      ... Mortality Parameter 2
  # r         ... Reproduction Rate
  # run       ... run number (to be used only for HPC)
  #
  #
  
  FF <- function(ini = 10, P = 10, K = 200, mu = 10, c1 = 3, c2 = 3, c3 = 0, c4 = 0,
                 k = 1, b = 0.5, sigma = 1, timesteps = 300, r = 0.05,
                 omega1 = 1.0, omega2 = 5, size = P, h = 1 , s = 1,
                 z = 1, run = 1)
  {
    # The following line prints the run number for HPC
    print(paste("RunNumber=",run,sep=""));
    
    
    #Initialise model
    tmp <- initialise(ini = ini, P = P, K = K); #initialise model
    groupSpace <- tmp$groupSpace; #extract groupSpace
    
    Raw <- cbind(groupSpace$R, groupSpace$C);
    RawMat <- matrix(0, nrow=length(groupSpace$R), ncol=timesteps);
    Raw <- cbind(Raw, RawMat);
    
    
    for (t in seq(timesteps)){
      
      # STEP1 Fitness Evaluation (computed by group):
      groupSpace <- evaluateFitness(groupSpace, mu = mu, b = b, sigma = sigma);   
      
      # STEP2 Reproduction & Death:
      groupSpace <- repDeath(groupSpace = groupSpace, mu = mu, r = r, omega1 = omega1,
                             omega2 = omega2);
      
      
      #Loophole in case of extinction: 
      if(sum(groupSpace$groupSize) == 0)
      {
        print("extinction!");
        return(Raw);
        break();
      }
      
      
      
      #STEP 3 FissionFusion
      groupSpacePre<-groupSpace
      groupSpace <- fissionfusion(groupSpace, k=k, c1=c1, 
                                  c2=c2, c3=c3, c4=c4, P=P, h=h, s=s, z=z, mu=mu);
      groupSpaceAfter<-groupSpace
      if (any(groupSpace$groupSize<0)){break()}
      
      
      #STEP 4 Record group size distribution
      Raw[,t+2] <- groupSpace$groupSize
      
      
      #RETURN ARGUMENTS    
      return(Raw)
      
    }
    \end{verbatim}
  }
  
  
  
  \subsubsection{initialise()}
  
  \scriptsize{
    \begin{verbatim}
    # Initialise function
    # Reads P (...square root of the Patch number), ini (... the initial number of agents),
    # and K (...the resource input size)
    # Ouputs:
    #
    # agentSet ... a data.frame with the number of rows corresponding to "ini"
    #              containing the following columns:
    #              R ... row coordinate
    #              C ... column coordinate
    #              fitness ... initial fitness (set to 0)
    #              contribution ... initial contribution (set to 0)
    #              groupID ... linker to specific groups
    
    # groupSpace...a data.frame with row number equal to P^2  with the following columns:
    #              R             ... row coordinate
    #              C             ... column coordinate
    #              occupied      ... 1=occupied; 0=not occupied 
    #              preoccupied   ... 1=previously occupied; 0=previously not occupied
    #              groupSize     ... current groupSize
    #              pregroupSize  ... previous groupSIze
    #              T             ... Total Group Contribution
    #              K             ... Resource Input Size
    #              fit           ... Individual Fitness
    
    # world    ...matrix of P by P representing the world
    
    
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
      
      #Define Group rank
      groupSpace$preoccupied = groupSpace$occupied;
      groupSpace$pregroupSize = groupSpace$groupSize;
      return(list(agentSet = agentSet, groupSpace = groupSpace, world = world))
    }
    
    \end{verbatim}
  }
  
  \subsubsection{evaluateFitness()}
  
  \scriptsize{
    \begin{verbatim}
    # evaluateFitness function
    # input: groupSpace, mu, b, sigma
    # outputs:groupSpace (updated)
    
    evaluateFitness <- function(groupSpace, mu, b, sigma)
    {
      index<-which(groupSpace$occupied==1);
      
      
      for (i in index)
      {
        g <- groupSpace[i, ]$groupSize; #collect group size
        groupSpace[i, ]$T =  sum(rnorm(n = g, mean = mu+(g^b)-1, sd = sigma)); 
        #compute group contribution
        
        if (groupSpace[i, ]$T>groupSpace[i, ]$K){
          groupSpace[i,]$T=groupSpace[i,]$K}
        #In case of overexploitation use K instead of T
        
        groupSpace[i, ]$fit = groupSpace[i,]$T/g;
        #compute individual fitness
      }
      
      return(groupSpace)
    }
    \end{verbatim}
  }
  
  \subsubsection{repDeath()}
  
  \scriptsize{
    \begin{verbatim}
    # repDeath function
    # Inputs: groupSpace, mu, r, omega1 and omega2
    # Exports:  groupSpace (updated)
    #
    #
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
    
    \end{verbatim}
  }
  
  \subsubsection{fissionfusion()}
  \scriptsize{
    \begin{verbatim}
    # inputs ...groupSpace, k,c1,c2,c3,s,h,P,z
    # outputs.... updated groupSpace
    # reverse dependencies ...FF()
    
    fissionfusion <- function(groupSpace, k, c1, c2, c3, c4, s, h, P, z, mu)
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
      
      
      
      #Create AgentSet with the following columns:
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
                if((myF>=modelF) & (myF>(mu-c1)))
                {
                  #STAY, DO NOTHING
                  agentSet[x, ]$moved <- 1
                }
                
                if (myF <= (mu-c1) & emptyPatches & ((myF>=modelF) | (modelF<=(mu-c1))))
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
                
                if ((myF <= (modelF-c2) | (myF <= (mu-c1))) & modelF > (mu-c1)){
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
                if ((myF < (modelF-c1) | (myF<=(mu-c1))) & emptyPatches)
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
                if ((myF > (modelF-c1) | (myF>(mu-c1))))
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
                if (myF <= (modelF-c3)){
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
                if (myF<(mu-c4) & modelF<(mu-c4) & agentSet[modelID, ]$moved == 0){  
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
            
            if (myF<(mu-c1) & emptyPatches){
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
    
    \end{verbatim}
  }
  
  
  \section{Predator-prey model}
  \normalsize{The predator-prey version of the model is based on a slightly modified version of the main function {\ttfamily FF()}, along with a modified version of {\ttfamily evaluateFitness()} and the addition of the new sub-function {\ttfamily regenerateResources()}.}
  
  \subsubsection{FF2()}
  
  \scriptsize{
    \begin{verbatim}
    #
    # *Predator-prey model model*
    # 
    #
    # Additional MODEL PARAMETERS:
    #  Kini		... Initial prey population size
    #  KMax		... Prey carrying capacity (kappa)
    #  gR		... Prey population growth rate (zeta)
    #  beta		... Prey population reslience
    
    
    FF2 <- function(ini = 10, P = 10, K = 200, mu = 10, c1 = 3, c2 = 3, c3 = 0, c4 = 0,
                    k = 1, b = 0.5, sigma = 1, timesteps = 300, r = 0.05,
                    omega1 = 1.0, omega2 = 5, size = P, h = 1 , s = 1,
                    z = 1, run = 1, Kini=200, KMax=200, gR=2, beta=0.3)
    {
      # The following line prints the run number for HPC
      print(paste("RunNumber=",run,sep=""));
      
      
      #Initialise model
      tmp <- initialise(ini = ini, P = P, K = K); #initialise model
      groupSpace <- tmp$groupSpace; #extract groupSpace
      
      Raw <- cbind(groupSpace$R, groupSpace$C);
      RawMat <- matrix(0, nrow=length(groupSpace$R), ncol=timesteps);
      Raw <- cbind(Raw, RawMat);
      
      
      for (t in seq(timesteps)){
        
        # STEP1 Fitness Evaluation (computed by group):
        groupSpace <- evaluateFitness2(groupSpace, mu = mu, b = b, sigma = sigma, beta=beta);   
        
        # STEP2 Reproduction & Death:
        groupSpace <- repDeath(groupSpace = groupSpace, mu = mu, r = r, omega1 = omega1,
                               omega2 = omega2);
        
        # STEP3 Resource regeneration:
        groupSpace <- regenerateResources(groupSpace=groupSpace, KMax=KMax, gR=gR);
        
        
        #Loophole in case of extinction: 
        if(sum(groupSpace$groupSize) == 0)
        {
          print("extinction!");
          return(Raw);
          break();
        }
        
        
        
        #STEP 4 FissionFusion
        groupSpacePre<-groupSpace
        groupSpace <- fissionfusion(groupSpace, k=k, c1=c1, 
                                    c2=c2, c3=c3, c4=c4, P=P, h=h, s=s, z=z, mu=mu);
        groupSpaceAfter<-groupSpace
        if (any(groupSpace$groupSize<0)){break()}
        
        
        #STEP 5 Record group size distribution
        Raw[,t+2] <- groupSpace$groupSize
        
        
        #RETURN ARGUMENTS    
        return(Raw)
        
      }
      \end{verbatim}
    }
    
    \subsubsection{evaluateFitness2()}
    
    \scriptsize{
      \begin{verbatim}
      # evaluateFitness function
      # input: groupSpace, mu, b, sigma, baseline
      # output: groupSpace (updated)
      #
      # 
      
      evaluateFitness <- function(groupSpace, mu, b, sigma, beta)
      {
        index<-which(groupSpace$occupied == 1);
        
        
        for (i in index)
        {
          startK <- groupSpace[i, ]$K;
          g <- groupSpace[i, ]$groupSize; #collect group size
          groupSpace[i, ]$T =  sum(rnorm(g, mean = mu+(g-1)^b, sd = sigma)); 
          #compute group contribution
          
          if (groupSpace[i, ]$T > (startK-startK*beta)){
            groupSpace[i,]$T <- startK-startK*beta}
          #In case of overexploitation use (K-K* beta) instead of T
          
          groupSpace[i, ]$fit = groupSpace[i,]$T/g
          #compute individual fitness
          
        }
        
        return(groupSpace)
      }
      
      \end{verbatim}
      
      
      \subsubsection{regenerateResources()}
      \scriptsize{
        \begin{verbatim}
        # regenerateResources function
        # input: groupSpace, KMax, gR
        # output: groupSpace (updated)
        
        regenerateResources <- function (groupSpace, KMax, gR)
        {
          diff <- (groupSpace$K - groupSpace$T);
          groupSpace$K <- diff + (diff * gR * (1- diff/KMax));
          groupSpace$K[which(groupSpace$K<0)]=1 ;
          #exit strategy in case there is complete depletion 
          #this should never happen, provided that beta>0
          groupSpace$T <- 0 ;
          return(groupSpace)  
        }
        \end{verbatim}
        
        
        
        
        
        \section{Exogenic Disturbance Model}
        \normalsize{The exogenic disturbance model is also based on a slightly modified version of the main function {\ttfamily FF()}.}
        
        
        \subsubsection{FF3()}
        
        \scriptsize{
          \begin{verbatim}
          #
          # *disturbance free model*
          # 
          #
          # Additional MODEL PARAMETERS:
          # Kseq	...Vector (with length timesteps) representing the change of K over time. 
          #
          #
          
          FF <- function(ini = 10, P = 10, K = 200, mu = 10, c1 = 3, c2 = 3, c3 = 0, c4 = 0,
                         k = 1, b = 0.5, sigma = 1, timesteps = 300, r = 0.05,
                         omega1 = 1.0, omega2 = 5, size = P, h = 1 , s = 1,
                         z = 1, run = 1,  Kseq=c(rep(200,299),
                                                 seq(from=200,t=100,length=5),rep(100,196)))
          {
            # The following line prints the run number for HPC
            print(paste("RunNumber=",run,sep=""));
            
            
            #Initialise model
            tmp <- initialise(ini = ini, P = P, K = K); #initialise model
            groupSpace <- tmp$groupSpace; #extract groupSpace
            
            Raw <- cbind(groupSpace$R, groupSpace$C);
            RawMat <- matrix(0, nrow=length(groupSpace$R), ncol=timesteps);
            Raw <- cbind(Raw, RawMat);
            
            
            for (t in seq(timesteps)){
              
              # STEP1 Environment Change
              groupSpace$K<-Kseq[t]
              
              # STEP2 Fitness Evaluation (computed by group):
              groupSpace <- evaluateFitness(groupSpace, mu = mu, b = b, sigma = sigma);   
              
              # STEP3 Reproduction & Death:
              groupSpace <- repDeath(groupSpace = groupSpace, mu = mu, r = r, omega1 = omega1,
                                     omega2 = omega2);
              
              
              #Loophole in case of extinction: 
              if(sum(groupSpace$groupSize) == 0)
              {
                print("extinction!");
                return(Raw);
                break();
              }
              
              #STEP 4 FissionFusion
              groupSpacePre<-groupSpace
              groupSpace <- fissionfusion(groupSpace, k=k, c1=c1, 
                                          c2=c2, c3=c3, c4=c4, P=P, h=h, s=s, z=z, mu=mu);
              groupSpaceAfter<-groupSpace
              if (any(groupSpace$groupSize<0)){break()}
              
              
              #STEP 5 Record group size distribution
              Raw[,t+2] <- groupSpace$groupSize
              
              
              #RETURN ARGUMENTS    
              return(Raw)
              
            }
            \end{verbatim}
          }
          
          
          