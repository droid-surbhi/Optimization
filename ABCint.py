import numpy as np
#import pdb

# In[2]:

def calculateFitness(fObjV):
    fFitness=np.zeros(fObjV.shape)
    ind=np.nonzero(fObjV>=0)
    if fObjV.shape == ():
        if fObjV >= 0:
            fFitness=1./(fObjV+1)
        #ind=np.nonzero(fObjV<0)
        if fObjV < 0:
            fFitness=1+abs(fObjV)
    else:
        fFitness[ind]=1./(fObjV[ind]+1)
        ind=np.nonzero(fObjV<0)
        fFitness[ind]=1+abs(fObjV[ind])
    return fFitness


# In[14]:

def runABC(objfun,NP, D, lb, ub, runtime, maxCycle, stopCost, init):
    """Control Parameters of ABC algorithm
    NP=The number of colony size (employed bees+onlooker bees)
    D=The number of parameters of the problem to be optimized
    low=lower bounds of the parameters
    up=upper bounds of the parameters
    runtime= number of runs (count of times the complete optimization has to be run)
    """
    FoodNumber = NP/2   #The number of food sources equals the half of the colony size
    limit = NP     #A food source which could not be improved through "limit" trials is abandoned by its employed bee
    #maxCycle = 70    #The number of cycles for foraging {a stopping criteria}
    epch = []
    #objfun=obj
    #ub=np.array(up)  #upper bounds of the parameters
    #lb=np.array(low)  #lower bounds of the parameters
    
    #runtime:Algorithm can be run many times in order to see its robustness
    
    """ 
    %Foods [FoodNumber][D]; /*Foods is the population of food sources. Each row of Foods matrix is a vector holding D parameters to be optimized. The number of rows of Foods matrix equals to the FoodNumber*/
    %ObjVal[FoodNumber];  /*f is a vector holding objective function values associated with food sources */
    %Fitness[FoodNumber]; /*fitness is a vector holding fitness (quality) values associated with food sources*/
    %trial[FoodNumber]; /*trial is a vector holding trial numbers through which solutions can not be improved*/
    %prob[FoodNumber]; /*prob is a vector holding probabilities of food sources (solutions) to be chosen*/
    %solution [D]; /*New solution (neighbour) produced by v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) j is a randomly chosen parameter and k is a randomlu chosen solution different from i*/
    %ObjValSol; /*Objective function value of new solution*/
    %FitnessSol; /*Fitness value of new solution*/
    %neighbour, param2change; /*param2change corrresponds to j, neighbour corresponds to k in equation v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})*/
    %GlobalMin; /*Optimum solution obtained by ABC algorithm*/
    %GlobalParams[D]; /*Parameters of the optimum solution*/
    %GlobalMins[runtime]; /*GlobalMins holds the GlobalMin of each run in multiple runs*/
    """
    GlobalMins=np.zeros([runtime,1])
    GlobalParameters=np.zeros([runtime,2])
    for r in range(0,runtime):
        """
        % /*All food sources are initialized */
        %/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */
        """ 
        #pdb.set_trace()
        Range = np.tile((ub-lb),[FoodNumber, 1])
        Lower = np.tile(lb, [FoodNumber, 1])
        Foods = np.random.random_integers(lb[0], ub[0], [FoodNumber,1])
        for i in range(1,D):
            Foods = np.concatenate([Foods, np.random.random_integers(lb[i], ub[i], [FoodNumber,1])], axis=1)
        Foods[0][:] = init
        ObjVal = 10*np.ones([1, FoodNumber])
        for t in range(0,FoodNumber):
            ObjVal[0][t] = objfun(Foods[t][:])
            if ObjVal[0][t] <= stopCost:
                GlobalMin = ObjVal[0][t]
                GlobalParams = Foods[t][:]
                break
        
        Fitness=calculateFitness(ObjVal)
        #%reset trial counters
        trial = np.zeros([1, FoodNumber])
        #pdb.set_trace()
        #The best food source is memorized
        BestInd=np.nonzero(ObjVal==min(min(ObjVal)))
        BestInd=BestInd[-1][-1]
        GlobalMin=ObjVal[0][BestInd].copy()
        GlobalParams=Foods[BestInd][:].copy()
        
        iters = 1
        
        while((iters!=maxCycle) and (GlobalMin > stopCost)):
            for i in range(0,FoodNumber):
                #The parameter to be changed is determined randomly
                Param2Change= int(np.fix(np.random.rand()*(D)))
            
                #A randomly chosen solution is used in producing a mutant solution of the solution i
                neighbour=int(np.fix(np.random.rand()*(FoodNumber)))
                
                #Randomly selected solution must be different from the solution i
                while(neighbour==i):
                    neighbour=int(np.fix(np.random.rand()*(FoodNumber)))
                
                sol=Foods[i][:].copy()
                
                #v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})
                sol[Param2Change]=Foods[i][Param2Change]+(Foods[i][Param2Change]-Foods[neighbour][Param2Change])*int(np.fix(np.random.rand()-0.5))*2
                #if generated parameter value is out of boundaries, it is shifted onto the boundaries
                ind=np.nonzero(sol<lb)
                sol[ind]=lb[ind].copy()
                ind=np.nonzero(sol>ub)
                sol[ind]=ub[ind].copy()
                
                
                # evaluate new solution
                ObjValSol=objfun(sol)
                if ObjValSol <= stopCost:
                    GlobalMin = ObjValSol
                    GlobalParams = sol
                    break
                FitnessSol=calculateFitness(ObjValSol)
                
                #a greedy selection is applied between the current solution i and its mutant*/
                if (FitnessSol>Fitness[0][i]):  #If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
                    Foods[i,:]=sol.copy()
                    Fitness[0][i]=FitnessSol.copy()
                    ObjVal[0][i]=ObjValSol.copy()
                    trial[0][i]=0
                else:
                    trial[0][i]=trial[0][i]+1  #if the solution i can not be improved, increase its trial counter
                    
            """  %%%%%%%%%%%%%%%%%%%%%%%% CalculateProbabilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %/* A food source is chosen with the probability which is proportioal to its quality*/
            %/*Different schemes can be used to calculate the probability values*/
            %/*For example prob(i)=fitness(i)/sum(fitness)*/
            %/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
            %/*probability values are calculated by using fitness values and normalized by dividing maximum fitness value*/"""
            #pdb.set_trace()
            prob=(0.9*Fitness/max(max(Fitness)))+0.1
            
            #ONLOOKER BEE PHASE
            i = 0
            t = 0
            while(t<FoodNumber-1 and GlobalMin > stopCost):
                if(np.random.rand()<prob[0][i]):
                    t=t+1
                    #The parameter to be changed is determined randomly
                    Param2Change=int(np.fix(np.random.rand()*(D)))
                    #A randomly chosen solution is used in producing a mutant solution of the solution i
                    neighbour=int(np.fix(np.random.rand()*FoodNumber))
                    #Randomly selected solution must be different from the solution i
                    while(neighbour==i):
                        neighbour=int(np.fix(np.random.rand()*FoodNumber))
                    sol=Foods[i,:].copy()
                    #v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij})
                    sol[Param2Change]=Foods[i][Param2Change]+(Foods[i][Param2Change]-Foods[neighbour][Param2Change])*int(np.fix(np.random.rand()-0.5))*2
                    #if generated parameter value is out of boundaries, it is shifted onto the boundaries
                    ind=np.nonzero(sol<lb);
                    sol[ind]=lb[ind].copy()
                    ind=np.nonzero(sol>ub);
                    sol[ind]=ub[ind].copy()
                    # evaluate new solution
                    ObjValSol=objfun(sol)
                    FitnessSol=calculateFitness(ObjValSol)
                    #a greedy selection is applied between the current solution i and its mutant*/
                    if (FitnessSol>Fitness[0][i]):  #If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
                        Foods[i,:]=sol.copy()
                        Fitness[0][i]=FitnessSol.copy()
                        ObjVal[0][i]=ObjValSol.copy()
                        trial[0][i]=0
                    else:
                        trial[0][i]=trial[0][i]+1  #if the solution i can not be improved, increase its trial counter
                    if t < FoodNumber:    
                        if ObjVal[0][t] <= stopCost:
                            GlobalMin = ObjVal[0][t]
                            GlobalParams = Foods[t][:]
                            break
                i = i+1
                if (i==FoodNumber): 
                    i=0
            
            #pdb.set_trace()
            #The best food source is memorized
            BestInd=np.nonzero(ObjVal==min(min(ObjVal)))
            BestInd=BestInd[-1][-1]
            if ObjVal[0][BestInd] < GlobalMin:
                GlobalMin=ObjVal[0][BestInd].copy()
                GlobalParams=Foods[BestInd][:].copy() 
            
            #SCOUT BEE PHASE
            #determine the food sources whose trial counter exceeds the "limit" value. 
#%In Basic ABC, only one scout is allowed to occur in each cycle*/
            ind=np.nonzero(trial==max(max(trial)))
            ind=ind[-1][-1]
            if (trial[0][ind]>limit):
                trial[0][ind]=0
                sol = np.random.random_integers(lb[0], ub[0], [1,1])
                for i in range(1,D):
                    sol = np.concatenate([sol, np.random.random_integers(lb[i], ub[i], [1,1])], axis=1)
                ObjValSol=objfun(sol)

                FitnessSol=calculateFitness(ObjValSol)
                Foods[ind,:]=sol.copy()
                Fitness[0][ind]=FitnessSol.copy()
                ObjVal[0][ind]=ObjValSol.copy()

            print 'iter='+str(iters)+" ObjVal="+str(GlobalMin)+" sol="+str(GlobalParams)
            iters=iters+1
        #pdb.set_trace()
        GlobalMins[r]=GlobalMin.copy()
        GlobalParameters[r][:] = GlobalParams.copy()
    return GlobalMins, GlobalParameters

