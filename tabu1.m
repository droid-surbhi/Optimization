function [globalSol] = tabu1(initSol, tabuSize, maxIters, noImprovLimit, neighborSearch, neighborSize, costLimit, objfunc, neighborSearchLimit)
% initSol = initial Solution
% tabuSize = maximum size of tabu list
% maxIters = maximum number of iterations
% noImprovLimit = maximum number of iterations for which if there is no
% improvement in the solution, then the search should terminate
% neighborSearch = minimum difference from selected candidate
% neighborSize = size of neighbor list
% costLimit = minimum cost for which search can terminate
% objfunc = objective function in string
% neighborSearchLimit = maximum difference from selected candidate

% Initialization
%neighborSearchLimit = 3
D = size(initSol, 2);
k = 0;
close all;
count = 0;
objFunc = str2func(objfunc);
globalSol = [initSol objFunc(initSol)];
tabuList = [];
figure(2);
fit_fig= plot([0],[0],'c');
hold on;
costVec= zeros(maxIters,1);
costVec(1) = globalSol(end);
set(fit_fig,'xDataSource','1:count');
set(fit_fig,'yDataSource','costVec(1:count)');
candidateList = [globalSol];


solSet = candidateList;

%global fitness_xend fitness_y fitness_avg fitness_yavg ptss
%fitness_y = []; fitness_xend = []; fitness_avg = []; fitness_yavg = []; ptss = [];

while ((k <=noImprovLimit) && (count <=maxIters)) && (max(candidateList(:,end)) > costLimit)
%define neighbourhood
    candidateList = [];
    neighborCandidate = globalSol(1:end-1);
    neighbor = [];
    for i = 1:D
        if (neighborCandidate(1,i)-neighborSearch >= 0)
            neighbor(:,i) = abs(randi([neighborCandidate(1,i)-neighborSearch, neighborCandidate(1,i)+neighborSearch],neighborSize,1));
        else
            neighbor(:,i) = abs(randi([0, neighborCandidate(1,i)+neighborSearch],neighborSize,1));
        end
    end
                  
      % make candidate list and add their fitness value in third row
      for i = 1:size(neighbor,1)
          sCandidate = neighbor(i,:);
          if size(tabuList) == 0
              check = 0;
          else
              scandidatelist = [];
              for j = 1:D
                  scandidatelist = [scandidatelist, sCandidate(j)*ones(size(tabuList, 1),1)];
              end
            check = sum(tabuList==scandidatelist, 2);
          end
          check(check~=D) = 0;
          if (sum(check) < D)
              sCandidate = [sCandidate, objFunc(sCandidate)];
              candidateList = [candidateList;sCandidate];
              solSet = [solSet;candidateList];
          end
      end
      [minCostVal, index] = min(candidateList(:,end));
% check if best fitness of candidate list > fitness of current best
% solution and update best solution and tabulist
    if(minCostVal < globalSol(:,end))&&(sum(candidateList(index,1:end-1)< 0)==0)
         tabuList = [tabuList;globalSol(1:end-1)];
        
        globalSol = candidateList(index,:);
        if (index == size(candidateList))
            addList = candidateList(1:end-1,1:end-1);
        else addList = [candidateList(1:index-1,1:end-1);candidateList(index+1:end,1:end-1)];
        end
        tabuList = [tabuList;addList];
        if neighborSearch > neighborSearchLimit
            neighborSearch = neighborSearch - 1;
        end
        count = count+1;
        k = 0;
    else
        tabuList = candidateList(:,1:end-1);
        count = count+1;
        k = k+1;
        neighborSearch = neighborSearch + 1;
    end
    
    % manage size of tabu list
    if size(tabuList,1) > tabuSize
        diff = size(tabuList,1) - tabuSize;
        tabuList = tabuList(diff+1:end,:);
    end
    costVec(count)= globalSol(end);
    refreshdata(fit_fig,'caller');
    drawnow;
    display([count, k])
end
%figure(4);
% obj_fun = patch(solSet(:,1),solSet(:,2),solSet(:,3),solSet(:,3),'FaceColor','interp');
