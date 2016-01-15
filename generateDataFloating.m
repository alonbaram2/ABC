function [cues, outcomes, vol, probs ] = generateDataFloating (totalTrials, deltaP)

stateBlockDuration = totalTrials/2;
%%
cues = randi(3,1,totalTrials) - 1; % 0<=>A, 1<=>B, 2<=>C
changedFlag = true;
while changedFlag
    indAC = find (cues==0 | cues==2);
    indAB = find (cues==0 | cues==1);
    indBC = find (cues==2 | cues==1);
    for i = 3:length(indAB)
        if cues(indAB(i-2))==cues(indAB(i-1)) & cues(indAB(i))==cues(indAB(i-1))
            if cues(indAB(i-2))==1
                if isempty(find(indAB(i) - indAC > 0))
                    cues(indAB(i))=0;
                    changedFlag=true;
                    break 
                else
                    cues(indAB(i)) = abs(cues(indAC(find(indAB(i) - indAC > 0,1,'last'))) - 2); % change to A or C - the one that didn't apear last.
                    changedFlag=true;
                    break             
                end
            elseif cues(indAB(i-2))==0
                if isempty(find(indAB(i) - indBC > 0))
                    cues(indAB(i))=2;
                    changedFlag=true;
                    break                 
                else
                    if cues(indBC(find(indAB(i) - indBC > 0,1,'last')))==1
                        cues(indAB(i)) = 2;  % change to B or C - the one that didn't apear last.
                        changedFlag=true;
                        break 
                    elseif cues(indBC(find(indAB(i) - indBC > 0,1,'last')))==2
                        cues(indAB(i)) = 1;
                        changedFlag=true;
                        break 
                    end
                end
            end
        end
    end
    if i==length(indAB)
        changedFlag=false;
    end

    indAC = find (cues==0 | cues==2);
    indAB = find (cues==0 | cues==1);
    indBC = find (cues==2 | cues==1);
    for i = 3:length(indAC)
        if cues(indAC(i-2))==cues(indAC(i-1)) & cues(indAC(i))==cues(indAC(i-1))
            if cues(indAC(i-2))==2
                if isempty(find(indAC(i) - indAB > 0))
                    cues(indAB(i))=0;
                    changedFlag=true;
                    break 
                else
                    cues(indAC(i)) = abs(cues(indAB(find(indAC(i) - indAB > 0,1,'last'))) - 1); % change to A or B - the one that didn't apear last.
                    changedFlag=true;
                    break             
                end
            elseif cues(indAC(i-2))==0
                if isempty(find(indAC(i) - indBC > 0))
                    cues(indAC(i))=2;
                    changedFlag=true;
                    break                 
                else
                    if cues(indBC(find(indAC(i) - indBC > 0,1,'last')))==1
                        cues(indAC(i)) = 2;  % change to B or C - the one that didn't apear last.
                        changedFlag=true;
                        break 
                    elseif cues(indBC(find(indAC(i) - indBC > 0,1,'last')))==2
                        cues(indAC(i)) = 1;
                        changedFlag=true;
                        break 
                    end
                end
            end
        end
    end
    if i==length(indAC)
        changedFlag=false;
    end
end

%%

outcomes = nan(1,totalTrials);
vol = nan(1,totalTrials);
probs = nan(1,totalTrials);

p1 = nan(1,totalTrials);
p2 = nan(1,totalTrials);

p1(1) = rand();
p2(1) = rand();

for i = 1:totalTrials
	blockNumber = ceil(i/stateBlockDuration); % starts at 1
    if (blockNumber==1) % AC corr 
        if (cues(i)==0 || cues(i) ==2)
            probs(i) = p1(i);
            if (i < totalTrials)
                signDeltaP = sign(rand()-0.5);
                p1(i+1) = p1(i) + signDeltaP*deltaP;
                if (p1(i+1) >=1 || p1(i+1)<=0)
                        p1(i+1) = p1(i) - signDeltaP*deltaP;
                end
                p2(i+1) = p2(i);
            end
        else
            probs(i) = p2(i);
            if (i < totalTrials)
                signDeltaP = sign(rand()-0.5);
                p2(i+1) = p2(i) + signDeltaP*deltaP;
                if (p2(i+1) >=1 || p2(i+1)<=0)
                        p2(i+1) = p2(i) - signDeltaP*deltaP;
                end                
                p1(i+1) = p1(i);                
            end
        end        
    elseif (blockNumber==2) % AB corr 
        if (cues(i)== 0 || cues(i) == 1)
            probs(i) = p1(i);
            if (i < totalTrials)
                signDeltaP = sign(rand()-0.5);
                p1(i+1) = p1(i) + signDeltaP*deltaP;
                if (p1(i+1) >=1 || p1(i+1)<=0)
                        p1(i+1) = p1(i) - signDeltaP*deltaP;
                end                
                p2(i+1) = p2(i);
            end            
        else
            probs(i) = p2(i);
            if (i < totalTrials)
                signDeltaP = sign(rand()-0.5);
                p2(i+1) = p2(i) + signDeltaP*deltaP;
                if (p2(i+1) >=1 || p2(i+1)<=0)
                        p2(i+1) = p2(i) - signDeltaP*deltaP;
                end                  
                p1(i+1) = p1(i);                
            end            
        end 
    end    
    outcomes(i) = (rand()<probs(i));
end
