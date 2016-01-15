function [cues, outcomes, vol, probs ] = generateDataBlocks (totalTrials, pCands,cuesArray)


stateBlockDuration = totalTrials/2;
%%
cues = cuesArray(randi(1000),:) ; % 0<=>A, 1<=>B, 2<=>C

%%

outcomes = nan(1,totalTrials);
vol = nan(1,totalTrials);
probs = nan(1,totalTrials);

lCandCorr = [5, 5, 8, 8, 12, 15];
lCandUnCorr = [5, 5, 7, 10];

l1 = lCandCorr(randi(length(lCandCorr)));
l2 = lCandUnCorr(randi(length(lCandUnCorr)));
p1 = pCands(randi(length(pCands)));
p2 = pCands(randi(length(pCands)));
l1Counter = 1;
l2Counter = 1;
for i = 1:totalTrials

	blockNumber = ceil(i/stateBlockDuration); % starts at 1
    if (blockNumber==1) % AC corr 
        if (cues(i)==0 || cues(i) ==2)
            probs(i) = p1;
            l1Counter = l1Counter+1;            
        else
            probs(i) = p2;
            l2Counter = l2Counter+1;            
        end        
    elseif (blockNumber==2) % AB corr 
        if (cues(i)== 0 || cues(i) == 1)
            probs(i) = p1;       
            l1Counter = l1Counter+1;                       
        else
            probs(i) = p2;
            l2Counter = l2Counter+1;                        
        end 
    end
    
    if l1Counter == l1
        l1 = lCandCorr(randi(length(lCandCorr)));
        p1 = pCands(randi(length(pCands)));
        l1Counter = 1;
    end
    if l2Counter == l2
        l2 = lCandUnCorr(randi(length(lCandUnCorr)));
        p2 = pCands(randi(length(pCands)));        
        l2Counter = 1;
    end
    outcomes(i) = (rand()<probs(i));
end
