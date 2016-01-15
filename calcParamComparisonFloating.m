function [scores] = calcParamComparisonFloating (evids, schedules, schedulesParams)


totalTrials = schedulesParams.totalTrials;
stateBlockDuration = schedulesParams.stateBlockDuration;
deltaPCand = schedulesParams.deltaPCand;
toGenerate = schedulesParams.toGenerate;

window = 5:5:15;

meanDiffEvidAC1=nan(length(deltaPCand),toGenerate);
meanDiffEvidAC2=nan(length(deltaPCand),toGenerate);
meanDiffEvidAB1=nan(length(deltaPCand),toGenerate);
meanDiffEvidAB2=nan(length(deltaPCand),toGenerate);

meanAC1=nan(length(deltaPCand),toGenerate);
meanAC2=nan(length(deltaPCand),toGenerate);
meanAB1=nan(length(deltaPCand),toGenerate);
meanAB2=nan(length(deltaPCand),toGenerate);

runAvgACindA = nan(length(deltaPCand),toGenerate,length(window),totalTrials);
runAvgACindC = nan(length(deltaPCand),toGenerate,length(window),totalTrials);
runAvgACindAC = nan(length(deltaPCand),toGenerate,length(window),totalTrials);

runAvgABindA = nan(length(deltaPCand),toGenerate,length(window),totalTrials);
runAvgABindB = nan(length(deltaPCand),toGenerate,length(window),totalTrials);
runAvgABindAB = nan(length(deltaPCand),toGenerate,length(window),totalTrials);

runAvgACflipAindA = nan(length(deltaPCand),toGenerate,length(window),totalTrials);
runAvgACflipAindC = nan(length(deltaPCand),toGenerate,length(window),totalTrials);
runAvgACflipAindAC = nan(length(deltaPCand),toGenerate,length(window),totalTrials);

runAvgABflipAindA = nan(length(deltaPCand),toGenerate,length(window),totalTrials);
runAvgABflipAindB = nan(length(deltaPCand),toGenerate,length(window),totalTrials);
runAvgABflipAindAB = nan(length(deltaPCand),toGenerate,length(window),totalTrials);

cumsumACindA = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumACindC = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumACindAC = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumABindA = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumABindB = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumABindAB = nan(length(deltaPCand),toGenerate,totalTrials);

cumsumACflipAindA = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumACflipAindC = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumACflipAindAC = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumABflipAindA = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumABflipAindB = nan(length(deltaPCand),toGenerate,totalTrials);
cumsumABflipAindAB = nan(length(deltaPCand),toGenerate,totalTrials);


for i=1:length(deltaPCand)
    for j=1:toGenerate
        cues = squeeze(schedules.cues(i,j,:));
        outcomes = squeeze(schedules.outcomes(i,j,:));
        vol = squeeze(schedules.vol(i,j,:));
        probs = squeeze(schedules.probs(i,j,:)); 
        
        evidA = squeeze(evids.A(i,j,:));
        evidB = squeeze(evids.B(i,j,:));
        evidC = squeeze(evids.C(i,j,:));        
        evidAC = squeeze(evids.AC(i,j,:));
        evidAB = squeeze(evids.AB(i,j,:));
        evidSepAC = squeeze(evids.sepAC(i,j,:));
        evidSepAB = squeeze(evids.sepAB(i,j,:));
        evidACflipA = squeeze(evids.ACflipA(i,j,:));
        evidABflipA = squeeze(evids.ABflipA(i,j,:));
        
        indA = find(cues==0);
        indB = find(cues==1);
        indC = find(cues==2);
        indAC = find (cues==0 | cues==2);
        indAB = find (cues==0 | cues==1);
        
        meanDiffEvidAC1(i,j) = max([mean(evidAC(indAC(find(indAC<=stateBlockDuration)))),mean(evidACflipA(indAC(find(indAC<=stateBlockDuration))))]) - mean(evidSepAC(indAC(find(indAC<=stateBlockDuration))));
        meanDiffEvidAC2(i,j) = max([mean(evidAC(indAC(find(indAC>stateBlockDuration)))),mean(evidACflipA(indAC(find(indAC>stateBlockDuration))))]) - mean(evidSepAC(indAC(find(indAC>stateBlockDuration))));
        meanDiffEvidAB1(i,j) = max([mean(evidAB(indAB(find(indAB<=stateBlockDuration)))),mean(evidABflipA(indAB(find(indAB<=stateBlockDuration))))]) - mean(evidSepAB(indAB(find(indAB<=stateBlockDuration))));
        meanDiffEvidAB2(i,j) = max([mean(evidAB(indAB(find(indAB>stateBlockDuration)))),mean(evidABflipA(indAB(find(indAB>stateBlockDuration))))]) - mean(evidSepAB(indAB(find(indAB>stateBlockDuration))));

        meanAC1(i,j) = max([mean(log(evidAC(indAC(find(indAC<=stateBlockDuration))))),mean(log(evidACflipA(indAC(find(indAC<=stateBlockDuration)))))]) - mean(log(evidSepAC(indAC(find(indAC<=stateBlockDuration)))));
        meanAC2(i,j) = max([mean(log(evidAC(indAC(find(indAC>stateBlockDuration))))),mean(log(evidACflipA(indAC(find(indAC>stateBlockDuration)))))]) - mean(log(evidSepAC(indAC(find(indAC>stateBlockDuration)))));
        meanAB1(i,j) = max([mean(log(evidAB(indAB(find(indAB<=stateBlockDuration))))),mean(log(evidABflipA(indAB(find(indAB<=stateBlockDuration)))))]) - mean(log(evidSepAB(indAB(find(indAB<=stateBlockDuration)))));
        meanAB2(i,j) = max([mean(log(evidAB(indAB(find(indAB>stateBlockDuration))))),mean(log(evidABflipA(indAB(find(indAB>stateBlockDuration)))))]) - mean(log(evidSepAB(indAB(find(indAB>stateBlockDuration)))));
        
        for r = 1:length(window)
            
            % note that the comparison is between 
            runAvgABindA(i,j,r,indA) = slidefun(@mean,window(r),log(evidAB(indA)) - log(evidA(indA)),'backward');
            runAvgABflipAindA(i,j,r,indA) = slidefun(@mean,window(r),log(evidABflipA(indA)) - log(evidA(indA)),'backward');
            runAvgABindB(i,j,r,indB) = slidefun(@mean,window(r),log(evidAB(indB)) - log(evidB(indB)),'backward');
            runAvgABflipAindB(i,j,r,indB) = slidefun(@mean,window(r),log(evidABflipA(indB)) - log(evidB(indB)),'backward');
            runAvgABindAB(i,j,r,indAB) = slidefun(@mean,window(r),log(evidAB(indAB)) - log(evidSepAB(indAB)),'backward');
            runAvgABflipAindAB(i,j,r,indAB) = slidefun(@mean,window(r),log(evidABflipA(indAB)) - log(evidSepAB(indAB)),'backward');
            runAvgACindA(i,j,r,indA) = slidefun(@mean,window(r),log(evidAC(indA)) - log(evidA(indA)),'backward');
            runAvgACflipAindA(i,j,r,indA) = slidefun(@mean,window(r),log(evidACflipA(indA)) - log(evidA(indA)),'backward');
            runAvgACindC(i,j,r,indC) = slidefun(@mean,window(r),log(evidAC(indC)) - log(evidC(indC)),'backward');
            runAvgACflipAindC(i,j,r,indC) = slidefun(@mean,window(r),log(evidACflipA(indC)) - log(evidC(indC)),'backward');
            runAvgACindAC(i,j,r,indAC) = slidefun(@mean,window(r),log(evidAC(indAC)) - log(evidSepAC(indAC)),'backward');
            runAvgACflipAindAC(i,j,r,indAC) = slidefun(@mean,window(r),log(evidACflipA(indAC)) - log(evidSepAC(indAC)),'backward');
        end
        cumsumACindA(i,j,indA) = cumsum(log(evidAC(indA)) - log(evidA(indA)));
        cumsumACindC(i,j,indC) = cumsum(log(evidAC(indC)) - log(evidC(indC)));
        cumsumACindAC(i,j,indAC) = cumsum(log(evidAC(indAC)) - log(evidSepAC(indAC)));
        cumsumABindA(i,j,indA) = cumsum(log(evidAB(indA)) - log(evidA(indA)));
        cumsumABindB(i,j,indB) = cumsum(log(evidAB(indB)) - log(evidB(indB)));
        cumsumABindAB(i,j,indAB) = cumsum(log(evidAB(indAB)) - log(evidSepAB(indAB)));
        
        cumsumACflipAindA(i,j,indA) = cumsum(log(evidACflipA(indA)) - log(evidA(indA)));
        cumsumACflipAindC(i,j,indC) = cumsum(log(evidACflipA(indC)) - log(evidC(indC)));
        cumsumACflipAindAC(i,j,indAC) = cumsum(log(evidACflipA(indAC)) - log(evidSepAC(indAC)));
        cumsumABflipAindA(i,j,indA) = cumsum(log(evidABflipA(indA)) - log(evidA(indA)));
        cumsumABflipAindB(i,j,indB) = cumsum(log(evidABflipA(indB)) - log(evidB(indB)));
        cumsumABflipAindAB(i,j,indAB) = cumsum(log(evidABflipA(indAB)) - log(evidSepAB(indAB)));
                
    end
end

scores = struct('meanAC1',meanDiffEvidAC1,'meanAC2',meanDiffEvidAC2,'meanDiffEvidsAB1',meanDiffEvidAB1,'meanDiffEvidsAB2',meanDiffEvidAB2,...
    'runAvgACindA',runAvgACindA,'runAvgACindC',runAvgACindC,'runAvgACindAC',runAvgACindAC,...
    'runAvgABindA',runAvgABindA,'runAvgABindB',runAvgABindB,'runAvgABindAB',runAvgABindAB,'window',window,...
    'runAvgACflipAindA',runAvgACflipAindA,'runAvgACflipAindC',runAvgACflipAindC,'runAvgACflipAindAC',runAvgACflipAindAC,...
    'runAvgABflipAindA',runAvgABflipAindA,'runAvgABflipAindB',runAvgABflipAindB,'runAvgABflipAindAB',runAvgABflipAindAB,...  
    'cumsumACindA',cumsumACindA,'cumsumACindC',cumsumACindC,'cumsumACindAC',cumsumACindAC,...
    'cumsumABindA',cumsumABindA,'cumsumABindB',cumsumABindB,'cumsumABindAB',cumsumABindAB,...
    'cumsumACflipAindA',cumsumACflipAindA,'cumsumACflipAindC',cumsumACflipAindC,'cumsumACflipAindAC',cumsumACflipAindAC,...
    'cumsumABflipAindA',cumsumABflipAindA,'cumsumABflipAindB',cumsumABflipAindB,'cumsumABflipAindAB',cumsumABflipAindAB);

save('floating.mat','evids','schedules','schedulesParams','scores','-v7.3')
