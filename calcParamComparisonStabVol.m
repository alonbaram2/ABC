function [scores] = calcParamComparisonStabVol (evids, schedules, schedulesParams)


totalTrials = schedulesParams.totalTrials;
stateBlockDuration = schedulesParams.stateBlockDuration;
flipVolCand = schedulesParams.flipVolCand;
flipStabCand = schedulesParams.flipStabCand;
highProbCand = schedulesParams.highProbCand;

meanDiffEvidAC1=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
meanDiffEvidAC2=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
meanDiffEvidAB1=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
meanDiffEvidAB2=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));

meanDiffLogEvidAC1=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
meanDiffLogEvidAC2=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
meanDiffLogEvidAB1=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
meanDiffLogEvidAB2=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
                
for i=1:length(highProbCand)
    for j=1:length(flipVolCand)
        for k=1:length(flipStabCand)
            cues = squeeze(schedules.cues(j,k,i,:));
            outcomes = squeeze(schedules.outcomes(j,k,i,:));
            vol = squeeze(schedules.vol(j,k,i,:));
            probs = squeeze(schedules.probs(j,k,i,:));
            
            evidAC = squeeze(evids.AC(j,k,i,:));
            evidAB = squeeze(evids.AB(j,k,i,:));
            evidSepAC = squeeze(evids.sepAC(j,k,i,:));
            evidSepAB = squeeze(evids.sepAB(j,k,i,:));
            evidACflipA = evids.ACflipA(j,k,i,:);
            evidABflipA = evids.ABflipA(j,k,i,:);
            indA = find(cues==0);
            indB = find(cues==1);
            indC = find(cues==2);
            indAC = find (cues==0 | cues==2);
            indAB = find (cues==0 | cues==1);
            
            meanDiffEvidAC1(j,k,i) = max([mean(evidAC(indAC(find(indAC<=stateBlockDuration)))),mean(evidACflipA(indAC(find(indAC<=stateBlockDuration))))]) - mean(evidSepAC(indAC(find(indAC<=stateBlockDuration))));
            meanDiffEvidAC2(j,k,i) = max([mean(evidAC(indAC(find(indAC>stateBlockDuration)))),mean(evidACflipA(indAC(find(indAC>stateBlockDuration))))]) - mean(evidSepAC(indAC(find(indAC>stateBlockDuration))));
            meanDiffEvidAB1(j,k,i) = max([mean(evidAB(indAB(find(indAB<=stateBlockDuration)))),mean(evidABflipA(indAB(find(indAB<=stateBlockDuration))))]) - mean(evidSepAB(indAB(find(indAB<=stateBlockDuration))));
            meanDiffEvidAB2(j,k,i) = max([mean(evidAB(indAB(find(indAB>stateBlockDuration)))),mean(evidABflipA(indAB(find(indAB>stateBlockDuration))))]) - mean(evidSepAB(indAB(find(indAB>stateBlockDuration))));

            meanDiffLogEvidAC1(j,k,i) = max([mean(log(evidAC(indAC(find(indAC<=stateBlockDuration))))),mean(log(evidACflipA(indAC(find(indAC<=stateBlockDuration)))))]) - mean(log(evidSepAC(indAC(find(indAC<=stateBlockDuration)))));
            meanDiffLogEvidAC2(j,k,i) = max([mean(log(evidAC(indAC(find(indAC>stateBlockDuration))))),mean(log(evidACflipA(indAC(find(indAC>stateBlockDuration)))))]) - mean(log(evidSepAC(indAC(find(indAC>stateBlockDuration)))));
            meanDiffLogEvidAB1(j,k,i) = max([mean(log(evidAB(indAB(find(indAB<=stateBlockDuration))))),mean(log(evidABflipA(indAB(find(indAB<=stateBlockDuration)))))]) - mean(log(evidSepAB(indAB(find(indAB<=stateBlockDuration)))));
            meanDiffLogEvidAB2(j,k,i) = max([mean(log(evidAB(indAB(find(indAB>stateBlockDuration))))),mean(log(evidABflipA(indAB(find(indAB>stateBlockDuration)))))]) - mean(log(evidSepAB(indAB(find(indAB>stateBlockDuration)))));
            
        end
    end
end

scores = struct('meanDiffEvidAC1',meanDiffEvidAC1,'meanDiffEvidAC2',meanDiffEvidAC2,'meanDiffEvidAB1',meanDiffEvidAB1,'meanDiffEvidAB2',meanDiffEvidAB2)