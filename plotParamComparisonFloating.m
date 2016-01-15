function [] = plotParamComparisonFloating(evids, schedules, schedulesParams, scores)

%% 

totalTrials = schedulesParams.totalTrials;
stateBlockDuration = schedulesParams.stateBlockDuration;
deltaPCand = schedulesParams.deltaPCand;
toGenerate = schedulesParams.toGenerate;
window = scores.window;

minRunAvgSwitchIndAC=nan(length(deltaPCand),toGenerate,length(window));
minCumsumSwitchIndAC=nan(length(deltaPCand),toGenerate);

for i=1:length(deltaPCand)
    for j=1:toGenerate        
        minCumsumSwitchIndAC(i,j) = min(find(scores.cumsumACindAC(i,j,stateBlockDuration+1:totalTrials)<0)); 
        for r = 1:length(window)
            minRunAvgSwitchIndAC(i,j,r) = min(find(scores.runAvgACindAC(i,j,r,stateBlockDuration+1:totalTrials)<0));             
        end    
    end
end

%%
% calculate the maximal diff of cumsum(log(evidence)) for each trial number accross all schedules
% 

maxcumsumACindACPerTrial = zeros(1,totalTrials);
maxcumsumACindACPerTrialI = zeros(1,totalTrials);
maxcumsumACindACPerTrialSub = zeros(totalTrials,2);
    
for k = 1:totalTrials
    tmp = scores.cumsumACindAC(:,:,k);
    tmp = tmp(:);
    [maxcumsumACindACPerTrial(k), maxcumsumACindACPerTrialI(k)] = max(tmp);
    [maxcumsumACindACPerTrialSub(k,1) maxcumsumACindACPerTrialSub(k,2)] = ind2sub([length(schedulesParams.deltaPCand),toGenerate],maxcumsumACindACPerTrialI(k));
    
end


