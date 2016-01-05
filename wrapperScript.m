close all
highProbToPlot = 0.9;

plotParamComparison(scores, schedules, schedulesParams,highProbToPlot)

flipStab = 40;
flipVol = 28;
highProb=1.0;

%% plot one schedule's estimation and scores

j = find(schedulesParams.flipVolCand==flipVol);
k = find(schedulesParams.flipStabCand==flipStab);
i = find(schedulesParams.highProbCand==highProb);

cues = squeeze(schedules.cues(j,k,i,:));
outcomes = squeeze(schedules.outcomes(j,k,i,:));
vol = squeeze(schedules.vol(j,k,i,:));
probs = squeeze(schedules.probs(j,k,i,:));
totalTrials=schedulesParams.totalTrials;
stateBlockDuration = schedulesParams.stateBlockDuration;

plotBayes=true;

plotModelSelection(cues, outcomes, vol, probs,totalTrials,stateBlockDuration,plotBayes)