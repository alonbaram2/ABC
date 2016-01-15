close all


totalTrials = 400;
stateBlockDuration = totalTrials/2;
toGenerate = 30000; % number of schedule to generate for each combination of highProb and pChange
pCands = [0.8, 0.5, 0.2];
AllCues = nan(toGenerate,totalTrials);
AllOutcomes = nan(toGenerate,totalTrials);
AllVol = nan(toGenerate,totalTrials);
AllProbs = nan(toGenerate,totalTrials);

AllEvidA = nan(toGenerate,totalTrials);
AllEvidB = nan(toGenerate,totalTrials);
AllEvidC = nan(toGenerate,totalTrials);
AllEvidAB = nan(toGenerate,totalTrials);
AllEvidAC = nan(toGenerate,totalTrials);
AllEvidSepAC = nan(toGenerate,totalTrials);
AllEvidSepAB = nan(toGenerate,totalTrials);
AllEvidAflipA = nan(toGenerate,totalTrials);
AllEvidABflipA = nan(toGenerate,totalTrials);
AllEvidACflipA = nan(toGenerate,totalTrials);
tic;
for i=1:toGenerate
    
    [cues, outcomes, vol, probs] = generateDataBlocks (totalTrials,pCands,cuesArray);
    
    AllCues(i,:)=cues;
    AllOutcomes(i,:)=outcomes;
    AllVol(i,:)=vol;
    AllProbs(i,:)=probs;
    
    indA = find(cues==0);
    indB = find(cues==1);
    indC = find(cues==2);
    indAC = find (cues==0 | cues==2);
    indAB = find (cues==0 | cues==1);
    
    probsFlipA = probs;
    probsFlipA(cues==0) = 1 - probsFlipA(cues==0);
    outcomesFlipA = outcomes;
    outcomesFlipA(cues==0) = outcomesFlipA(cues==0).*(-1) + 1;
    
    %% run different models
    
    [evidAC,prior_p_qHAC,~,~] = Bayes(probs(cues==0 | cues==2),outcomes(cues==0 | cues==2),vol(cues==0 | cues==2),indAC);
    [evidAB,prior_p_qHAB,~,~] = Bayes(probs(cues==0 | cues==1),outcomes(cues==0 | cues==1),vol(cues==0 | cues==1),indAB);
    [evidA,prior_p_qHA,~,~] = Bayes(probs(cues==0),outcomes(cues==0),vol(cues==0),indA);
    [evidB,prior_p_qHB,~,~] = Bayes(probs(cues==1),outcomes(cues==1),vol(cues==1),indB);
    [evidC,prior_p_qHC,~,~] = Bayes(probs(cues==2),outcomes(cues==2),vol(cues==2),indC);
    [evidAflipA,prior_p_qHAflipA,~,~] = Bayes(probsFlipA(cues==0),outcomesFlipA(cues==0),vol(cues==0),indA);
    [evidACflipA,prior_p_qHACflipA,~,~] = Bayes(probsFlipA(cues==0 | cues==2),outcomesFlipA(cues==0 | cues==2),vol(cues==0 | cues==2),indAC);
    [evidABflipA,prior_p_qHABflipA,~,~] = Bayes(probsFlipA(cues==0 | cues==1),outcomesFlipA(cues==0 | cues==1),vol(cues==0 | cues==1),indAB);
    
    AllEvidA(i,indA) = evidA;
    AllEvidB(i,indB) = evidB;
    AllEvidC(i,indC) = evidC;
    AllEvidAB(i,indAB) = evidAB;
    AllEvidAC(i,indAC) = evidAC;
    AllEvidAflipA(i,indA) = evidAflipA;
    AllEvidABflipA(i,indAB) = evidABflipA;
    AllEvidACflipA(i,indAC) = evidACflipA;
    
    
    AllEvidSepAC(i,indA)=evidA;
    AllEvidSepAC(i,indC)=evidC;
    AllEvidSepAB(i,indA)=evidA;
    AllEvidSepAB(i,indB)=evidB;
    
end
timeSpent = toc;
evids = struct('A',AllEvidA,'B',AllEvidB,'C',AllEvidC,'AB',AllEvidAB,'AC',AllEvidAC,'sepAC',AllEvidSepAC,'sepAB',AllEvidSepAB,'AflipA',AllEvidAflipA,'ABflipA',AllEvidABflipA,'ACflipA',AllEvidACflipA);
schedules = struct('cues',AllCues,'probs',AllProbs,'vol',AllVol,'outcomes',AllOutcomes);
schedulesParams = struct('toGenerate',toGenerate,'totalTrials',totalTrials,'stateBlockDuration',stateBlockDuration,'pCands',pCands);

clear All* evidA* evidB evidC evidAB evidAC ind* blockChange* prior* MeanDiff* cues outcomes* probs* vol i j k flip* evidSep* q_can* H_can* deltaP* totalTrials stateBlock*
save('blocks.mat','evids','schedules','schedulesParams')
