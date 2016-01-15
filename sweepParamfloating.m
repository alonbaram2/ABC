clear all
close all


totalTrials = 400;
stateBlockDuration = totalTrials/2;
toGenerate = 2; % number of schedule to generate for each combination of highProb and pChange
deltaPCand = 0.01:0.05:0.1;

AllCues = nan(length(deltaPCand),toGenerate,totalTrials);
AllOutcomes = nan(length(deltaPCand),toGenerate,totalTrials);
AllVol = nan(length(deltaPCand),toGenerate,totalTrials);
AllProbs = nan(length(deltaPCand),toGenerate,totalTrials);

AllEvidA = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidB = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidC = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidAB = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidAC = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidSepAC = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidSepAB = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidAflipA = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidABflipA = nan(length(deltaPCand),toGenerate,totalTrials);
AllEvidACflipA = nan(length(deltaPCand),toGenerate,totalTrials);


for i=1:length(deltaPCand)
    for j=1:toGenerate
        
        [cues, outcomes, vol, probs] = generateDataFloating (totalTrials,deltaPCand(i));
        
        AllCues(i,j,:)=cues;
        AllOutcomes(i,j,:)=outcomes;
        AllVol(i,j,:)=vol;
        AllProbs(i,j,:)=probs;
        
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
        
        [evidAC,prior_p_qHAC,q_candidates,H_candidates] = Bayes(probs(cues==0 | cues==2),outcomes(cues==0 | cues==2),vol(cues==0 | cues==2),indAC);
        [evidAB,prior_p_qHAB,q_candidates,H_candidates] = Bayes(probs(cues==0 | cues==1),outcomes(cues==0 | cues==1),vol(cues==0 | cues==1),indAB);
        [evidA,prior_p_qHA,q_candidates,H_candidates] = Bayes(probs(cues==0),outcomes(cues==0),vol(cues==0),indA);
        [evidB,prior_p_qHB,q_candidates,H_candidates] = Bayes(probs(cues==1),outcomes(cues==1),vol(cues==1),indB);
        [evidC,prior_p_qHC,q_candidates,H_candidates] = Bayes(probs(cues==2),outcomes(cues==2),vol(cues==2),indC);
        [evidAflipA,prior_p_qHAflipA,q_candidates,H_candidates] = Bayes(probsFlipA(cues==0),outcomesFlipA(cues==0),vol(cues==0),indA);
        [evidACflipA,prior_p_qHACflipA,q_candidates,H_candidates] = Bayes(probsFlipA(cues==0 | cues==2),outcomesFlipA(cues==0 | cues==2),vol(cues==0 | cues==2),indAC);
        [evidABflipA,prior_p_qHABflipA,q_candidates,H_candidates] = Bayes(probsFlipA(cues==0 | cues==1),outcomesFlipA(cues==0 | cues==1),vol(cues==0 | cues==1),indAB);       

        AllEvidA(i,j,indA) = evidA;
        AllEvidB(i,j,indB) = evidB;
        AllEvidC(i,j,indC) = evidC;
        AllEvidAB(i,j,indAB) = evidAB;
        AllEvidAC(i,j,indAC) = evidAC;       
        AllEvidAflipA(i,j,indA) = evidAflipA;
        AllEvidABflipA(i,j,indAB) = evidABflipA;
        AllEvidACflipA(i,j,indAC) = evidACflipA;

        
        AllEvidSepAC(i,j,indA)=evidA;
        AllEvidSepAC(i,j,indC)=evidC;
        AllEvidSepAB(i,j,indA)=evidA;
        AllEvidSepAB(i,j,indB)=evidB;        
        
    end
end

evids = struct('A',AllEvidA,'B',AllEvidB,'C',AllEvidC,'AB',AllEvidAB,'AC',AllEvidAC,'sepAC',AllEvidSepAC,'sepAB',AllEvidSepAB,'AflipA',AllEvidAflipA,'ABflipA',AllEvidABflipA,'ACflipA',AllEvidACflipA);
schedules = struct('cues',AllCues,'probs',AllProbs,'vol',AllVol,'outcomes',AllOutcomes);
schedulesParams = struct('deltaPCand',deltaPCand,'toGenerate',toGenerate,'totalTrials',totalTrials,'stateBlockDuration',stateBlockDuration);

clear All* evidA* evidB evidC evidAB evidAC ind* blockChange* prior* MeanDiff* cues outcomes* probs* vol i j k flip* evidSep* q_can* H_can* deltaP* totalTrials stateBlock*
save('floatingPlayground.mat','evids','schedules','schedulesParams')
