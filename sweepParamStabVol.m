clear all
close all

totalTrials = 400;
stateBlockDuration = totalTrials/2;
highProbCand=1%0.8:0.05:1.0;
flipVolCand=[25, 28]%10:3:70;
flipStabCand=[25, 80]%10:3:150;


AllCues = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllOutcomes = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllVol = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllProbs = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);

AllEvidA = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidB = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidC = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidAB = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidAC = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidSepAB = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidSepAC = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidAflipA = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidABflipA = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllEvidACflipA = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);

for i=1:length(highProbCand)
    for j=1:length(flipVolCand)
        for k=1:length(flipStabCand)
            
            [cues, outcomes, vol, probs] = generateDataStabVol (stateBlockDuration, totalTrials, flipVolCand(j), flipStabCand(k), highProbCand(i));
            
            AllCues(j,k,i,:)=cues;
            AllOutcomes(j,k,i,:)=outcomes;
            AllVol(j,k,i,:)=vol;
            AllProbs(j,k,i,:)=probs;
            
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
            
            AllEvidA(j,k,i,indA) = evidA;
            AllEvidB(j,k,i,indB) = evidB;
            AllEvidC(j,k,i,indC) = evidC;
            AllEvidAB(j,k,i,indAB) = evidAB;
            AllEvidAC(j,k,i,indAC) = evidAC;
            AllEvidAflipA(j,k,i,indA) = evidAflipA;
            AllEvidABflipA(j,k,i,indAB) = evidABflipA;
            AllEvidACflipA(j,k,i,indAC) = evidACflipA;
            
            
            AllEvidSepAC(j,k,i,indA)=evidA;
            AllEvidSepAC(j,k,i,indC)=evidC;
            AllEvidSepAB(j,k,i,indA)=evidA;
            AllEvidSepAB(j,k,i,indB)=evidB;                                  
                     
            
            
        end
    end
end


evids = struct('A',AllEvidA,'B',AllEvidB,'C',AllEvidC,'AB',AllEvidAB,'AC',AllEvidAC,'sepAC',AllEvidSepAC,'sepAB',AllEvidSepAB,'AflipA',AllEvidAflipA,'ABflipA',AllEvidABflipA,'ACflipA',AllEvidACflipA);
schedules = struct('cues',AllCues,'probs',AllProbs,'vol',AllVol,'outcomes',AllOutcomes);
schedulesParams = struct('flipVolCand',flipVolCand,'flipStabCand',flipStabCand,'highProbCand',highProbCand,'totalTrials',totalTrials,'stateBlockDuration',stateBlockDuration);

clear All* evidA* evidB evidC evidAB evidAC ind* blockChange* prior* MeanDiff* cues outcomes* probs* vol i j k flip* evidSep* q_can* H_can* highProbCand totalTrials stateBlock*

save('total600block300.mat','evids','schedules','schedulesParams')