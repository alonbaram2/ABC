clear all
close all

totalTrials = 600;
stateBlockDuration = totalTrials/2;
highProbCand=0.8:0.05:1.0;
flipVolCand=10:3:70;
flipStabCand=10:3:150;

MeanDiffScoresAC1=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
MeanDiffScoresAC2=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
MeanDiffScoresAB1=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
MeanDiffScoresAB2=nan(length(flipVolCand),length(flipStabCand),length(highProbCand));
AllCues = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllOutcomes = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllVol = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllProbs = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);



AllScoreA = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllScoreB = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllScoreC = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllScoreAB = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);
AllScoreAC = nan(length(flipVolCand),length(flipStabCand),length(highProbCand),totalTrials);




for i=1:length(highProbCand)
    for j=1:length(flipVolCand)
        for k=1:length(flipStabCand)
            
            [cues, outcomes, vol, probs] = generateData (stateBlockDuration, totalTrials, flipVolCand(j), flipStabCand(k), highProbCand(i));
            
            AllCues(j,k,i,:)=cues;
            AllOutcomes(j,k,i,:)=outcomes;
            AllVol(j,k,i,:)=vol;
            AllProbs(j,k,i,:)=probs;
            
            indA = find(cues==0);
            indB = find(cues==1);
            indC = find(cues==2);
            indAC = find (cues==0 | cues==2);
            indAB = find (cues==0 | cues==1);
            
            %% run different models
            
            [scoreAC,prior_p_qHAC,q_candidates,H_candidates] = Bayes(probs(cues==0 | cues==2),outcomes(cues==0 | cues==2),vol(cues==0 | cues==2),indAC);
            [scoreAB,prior_p_qHAB,q_candidates,H_candidates] = Bayes(probs(cues==0 | cues==1),outcomes(cues==0 | cues==1),vol(cues==0 | cues==1),indAB);            
            [scoreA,prior_p_qHA,q_candidates,H_candidates] = Bayes(probs(cues==0),outcomes(cues==0),vol(cues==0),indA);
            [scoreB,prior_p_qHB,q_candidates,H_candidates] = Bayes(probs(cues==1),outcomes(cues==1),vol(cues==1),indB);
            [scoreC,prior_p_qHC,q_candidates,H_candidates] = Bayes(probs(cues==2),outcomes(cues==2),vol(cues==2),indC);
            
            AllScoreA(j,k,i,indA) = scoreA;
            AllScoreB(j,k,i,indB) = scoreB;
            AllScoreC(j,k,i,indC) = scoreC;
            AllScoreAB(j,k,i,indAB) = scoreAB;
            AllScoreAC(j,k,i,indAC) = scoreAC;
            
            
            scoreSepAC(1:totalTrials)=0;
            scoreSepAC(indA)=scoreA;
            scoreSepAC(indC)=scoreC;
            scoreSepAC(scoreSepAC==0)=[];
            
            scoreSepAB(1:totalTrials)=0;
            scoreSepAB(indA)=scoreA;
            scoreSepAB(indB)=scoreB;
            scoreSepAB(scoreSepAB==0)=[];            
            
            
%             
%             figure
%             str = sprintf('highProb=%g,flipVol=%g,stabVol=%g',highProbCand(i),flipVolCand(j),flipStabCand(k));
%             suptitle(str);
%             subplot(3,1,1)
%             hold on
%             plot(indAC,scoreAC,'.')
%             plot(indA,scoreA,'r.')
%             plot(indC,scoreC,'g.')
%             plot(indA,outcomes(cues==0)/5+1,'r.','MarkerSize',8)
%             plot(indC,outcomes(cues==2)/5+1,'c.','MarkerSize',8)
%             plot([stateBlockDuration*0.5 stateBlockDuration*0.5],[0 1],'--k')
%             plot([stateBlockDuration stateBlockDuration],[0 1],'--k',[stateBlockDuration*1.5 stateBlockDuration*1.5],[0 1],'--k')
%             legend('AC','A','C')
%             
%             subplot(3,1,2)
%             hold on
%             plot(indAC,scoreAC-scoreSepAC,'-')
%             plot([stateBlockDuration*0.5 stateBlockDuration*0.5],[0 1],'--k')
%             plot([stateBlockDuration stateBlockDuration],[0 1],'--k',[stateBlockDuration*1.5 stateBlockDuration*1.5],[0 1],'--k')
%             plot([0 totalTrials],[0 0])
%             
%             subplot(3,1,3)
%             hold on
%             plot(indA,probs(cues==0),'r-.')
%             plot(indC,probs(cues==2),'c-.') %indB,probs(cues==1),'b-.^',
%             plot([stateBlockDuration*0.5 stateBlockDuration*0.5],[0 1],'--k',[stateBlockDuration stateBlockDuration],[0 1],'--k',...
%                 [stateBlockDuration*1.5 stateBlockDuration*1.5],[0 1],'--k');
%             legend('A','C')
            
            [~, blockChangeIndAC] = min(abs(indAC-stateBlockDuration));
            MeanDiffScoresAC1(j,k,i) = mean(scoreAC(1:blockChangeIndAC)-scoreSepAC(1:blockChangeIndAC));
            MeanDiffScoresAC2(j,k,i) = mean(scoreAC(blockChangeIndAC:end)-scoreSepAC(blockChangeIndAC:end));
            
            [~, blockChangeIndAB] = min(abs(indAB-stateBlockDuration));
            MeanDiffScoresAB1(j,k,i) = mean(scoreAB(1:blockChangeIndAB)-scoreSepAB(1:blockChangeIndAB));
            MeanDiffScoresAB2(j,k,i) = mean(scoreAB(blockChangeIndAB:end)-scoreSepAB(blockChangeIndAB:end));            
            
            
        end
    end
end


scores = struct('A',AllScoreA,'B',AllScoreB,'C',AllScoreC,'AB',AllScoreAB,'AC',AllScoreAC,...
    'meanDiffAC1',MeanDiffScoresAC1,'meanDiffAC2',MeanDiffScoresAC2,...
    'meanDiffAB1',MeanDiffScoresAB1,'meanDiffAB2',MeanDiffScoresAB2);
schedules = struct('cues',AllCues,'probs',AllProbs,'vol',AllVol,'outcomes',AllOutcomes);
schedulesParams = struct('flipVolCand',flipVolCand,'flipStabCand',flipStabCand,'highProbCand',highProbCand,'totalTrials',totalTrials,'stateBlockDuration',stateBlockDuration);

clear All* scoreA scoreB scoreC scoreAB scoreAC ind* blockChange* prior* MeanDiff* cues outcomes probs vol i j k flip* scoreSep* q_can* H_can* highProbCand totalTrials stateBlock*

save('total600block300.mat','scores','schedules','schedulesParams')