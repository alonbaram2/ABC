%%

function plotModelSelection(cues, outcomes, vol, probs,totalTrials,stateBlockDuration,plotBayes)
    
if nargin<5
    totalTrials=length(cues);
    stateBlockDuration=totalTrials/2; % assuming only two state blocks
end
if nargin<7
    plotBayes=0;
end

cues=squeeze(cues); outcomes=squeeze(outcomes); vol=squeeze(vol); probs=squeeze(probs);

indA = find(cues==0); indB = find(cues==1); indC = find(cues==2);
indAC = find (cues==0 | cues==2); indAB = find (cues==0 | cues==1);

% figure
% plot(indA,probs(cues==0),'r-.^',indC,probs(cues==2),'c-.^',... %indB,probs(cues==1),'b-.^',
%     [stateBlockDuration*0.5 stateBlockDuration*0.5],[0 1],'--k',[stateBlockDuration stateBlockDuration],[0 1],'--k',...
%     [stateBlockDuration*1.5 stateBlockDuration*1.5],[0 1],'--k');
% legend('A','C')

%% run different models

[scoreAC,prior_p_qHAC,q_candidates,H_candidates] = Bayes(probs(cues==0 | cues==2),outcomes(cues==0 | cues==2),vol(cues==0 | cues==2),indAC);
[scoreA,prior_p_qHA,q_candidates,H_candidates] = Bayes(probs(cues==0),outcomes(cues==0),vol(cues==0),indA);
[scoreC,prior_p_qHC,q_candidates,H_candidates] = Bayes(probs(cues==2),outcomes(cues==2),vol(cues==2),indC);
[scoreB,prior_p_qHB,q_candidates,H_candidates] = Bayes(probs(cues==1),outcomes(cues==1),vol(cues==1),indB);
[scoreAB,prior_p_qHAB,q_candidates,H_candidates] = Bayes(probs(cues==0 | cues==1),outcomes(cues==0 | cues==1),vol(cues==0 | cues==1),indAB);

if plotBayes
    EqAC=nan(1,length(indAC)); 
    EHAC=nan(1,length(indAC)); 
    EqA=nan(1,length(indA));
    EHA=nan(1,length(indA)); 
    EqC=nan(1,length(indC)); 
    EHC=nan(1,length(indC)); 
    
    EqAB=nan(1,length(indAB)); 
    EHAB=nan(1,length(indAB)); 
    EqB=nan(1,length(indB)); 
    EHB=nan(1,length(indB));
    
    [qq,HH]=ndgrid(q_candidates,H_candidates);
    for i=1:length(outcomes(cues==0 | cues==2))
        EqAC(i)=sum(myvect(prior_p_qHAC(:,:,i).*qq));
        EHAC(i)=sum(myvect(prior_p_qHAC(:,:,i).*HH));
    end
    for i=1:length(outcomes(cues==0))
        EqA(i)=sum(myvect(prior_p_qHA(:,:,i).*qq));
        EHA(i)=sum(myvect(prior_p_qHA(:,:,i).*HH)); 
    end
    for i=1:length(outcomes(cues==2))
        EqC(i)=sum(myvect(prior_p_qHC(:,:,i).*qq));
        EHC(i)=sum(myvect(prior_p_qHC(:,:,i).*HH));       
    end
    for i=1:length(outcomes(cues==0 | cues==1))
        EqAB(i)=sum(myvect(prior_p_qHAB(:,:,i).*qq));
        EHAB(i)=sum(myvect(prior_p_qHAB(:,:,i).*HH));
    end

    for i=1:length(outcomes(cues==1))
        EqB(i)=sum(myvect(prior_p_qHB(:,:,i).*qq));
        EHB(i)=sum(myvect(prior_p_qHB(:,:,i).*HH));       
    end
    
    % plot estimates based on marginal means of distribution - AC
    figure; hold on;   
    title('model estimation - AC')
    plot(indA,probs(cues==0),'r--',indC,probs(cues==2),'c--'); plot(indAC,EqAC,'k'); plot(indAC,vol(cues==0 | cues==2),'k--'); plot(indAC,EHAC,'m');  
    plot(indA,outcomes(cues==0),'r.','MarkerSize',8)
    plot(indC,outcomes(cues==2),'c.','MarkerSize',8); % plot the individual data points 
    plot([stateBlockDuration*0.5 stateBlockDuration*0.5],[-0.2 1.2],'--k',[stateBlockDuration stateBlockDuration],[-0.2 1.2],'--k',...
    [stateBlockDuration*1.5 stateBlockDuration*1.5],[-0.2 1.2],'--k');
    set(gca,'YLim',[-0.2 1.2]); set(gca,'XLim',[0 totalTrials]); legend('A-probs','C-probs','estimated q-AC','true H-AC','estimated H-AC','A-data','C-data');
    xlabel('trial number','FontSize',16); ylabel('model (joined AC) estimate of q','FontSize',16); set(gcf,'color','w'); set(gca,'FontSize',16);

    % plot estimates based on marginal means of distribution - AB
    figure; hold on;   
    title('model estimation - AB')
    plot(indA,probs(cues==0),'r--',indB,probs(cues==1),'g--'); plot(indAB,EqAB,'k'); plot(indAB,vol(cues==0 | cues==1),'k--'); plot(indAB,EHAB,'m');  
    plot(indA,outcomes(cues==0),'r.','MarkerSize',8)
    plot(indB,outcomes(cues==1),'g.','MarkerSize',8); % plot the individual data points 
    set(gca,'YLim',[-0.2 1.2]); set(gca,'XLim',[0 totalTrials]); legend('A-probs','B-probs','estimated q-AB','true H - AB','estimated H-AB','A-data','B-data');
    xlabel('trial number','FontSize',16); ylabel('model (joined AB) estimate of q','FontSize',16); set(gcf,'color','w'); set(gca,'FontSize',16);    
    
    % plot the prior over q on each trial (marginalizing over H) - only AC
    figure; hold on;
    prior_p_qAC=squeeze(sum(prior_p_qHAC,2)); % marginalize over H
    priorDispAC=nan(size(prior_p_qAC,1),totalTrials);
    priorDispAC(:,indAC)=prior_p_qAC;
    imagesc(1:length(outcomes),q_candidates,priorDispAC); %plot probability distribution over q on each trial
    colorbar
    plot(indA,probs(cues==0),'r',indC,probs(cues==2),'c'); % plot actual value of pL on each trial as a white line
    plot(indAC,EqAC,'w--'); 
    plot(indA,outcomes(cues==0),'r.','MarkerSize',14)
    plot(indC,outcomes(cues==2),'c.','MarkerSize',14); % plot the individual data points     set(gcf,'Color','w');
    set(gca,'YLim',[0 1],'XLim',[0 totalTrials])
    xlabel('trial number','FontSize',16)
end

% plot model avidence for joint and separate models - AC

figure
suptitle('Model evidence - AC')

subplot(4,1,1)
hold on
plot(indAC,scoreAC,'.') 
plot(indA,scoreA,'r.')
plot(indC,scoreC,'g.')
plot(indA,outcomes(cues==0)/5+1.1,'r.','MarkerSize',8)
plot(indC,outcomes(cues==2)/5+1.1,'c.','MarkerSize',8)
plot([stateBlockDuration*0.5 stateBlockDuration*0.5],[0 1.2],'--k')
plot([stateBlockDuration stateBlockDuration],[0 1.2],'--k',[stateBlockDuration*1.5 stateBlockDuration*1.5],[0 1.2],'--k')
set(gca,'YLim',[0 1.3],'XLim',[0 totalTrials])
legend('AC','A','C')
ylabel('Model Evidence')

subplot(4,1,2)
hold on
% scoreSepAC(1:totalTrials)=0;
% scoreSepAC(indA)=scoreA;
% scoreSepAC(indC)=scoreC;
% scoreSepAC(scoreSepAC==0)=[];
% plot(indAC,scoreAC-scoreSepAC,'-')
% for plotting A and B differences on separate lines 
plot(indA,log(scoreAC(find(ismember(indAC,indA))))-log(scoreA),'r-')
plot(indC,log(scoreAC(find(ismember(indAC,indC))))-log(scoreC),'c-')

plot(indA,probs(cues==0)+1.5,'r-.')
plot(indC,probs(cues==2)+1.5,'c-.')
plot([stateBlockDuration*0.5 stateBlockDuration*0.5],[-0.7 1],'--k')
plot([stateBlockDuration stateBlockDuration],[-0.7 1],'--k',[stateBlockDuration*1.5 stateBlockDuration*1.5],[-0.7 1],'--k')
set(gca,'YLim',[-1.5 2.5],'XLim',[0 totalTrials])
plot([0 totalTrials],[0 0])
title('diff of log(Evidence)')

subplot(4,1,3)
hold on
plot(indA,cumsum(log(scoreAC(find(ismember(indAC,indA)))))-cumsum(log(scoreA)),'r')
plot(indC,cumsum(log(scoreAC(find(ismember(indAC,indC))))) - cumsum(log(scoreC)),'c')
plot([0 totalTrials],[0 0],'k-')
legend('AC - A (only indA)','AC - C (only indC')
title({'diff ofcumsum(log(evidence))'})
 
subplot(4,1,4)
w = 10; % window size for calculating running average of the difference between log scores 
hold on
plot(indA,slidefun(@mean,w,log(scoreAC(find(ismember(indAC,indA)))) - log(scoreA)),'r')
plot(indC,slidefun(@mean,w,log(scoreAC(find(ismember(indAC,indC)))) - log(scoreC)),'c')
plot([0 totalTrials],[0 0],'k-')
title(sprintf('running avg - diff of log(evidence), window size = %g',w))



% plot model avidence for joint and separate models - AB

figure
suptitle('Model evidence - AB')
subplot(4,1,1)
hold on
plot(indAB,scoreAB,'.') 
plot(indA,scoreA,'r.')
plot(indB,scoreB,'g.')
plot(indA,outcomes(cues==0)/5+1,'r.','MarkerSize',8)
plot(indB,outcomes(cues==1)/5+1,'g.','MarkerSize',8)
plot([stateBlockDuration*0.5 stateBlockDuration*0.5],[0 1.2],'--k')
plot([stateBlockDuration stateBlockDuration],[0 1.2],'--k',[stateBlockDuration*1.5 stateBlockDuration*1.5],[0 1.2],'--k')
set(gca,'YLim',[0 1.2],'XLim',[0 totalTrials])
legend('AB','A','B')
ylabel('Model Evidence')

subplot(4,1,2)
hold on
% %for ploting both A and B differences on the same line:
% scoreSepAB(1:totalTrials)=0;
% scoreSepAB(indA)=scoreA;
% scoreSepAB(indB)=scoreB;
% scoreSepAB(scoreSepAB==0)=[];
% plot(indAB,scoreAB-scoreSepAB,'-')

% for plotting A and B differences on separate lines 
plot(indA,log(scoreAB(find(ismember(indAB,indA))))-log(scoreA),'r-')
plot(indB,log(scoreAB(find(ismember(indAB,indB))))-log(scoreB),'g-')

plot(indA,probs(cues==0)+1.5,'r-.')
plot(indB,probs(cues==1)+1.5,'g-.')
plot([stateBlockDuration*0.5 stateBlockDuration*0.5],[-0.7 1],'--k')
plot([stateBlockDuration stateBlockDuration],[-0.7 1],'--k',[stateBlockDuration*1.5 stateBlockDuration*1.5],[-0.7 1],'--k')
set(gca,'YLim',[-1.5 2.5],'XLim',[0 totalTrials])
plot([0 totalTrials],[0 0])
title('diff of log(Evidence)')

subplot(4,1,3)
hold on
plot(indA,cumsum(log(scoreAB(find(ismember(indAB,indA)))))-cumsum(log(scoreA)),'r')
plot(indB,cumsum(log(scoreAB(find(ismember(indAB,indB))))) - cumsum(log(scoreB)),'g')
plot([0 totalTrials],[0 0],'k-')
legend('AB - A (only indA)','AB - B (only indB')
title({'diff of cumsum(log(evidence))'})

subplot(4,1,4)
hold on
plot(indA,slidefun(@mean,w,log(scoreAB(find(ismember(indAB,indA)))) - log(scoreA)),'r')
plot(indB,slidefun(@mean,w,log(scoreAB(find(ismember(indAB,indB)))) - log(scoreB)),'g')
plot([0 totalTrials],[0 0],'k-')
title(sprintf('running avg - diff of log(evidence), window size = %g',w))



% figure
% plot(indA,scoreA,'o',[stateBlockDuration*0.5 stateBlockDuration*0.5],[0 1],'--k',...
%     [stateBlockDuration stateBlockDuration],[0 1],'--k',[stateBlockDuration*1.5 stateBlockDuration*1.5],[0 1],'--k',...
%     indA,outcomes(cues==0),'r.','MarkerSize',20)
