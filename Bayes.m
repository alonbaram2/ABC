% Alon Baram (based on code from Jill O'Rielly
% ------------------------------------------------------------------------

function [score,prior_p_qH,q_candidates,H_candidates] = Bayes(q,y,H,ind)

% q=load('data/schedule_volatile.txt');
% y=load('data/trials_volatile.txt');

%% now set up the state space
% possible values for p and j
q_candidates=(0.01:0.01:0.99)';
H_candidates=exp(log(0.01):(log(0.3)-log(0.01))/20:log(0.3))';
score = NaN(1,length(y));

% grids
[qq,HH]=ndgrid(q_candidates,H_candidates);

% transition function
transfunc=  (reshape(repmat(eye(length(q_candidates)),1,length(H_candidates)),length(q_candidates),length(q_candidates),length(H_candidates))... % p(pL(t)| no jump occurred) * p(no jump occurred)
            .* permute(reshape(repmat(1-H_candidates,length(q_candidates),length(q_candidates)),length(H_candidates),length(q_candidates),length(q_candidates)),[2 3 1]))...
            + ((ones(length(q_candidates),length(q_candidates),length(H_candidates))./length(q_candidates))... % + p(pL(t)| jump occurred) * p(jump_occurred)
            .* permute(reshape(repmat(H_candidates,length(q_candidates),length(q_candidates)),length(H_candidates),length(q_candidates),length(q_candidates)),[2 3 1]));

prior_p_qH=NaN(size(qq,1),size(HH,2),length(y)); %initiate prior and post matrices
post_p_qH=NaN(size(qq,1),size(HH,2),length(y));
% uniform prior to start with
prior_p_qH(:,:,1)=ones(size(qq))./length(qq(:));

%% fit the model
for i=1:length(y)    

    % 'leak' or apply transition function
    if i>1
        unpacked_post=permute(reshape(repmat(post_p_qH(:,:,i-1),1,length(q_candidates)),length(q_candidates),length(H_candidates),length(q_candidates)),[1 3 2]);
        unpacked_prior=unpacked_post.*transfunc;
        prior_p_qH(:,:,i)=squeeze(sum(unpacked_prior,1));
        prior_p_qH(:,:,i)=prior_p_qH(:,:,i)./sum(sum(prior_p_qH(:,:,i)));
    end
        
    % p(y|q,H)
    if y(i)==1
        % if y(i) is a hit, p(y(i)|q) is simply q
        py_given_qH=qq;
    else
        % if y(i) is a miss, p(y(i)|q) is 1-q
        py_given_qH=1-qq;
    end
    
    % Bayes' theorem: p(q,H|y_1:i)=p(y|q_i,H)p(q_1:i-1,H);
    post_p_qH(:,:,i)=py_given_qH.*prior_p_qH(:,:,i);
    
    %normalise posterior
    score(i) = sum(sum(post_p_qH(:,:,i)));
    post_p_qH(:,:,i)=post_p_qH(:,:,i)./sum(sum(post_p_qH(:,:,i)));
    % get joint maximum of prior for q and H
    [mx,ix]=max(myvect(prior_p_qH(:,:,i)));
    [ixq,ixH]=ind2sub(size(qq),ix);
    joint_mode(i,:)=[q_candidates(ixq), H_candidates(ixH)];
    
    % get marginal expected values for q and H
    Eq(i)=sum(myvect(prior_p_qH(:,:,i).*qq));
    EH(i)=sum(myvect(prior_p_qH(:,:,i).*HH));
    
end

    
%     %% plot true q and model's estimate
%     figure; hold on;   
%     
%     % plot estimates based on mode/peak of distribution
%     %plot(joint_mode(:,1),'r'); plot([0 500],[1/25 1/25],'c--'); plot(joint_mode(:,2),'g'); plot(q,'b'); 
%      
%     % plot estimates based on marginal means of distribution
%     plot(ind,q,'b'); plot(ind,Eq,'r'); plot(ind,H,'c--'); plot(ind,EH,'g');  
%     
%     plot(ind,y,'k.','MarkerSize',8); % plot the individual data points 
%     set(gca,'YLim',[-0.2 1.2]); legend('true q','estimated q','true H','estimated H','Data points');
%     xlabel('trial number','FontSize',16); ylabel('model estimate of q','FontSize',16); set(gcf,'color','w'); set(gca,'FontSize',16);
% 
%     %% Plot joint probability distribution over q and H
%     figure; hold on;
%     for i=1:60
%             subplot(6,10,i);
%             imagesc(H_candidates,q_candidates,prior_p_qH(:,:,i)); ylabel('value of q'); xlabel('value of H'); 
%             title(['trial ' int2str(i)]);
%     end
%     set(gcf,'color','w');
% 
%     %% plot the prior over q on each trial (marginalizing over H)
%     figure; hold on;
%     prior_p_q=squeeze(sum(prior_p_qH,2)); % marginalize over H
%     imagesc(1:length(y),q_candidates,prior_p_q); %plot probability distribution over q on each trial
%     plot(ind,q,'w'); % plot actual value of pL on each trial as a white line
%     plot(ind,Eq,'k'); plot(ind,Eq,'w--');% plot peak value for q on each trial as black and white dashed line
%     
%     
%     plot(y,'k.','MarkerSize',20);
%     set(gcf,'Color','w');



