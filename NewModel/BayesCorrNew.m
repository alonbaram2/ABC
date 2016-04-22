
% close all

%%
% 
% AntiCorrStabTotal.mat CorrStabTotal.mat AntiCorrVol10.mat CorrVol10Long.mat h0m0AOnlyVol.mat
sesName = 'AC2C8-2ContM.mat';
% load(strcat('toLoad',sesName));
unmodeledCue = 2; % a numbers from {0,1,2}.
del = find(cues==unmodeledCue);
cues(del)=[];
outcomes(del)=[];
probs(del)=[];

tic;

%%

mLogModel = 0; % 0 is -1 < m < 1 in normal space, 1 is 0 < m < 1 in log space.

% q-axis (reward rate)
qVec = .01:.1:.99;
qSize = length(qVec);

% hLog axis
hLog = log(1/1000):0.8:log(1/2);
hSize = length(hLog);

% m-axis  %to have the sime mSize=14 in both options choose resolutions of
% .65 and .15 respectively.
if mLogModel==1
    mVec = log(1/10000):0.65:log(1/2);
else
    mVec = -.99:.15:.99;
end 
mSize = length(mVec);

% klog axis
kLog = log(3e-2):0.9:log(0.5);
kSize = length(kLog);

% qAqBhAhBmk: The joint 6D distribution that we're always updating p(qA, qB, hA, hB, m, k | every
% observation to date).
qAqBhAhBmk = ones(qSize, qSize, hSize, hSize, mSize, kSize) ./ (qSize*qSize*hSize*hSize*mSize*kSize);

qADist = zeros(length(outcomes),qSize);
qBDist = zeros(length(outcomes),qSize);
hADist = zeros(length(outcomes),hSize);
hBDist = zeros(length(outcomes),hSize);
mDist = zeros(length(outcomes),mSize);
kDist = zeros(length(outcomes),kSize);
qAEst = zeros(1,length(outcomes));
qBEst = zeros(1,length(outcomes));
hAEst = zeros(1,length(outcomes));
hBEst = zeros(1,length(outcomes));
mEst = zeros(1,length(outcomes));
kEst = zeros(1,length(outcomes));
hAEstExp = zeros(1,length(outcomes));
hBEstExp = zeros(1,length(outcomes));
mEstExp = zeros(1,length(outcomes));
kEstExp = zeros(1,length(outcomes));

% precompute p(mi+1|mi,k) for every m_{i+1}, m_i, k
mp1gmk = zeros(mSize, mSize, kSize); % 3d array for m transition.
tmpM = zeros(mSize, mSize);  % to keep N(m_{i},k2)
for k = 1 : kSize
    for m = 1:mSize
        for mp1 = 1:mSize
            %  N(m_{i},k2)
            var = exp(kLog(k)*2); % k is stdev
            tmpM(mp1, m) = (exp(-power((mVec(m) - mVec(mp1)),2)/(2*var))) / (sqrt(2*pi*var));
        end
        % normalise so p(m_i+1|m_i,k) sums to 1;
        tmpM(:,m) = tmpM(:,m)./sum(tmpM(:,m));
    end
    
    mp1gmk(:,:,k) = tmpM; % place tmpm in it.
end

% precompute p(qA_{i+1}|qA_i,hA, T_{i+1) = 1)  

qAp1gqAhAT1 = zeros(qSize, qSize, hSize);

for hA = 1:hSize
    jump = (ones(qSize)./qSize).*(exp(hLog(hA)));
    noJump = eye(qSize).*(1-exp(hLog(hA)));
    qAp1gqAhAT1(:,:,hA) = jump + noJump;
end
% precompute p(qB_{i+1}|qB_i,hB, T_{i+1) = 2)  
qBp1gqBhBT2 = qAp1gqAhAT1;

% precompute p(qA_{i+1}|qA_i,qB_i,m_{i+1}, T_{i+1) = 2)
qAp1gqAqBmp1T2 = zeros(qSize, qSize, qSize, mSize);
if mLogModel == 1
    for mp1 = 1:mSize
        for qB = 1:qSize
            qAp1gqAqBmp1T2(:,:,qB,mp1) = eye(qSize).*(1-exp(mVec(mp1)));
        end
        for qA = 1:qSize
            qAp1gqAqBmp1T2(:,qA,:,mp1) = qAp1gqAqBmp1T2(:,qA,:,mp1) + permute(eye(qSize).*(exp(mVec(mp1))),[1,3,2]);
        end
    end
elseif mLogModel == 0 
    for mp1 = 1:mSize
        for qA = 1:qSize
            for qB = 1:qSize
                tmp = (1-abs(mVec(mp1)))*(qVec(qA)-0.5) + (mVec(mp1))*(qVec(qB) - 0.5) + 0.5;
                [~, qAp1] = min(power((qVec-tmp),2));
                qAp1gqAqBmp1T2(qAp1,qA,qB,mp1) = 1;
            end
        end
    end
end

% precompute p(qB_{i+1}|qB_i,qA_i,m_{i+1}, T_{i+1) = 1)
qBp1gqBqAmp1T1 = qAp1gqAqBmp1T2;

%%

display(sprintf('start Bayesian updating'))
toc

%
% BAYESIAN UPDATE
%

for t = 1:length(outcomes)
    if (cues(t)==0 && outcomes(t) ==1)
        T = 1;
        for qA=1:qSize
            qAqBhAhBmk(qA,:,:,:,:,:) = qAqBhAhBmk(qA,:,:,:,:,:)*qVec(qA);
        end
    elseif (cues(t)==0 && outcomes(t) ==0)
        T = 1;
        for qA=1:qSize
            qAqBhAhBmk(qA,:,:,:,:,:) = qAqBhAhBmk(qA,:,:,:,:,:)*(1-qVec(qA));
        end
    elseif (cues(t)==1 && outcomes(t)==1)
        T = 2;
        for qB=1:qSize
            qAqBhAhBmk(:,qB,:,:,:,:) = qAqBhAhBmk(:,qB,:,:,:,:)*qVec(qB);
        end
    elseif (cues(t)==1 && outcomes(t)==0)
        T = 2;
        for qB=1:qSize
            qAqBhAhBmk(:,qB,:,:,:,:) = qAqBhAhBmk(:,qB,:,:,:,:)*(1-qVec(qB));
        end
    end
    
    % now do normalization
    qAqBhAhBmk = qAqBhAhBmk ./ sum(sum(sum(sum(sum(sum(qAqBhAhBmk))))));
    
    
    % GET MARGINALS
    
    % qA
    qADist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),5),4),3),2);
    qAEst(t)  = sum(qADist(t,:).*qVec);
    % qB
    qBDist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),5),4),3),1);
    qBEst(t)  = sum(qBDist(t,:).*qVec);
    % hA
    hADist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),5),4),2),1);
    hAEst(t)  = sum(hADist(t,:).*hLog);
    hAEstExp(t)  = sum(hADist(t,:).*exp(hLog));
    % hB
    hBDist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),5),3),2),1);
    hBEst(t)  = sum(hBDist(t,:).*hLog);
    hBEstExp(t)  = sum(hBDist(t,:).*exp(hLog));
    % m
    mDist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),4),3),2),1);
    mEstExp(t) = sum(mDist(t,:).*exp(mVec));
    mEst(t)  = sum(mDist(t,:).*mVec);
    % k
    kDist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,5),4),3),2),1);
    kEst(t)  = sum(kDist(t,:).*kLog);
    kEstExp(t) = sum(kDist(t,:).*exp(kLog));
    
    
    % INFORMATION LEAK (increase the variance in the joint distribution)
    %
    
    %     sprintf('Information Leak, trial = %d',t)
    %     toc
    %
    % I) multiply qAqBhAhBmk by mp1gmk, and integrate out m. This will give qAqBhAhBmp1k.
    
    for k = 1:kSize
        qAqBhAhBmp1k = zeros(qSize,qSize,hSize,hSize,mSize);
        for mp1 = 1:mSize
            for hA = 1:hSize
                for hB = 1:hSize
                    for qA = 1:qSize
                        for qB = 1:qSize
                            qAqBhAhBmp1k(qA,qB,hA,hB,mp1) = sum(reshape(mp1gmk(mp1,:,k),mSize,1).*reshape(qAqBhAhBmk(qA,qB,hA,hB,:,k),mSize,1));
                        end
                    end
                end
            end
        end
        
        %         sprintf('end stage 1, k = %d, trial = %d',k,t)
        %         toc
        %
        % II) multiply qAqBhAhBmp1k (pIp1k) by qAp1gqAqBhAmp1T (pp1gpIp1) and integrate out qA,
        % This will give qAp1qBhAhBmp1k (pp1Ip1k).
        
        if T==1
            tmp = zeros(qSize,qSize,qSize,qSize,hSize,hSize,mSize);           
            qAp1qBp1hAhBmp1k = zeros(qSize,qSize,hSize,hSize,mSize);
            for qAp1 = 1:qSize
                for qBp1 = 1:qSize
                    for hA = 1:hSize
                        for mp1 = 1:mSize
                            for qA = 1:qSize
                                for qB = 1:qSize
                                    tmp(qAp1,qBp1,qA,qB,hA,:,mp1) = qAqBhAhBmp1k(qA,qB,hA,:,mp1)*...
                                        qAp1gqAhAT1(qAp1,qA,hA)*qBp1gqBqAmp1T1(qBp1,qB,qA,mp1);
                                end
                            end
                            qAp1qBp1hAhBmp1k(qAp1,qBp1,hA,:,mp1) = sum(sum(tmp(qAp1,qBp1,:,:,hA,:,mp1),4),3);
                        end
                    end
                end
            end
            qAqBhAhBmk(:,:,:,:,:,k) = qAp1qBp1hAhBmp1k;
            
        elseif T==2
            
            tmp = zeros(qSize,qSize,qSize,qSize,hSize,hSize,mSize);            
            qAp1qBp1hAhBmp1k = zeros(qSize,qSize,hSize,hSize,mSize);        
            for qAp1 = 1:qSize
                for qBp1 = 1:qSize
                    for hB = 1:hSize
                        for mp1 = 1:mSize
                            for qA = 1:qSize
                                for qB = 1:qSize
                                    tmp(qAp1,qBp1,qA,qB,:,hB,mp1) = qAqBhAhBmp1k(qA,qB,:,hB,mp1)*qBp1gqBhBT2(qBp1,qB,hB)*qAp1gqAqBmp1T2(qAp1,qA,qB,mp1);
                                end
                            end
                            qAp1qBp1hAhBmp1k(qAp1,qBp1,:,hB,mp1) = sum(sum(tmp(qAp1,qBp1,:,:,:,hB,mp1),4),3);
                        end
                    end
                end
            end
            qAqBhAhBmk(:,:,:,:,:,k) = qAp1qBp1hAhBmp1k;

            
            
%             tmp = zeros(qSize,qSize,qSize,qSize,hSize,mSize);            
%             qAp1qBp1hBmp1k = zeros(qSize,qSize,hSize,mSize);        
%             qAqBhBmp1k = reshape(sum(qAqBhAhBmp1k,3),[qSize,qSize,hSize,mSize]);
%             for qAp1 = 1:qSize
%                 for qBp1 = 1:qSize
%                     for hB = 1:hSize
%                         for mp1 = 1:mSize
%                             for qA = 1:qSize
%                                 for qB = 1:qSize
%                                     tmp(qAp1,qBp1,qA,qB,hB,mp1) = qAqBhBmp1k(qA,qB,hB,mp1)*qBp1gqBhBT2(qBp1,qB,hB)*qAp1gqAqBmp1T2(qAp1,qA,qB,mp1);
%                                 end
%                             end
%                             qAp1qBp1hBmp1k(qAp1,qBp1,hB,mp1) = sum(sum(tmp(qAp1,qBp1,:,:,hB,mp1),4),3);
%                         end
%                     end
%                 end
%             end
%             for hA = 1:hSize
%                 qAqBhAhBmk(:,:,hA,:,:,k) = permute(qAp1qBp1hBmp1k,[1 2 5 3 4 6]);
%             end
        end
        
        
        
        %         sprintf('end stage 2, k = %d, trial = %d',k,t)
        %         toc
        %
        %         % III) multiply qAp1qBp1hAhBmp1k by qBp1gqBqAhBmp1T and integrate out qB,
        %         % This will give qAp1qBp1hAhBmp1k (pp1Ip1k).
        %         qAp1qBp1hAhBmp1k = zeros(qSize,qSize,hSize,hSize,mSize);
        %         for mp1 = 1:mSize
        %             for hA = 1:hSize
        %                 for hB = 1:hSize
        %                     for qAp1 = 1:qSize
        %                         for qBp1 = 1:qSize
        %                             qAp1qBp1hAhBmp1k(qAp1,qBp1,hA,hB,mp1) = sum(squeeze(qAp1qBhAhBmp1k(qAp1,:,hA,hB,mp1)).*squeeze(qBp1gqBqAhBmp1T(qBp1,:,qAp1,hB,mp1,T))); % the last qApq should actually be qA but it's the same index so it doesn't matter.
        %                         end
        %                     end
        %                 end
        %             end
        %         end
        %
        %     sprintf('end stage 3, k = %d, trial = %d',k,t)
        %         toc
        % IV) Place qAp1qBp1hAhBmp1k into qAqBhAhBmqk (belief that is carried to the next
        % trial).
        
    end
    display(sprintf('end Bayesian update, trial = %d',t))
    toc
end

timeTotal = toc

modelStr = sprintf('h0m%d',mLogModel); % h0 is for no jump on oppite trial, m0 is for normal space prior for m with the option rof anti-correlation. 

save(strcat(modelStr,sesName),'qADist','qBDist','hADist','hBDist','mDist','kDist',...
    'qAEst','qBEst','hAEst','hBEst','vAEst','vBEst','mEst','kEst',...
    'hAEstExp','hBEstExp','vAEstExp','vBEstExp','mEstExp','kEstExp','timeTotal',...
    'qVec','hLog','mVec','kLog','cues','outcomes','probs','mLogModel','sesName');
% 
% 


% % make some plots
% X = 1:length(reward);
%
% figure
% suptitle(str)
% subplotNum = 6;
% subplot(subplotNum,1,1)
% plot (X(reward==0),IEst(reward==0),'.r',X(reward==1),IEst(reward==1),'.b') , title('IEst')
% hold on
% subplot(subplotNum,1,2)
% plot(X(reward==0),1./IEst(reward==0),'.r',X(reward==1),1./IEst(reward==1),'.b')  , title('1./IEst')  % volatility
% subplot(subplotNum,1,3)
% plot(X(reward==0),IEstExp(reward==0),'.r',X(reward==1),IEstExp(reward==1),'.b')  , title('IEstExp')
% subplot(subplotNum,1,4)
% plot(X(reward==0),1./IEst(reward==0),'.r',X(reward==1),1./IEst(reward==1),'.b') , title('1./IEstExp')
% subplot(subplotNum,1,5)
% plot(X(reward==0),vEst(reward==0),'.r',X(reward==1),vEst(reward==1),'.b') , title('vEst')
% subplot(subplotNum,1,6)
% plot(X(reward==0),vEstExp(reward==0),'.r',X(reward==1),vEstExp(reward==1),'.b') , title('vEstExp (the important one)')
%
% figure
% suptitle(str)
% hold on
% subplot(2,1,1)
% plot (kEst) , title('kEst')
% subplot(2,1,2)
% plot(kEstExp), title('kEstExp')
%
%
% probSched = [ones(1,120)'*.75; ones(1,40)'*.20; ones(1,40)'*.80; ones(1,40)'*.20; ones(1,40)'*.80]; % experimental setup
% %probSchedAdvice = [ones(1,30)'*.75; ones(1,10)'*.20; ones(1,10)'*.80; ones(1,10)'*.20; ones(1,10)'*.80; ones(1,50)'*.15];
% figure, plot(X(reward==0),pEst(reward==0),'.r',X(reward==1),pEst(reward==1),'.b') % reward rate
% suptitle(str)
% hold on, plot(probSched,'r')
% ylim([0 1])
%
%
%
%
