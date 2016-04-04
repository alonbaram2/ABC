  
close all

%% 

load('toLoadTest.mat')

tic;

%%

% q-axis (reward rate)
qVec = .01:.02:.99;
qSize = length(qVec);

% hLog axis (h = 1/v)
hLog = log(2):0.3:log(10000);
hSize = length(hLog);

% mlog axis
mVec = .01:.05:.99;
mSize = length(mVec);

% klog axis
kLog = log(5e-4):0.3:log(20);
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
vAEst = zeros(1,length(outcomes));
vBEst = zeros(1,length(outcomes));
mEst = zeros(1,length(outcomes));
kEst = zeros(1,length(outcomes));
hAEstExp = zeros(1,length(outcomes));
hBEstExp = zeros(1,length(outcomes));
vAEstExp = zeros(1,length(outcomes));
vBEstExp = zeros(1,length(outcomes));
kEstExp = zeros(1,length(outcomes));


% precompute p(mi+1|mi,k) for every m_{i+1}, m_i, k
mp1gmk = zeros(mSize, mSize, kSize); % 3d array for m transition.
tmpM = zeros(mSize, mSize);  % to keep N(m_{i},k2)
for k = 1 : kSize
    for m = 1:mSize
        for mp1 = 1:mSize
            %  N(m_{i},k2)
            var = exp(kLog(k)*2); % k is stdev
            tmpM(mp1, m) = (exp(-power((mVec(m) - mVec(mp1)),2))/(2*var)) / (sqrt(2*pi*var));
        end
        % normalise so p(m_i+1|m_i,k) sums to 1;
        tmpM(:,m) = tmpM(:,m)./sum(tmpM(:,m));
    end
    mp1gmk(:,:,k) = tmpM; % place tmpm in it.
end

% precompute p(qA_{i+1}|qA_i,qB_i,hA, m_{i+1}, T_i) 

qAp1gqAqBhAmp1T = zeros(qSize, qSize, qSize, hSize, mSize, 2);
qBp1gqBqAhBmp1T = zeros(qSize, qSize, qSize, hSize, mSize, 2);
for  hA = 1:hSize
    for mp1 = 1:mSize
        jump = (ones(qSize)./qSize).*(exp(-hLog(hA)));
        for qB=1:qSize            
            T = 1; % Current cue is A
            noJump = eye(qSize).*(1-exp(-hLog(hA)));
            qAp1gqAqBhAmp1T(:,:,qB,hA,mp1,T) = jump + noJump;
            
            T = 2; % Current cue is B
            noJumpCorrWeighted = eye(qSize).*(1-exp(-hLog(hA))).*(1-mVec(mp1));
            qAp1gqAqBhAmp1T(:,:,qB,hA,mp1,T) = jump + noJumpCorrWeighted;
        end
        T = 2;
        for qA=1:qSize
            corrInfluence = eye(qSize).*(1-exp(-hLog(hA))).*(mVec(mp1));            
            qAp1gqAqBhAmp1T(:,qA,:,hA,mp1,T) = qAp1gqAqBhAmp1T(:,qA,:,hA,mp1,T) + permute(corrInfluence,[1,3,2]);
        end
    end
end

 %p(qB_{i+1}|qB_i,qA_i,hB, m_{i+1}, T_i) is symmetrical for swithcing A and B 
qBp1gqBqAhBmp1T(:,:,:,:,:,2) = qAp1gqAqBhAmp1T(:,:,:,:,:,1);
qBp1gqBqAhBmp1T(:,:,:,:,:,1) = qAp1gqAqBhAmp1T(:,:,:,:,:,2);
   
%%

display('start Bayesian updating')
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
            qAqBhAhBmk(:,qB,:,:,:,:) = qAqBhAhBmk(:,qB,:,:,:,:)*(1-qVec(qA));
        end              
    end
        
    % now do normalization
    qAqBhAhBmk = qAqBhAhBmk ./ sum(sum(sum(sum(sum(sum(qAqBhAhBmk))))));
    
        % GET mARGINALS
    %
    
    % qA
    qADist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),5),4),3),2);
    qAEst(t)  = sum(qADist(t,:).*qVec);
    
    
    % qB
    qBDist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),5),4),3),1);
    qBEst(t)  = sum(qBDist(t,:).*qVec);
        
    % hA
    hADist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),5),4),2),1);
    hAEst(t)  = sum(hADist(t,:).*hLog);
    vAEst(t)  = sum(hADist(t,:).*(-hLog));    
    hAEstExp(t)  = sum(hADist(t,:).*exp(hLog));
    vAEstExp(t)  = sum(hADist(t,:).*exp(-hLog));
    
    % hB
    hBDist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),5),3),2),1);
    hBEst(t)  = sum(hBDist(t,:).*hLog);
    vBEst(t)  = sum(hBDist(t,:).*(-hLog));    
    hBEstExp(t)  = sum(hBDist(t,:).*exp(hLog));
    vBEstExp(t)  = sum(hBDist(t,:).*exp(-hLog));    
    
    
    % m
    mDist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,6),4),3),2),1);
    mEst(t)  = sum(mDist(t,:).*mVec);  
    
    %k
    kDist(t,:) = sum(sum(sum(sum(sum(qAqBhAhBmk,5),4),3),2),1);
    kEst(t)  = sum(kDist(t,:).*kLog);
    kEstExp(t) = sum(kDist(t,:).*exp(kLog));
        
    %
    % INFORMATION LEAK (increase the variance in the joint distribution)
    %
    
    display('Information Leak, trial = %d',t)
    toc
    
    % I) multiply qAqBhAhBmk by mp1gmk, and integrate out m. This will give qAqBhAhBmp1k.
    for k = 1:kSize
        qAqBhAhBmp1k = zeros(qSize,qSize,hSize,hSize,mSize);
        for mp1 = 1:mSize
            for hA = 1:hSize
                for hB = 1:hSize
                    for qA = 1:qSize
                        for qB = 1:qSize                            
                             qAqBhAhBmp1k(qA,qB,hA,hB,mp1) = sum(squeeze(mp1gmk(mp1,:,k)).*squeeze(qAqBhAhBmk(qA,qB,hA,hB,:,k))');
                        end
                    end
                end
            end
        end
        
        sprintf('end stage 1, k = %d',k)
        toc
        
        % II) multiply qAqBhAhBmp1k (pIp1k) by qAp1gqAqBhAmp1T (pp1gpIp1) and integrate out qA, 
        % This will give qAp1qBhAhBmp1k (pp1Ip1k).
        qAp1qBhAhBmp1k = zeros(qSize,qSize,hSize,hSize,mSize);
        for mp1 = 1:mSize
            for hA = 1:hSize
                for hB = 1:hSize
                    for qB = 1:qSize
                        for qAp1 = 1:qSize                          
                            qAp1qBhAhBmp1k(qAp1,qB,hA,hB,mp1) = sum(squeeze(qAqBhAhBmp1k(:,qB,hA,hB,mp1)).*squeeze(qAp1gqAqBhAmp1T(qAp1,:,qB,hA,mp1,T))');
                        end
                    end
                end
            end
        end
    sprintf('end stage 2, k = %d',k)
    toc       
        
        % III) multiply qAp1qBp1hAhBmp1k by qBp1gqBqAhBmp1T and integrate out qB, 
        % This will give qAp1qBp1hAhBmp1k (pp1Ip1k).  
        qAp1qBp1hAhBmp1k = zeros(qSize,qSize,hSize,hSize,mSize);        
        for mp1 = 1:mSize
            for hA = 1:hSize
                for hB = 1:hSize
                    for qAp1 = 1:qSize
                        for qBp1 = 1:qSize         
                            qAp1qBp1hAhBmp1k(qAp1,qBp1,hA,hB,mp1) = sum(squeeze(qAp1qBhAhBmp1k(qAp1,:,hA,hB,mp1)).*squeeze(qBp1gqBqAhBmp1T(qBp1,:,qAp1,hB,mp1,T))); % the last qApq should actually be qA but it's the same index so it doesn't matter.
                        end
                    end
                end
            end
        end        
        
    sprintf('end stage 3, k = %d',k)
    toc       
        % IV) Place qAp1qBp1hAhBmp1k into qAqBhAhBmp1k (belief that is carried to the next
        % trial).
        qAqBhAhBmp1k(:,:,:,:,:,k) = qAp1qBp1hAhBmp1k;
             
    end
    sprintf('end Bayesian update, trial = %d',t)
    toc       
 end

timeOverall = toc

blocks201 = struct('qADist',qADist,'qBDist',qBDist,'hADist',hADist,'hBDist',hBDist,'mDist',mDist,'kDist',kDist,...
                    'qAEst',qAEst,'qBEst',qBEst,'hAEst',hAEst,'hBEst',hBEst,'vAEst',vAEst,'vBEst',vBEst,'mEst',mEst,'kEst',kEst,...
                    'hAEstExp',hAEstExp,'hBEstExp',hBEstExp,'vAEstExp',vAEstExp,'vBEstExp',vBEstExp,'kEstExp',kEstExp,'elapsedTime',elapsedTime);

save('blocks201.mat','blocks201')               
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
