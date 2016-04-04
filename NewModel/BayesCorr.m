% Hierarichical Bayesian Learner
%
% 09-08-2013
% MPI.
%
% Adapted by Erie Boorman.
% Also see the original .cpp code by Tim Behrens.
clearvars -except schedules schedulesParams scores
close all

%% 

sesToLoad = 201; % from 'blocks.mat'

cues = squeeze(schedules.cues(sesToLoad,:));
outcomes = squeeze(schedules.outcomes(sesToLoad,:));
probs = squeeze(schedules.probs(sesToLoad,:));
totalTrials=schedulesParams.totalTrials;
stateBlockDuration = schedulesParams.stateBlockDuration;

tic;

%%

% q-axis (reward rate)
qVec = .01:.02:.99;
qSize = length(qVec);

% hLog axis (H = 1/v)
HLog = log(2):0.3:log(10000);
HSize = length(HLog);

% Mlog axis
MVec = .01:.05:.99;
MSize = length(MVec);

% klog axis
kLog = log(5e-4):0.3:log(20);
kSize = length(kLog);

% load tmp into pIk for each k
% pIk: 3D distribution that we're always updating p(p_i, I_i, k | every
% observation to date).
qAqBhAhBmk = ones(qSize, qSize, HSize, HSize, MSize, kSize) ./ (qSize*qSize*HSize*HSize*MSize*kSize);

% precompute p(Mi+1|Mi,k) for every M_{i+1}, M_i, k
Mp1gMk = zeros(MSize, MSize, kSize); % 3d array for M transition.
tmpM = zeros(MSize, MSize);  % to keep N(M_{i+1},k2) % Alon: it's actually N(M_{i},k2)
for k = 1 : kSize
    for M = 1:MSize
        for Mp1 = 1:MSize
            % N(M_{i+1},k2) % Alon: actually N(M_{i},k2)
            var = exp(kLog(k)*2); % k is stdev
            tmpM(Mp1, M) = (exp(-power((MVec(M) - MVec(Mp1)),2))/(2*var)) / (sqrt(2*pi*var));
        end
        % normalise so p(M_i+1|M_i,k) sums to 1;
        tmpM(:,M) = tmpM(:,M)./sum(tmpM(:,M));
    end
    Mp1gMk(:,:,k) = tmpM; % place tmpM in it.
end
 
%%
% precompute p(qA_{i+1}|qA_i,qB_i,hA, M_{i+1}, T_i) 

qAp1gqAqBhAT = zeros(qSize, qSize, qSize, HSize, MSize, 2);
for  HA = 1:HSize
    for Mp1 = 1:MSize              
        for qA = 1:qSize
            for qB = 1:qSize
                for T = 1:2
                    for qAp1 = 1:qSize
                        if T==1 
                        qAp1gqAqBhAT(qAp1,qA,qB,HA,Mp1,T) = exp(HLog(HA))*(1/qSize) + (1-exp(HLog(HA)))*(qVec(qAp1)==qVec(qA));
                        else
                        qAp1gqAqBhAT(qAp1,qA,qB,HA,Mp1,T) = exp(HLog(HA))*(1/qSize) + (1-exp(HLog(HA)))*((1-MVec(Mp1))*(qVec(qAp1)==qVec(qA)) + MVec(Mp1)*(qVec(qAp1)==qVec(qB)));
                        end
                    end
                end
            end
        end
    end
end

keyboard
                            
                           
                    





% tmpp = zeros(pSize, pSize);
% for Ip1 = 1:Isize
%     for p = 1:pSize
%         for pp1 = 1:pSize
%             % Beta(x; a,b)
%             a = 1 + (exp(Ilog(Ip1)))*pvec(p);
%             b = 1 + (exp(Ilog(Ip1)))*(1-pvec(p));
%             x = pvec(pp1);
%             if ~(x==0)||~(x==1)
%                 logkerna = (a-1)*log(x);
%                 logkernb = (b-1)*log(1-x);
%                 betaln_ab = gammaln(a) + gammaln(b) - gammaln(a+b);
%                 tmpp(pp1, p) = exp(logkerna + logkernb - betaln_ab);
%             else
%                 tmpp(pp1, p) = 0;
%             end
%             
%         end
%         tmpp(:,p) = tmpp(:,p)./sum(tmpp(:,p));
%     end
%     pp1gpIp1(:,:,Ip1) = tmpp; % place tmpp in it.
% end

%% 

%
% BAYESIAN UPDATE
%

for trial = 1:length(reward)
    if (reward(trial)==1)
        for k = 1:kSize
            for p=1:pSize
                pIk(p,:,k) = pIk(p,:,k)*pvec(p);
            end
        end
    else
        for k = 1:kSize
            for p=1:pSize
                pIk(p,:,k) = pIk(p,:,k)*(1-pvec(p));
            end
        end
        
    end
    
    % now do normalization
    pIk = pIk ./ sum(sum(sum(pIk)));
    
    
    % GET MARGINALS
    %
    
    % p(reward)
    pI = sum(pIk,3);
    pDist(trial,:) = sum(pI,2); % sum over each row.
    pEst(trial,:)  = sum(pDist(trial,:).*pvec);
    
    % p(volatility)
    pI = sum(pIk,3);
    vDist(trial,:) = sum(pI); % column-wise sum.
    vEst(trial,:)  = sum(vDist(trial,:).*Ilog);    
    
    
    %
    % INFORMATION LEAK (increase the variance in the joint distribution)
    %
    
    % I) multiply pIk by Ip1gIk, and integrate out I. This will give pIp1k.
    for k = 1:kSize
        pIp1k = zeros(pSize,Isize);
        for Ip1 = 1:Isize
            for p = 1:pSize
                pIp1k(p,Ip1) = sum(Ip1gIk(Ip1,:,k).*pIk(p,:,k));
            end
        end
        % II) multiply pIp1k by pp1gpIp1, and integrate out p. This will give pp1Ip1k.
        pp1Ip1k = zeros(pSize,Isize);
        for Ip1 = 1:Isize
            for pp1 = 1:pSize
                pp1Ip1k(pp1,Ip1) = sum(pIp1k(:,Ip1).*pp1gpIp1(pp1,:,Ip1)');
            end
        end
        % III) Place pp1Ip1k into pIk (belief that is carried to the next
        % trial.
        pIk(:,:,k) = pp1Ip1k;
    end
       
 end



t=toc;


% make some plots
plot(1./vEst)   % volatility
probSched = [ones(1,60)'*.75; ones(1,20)'*.20; ones(1,20)'*.80; ones(1,20)'*.20]; % experimental setup
%probSchedAdvice = [ones(1,30)'*.75; ones(1,10)'*.20; ones(1,10)'*.80; ones(1,10)'*.20; ones(1,10)'*.80; ones(1,50)'*.15];
figure, plot(pEst) % reward rate
hold on, plot(probSched,'r')
ylim([0 1])




