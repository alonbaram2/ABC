% Hierarichical Bayesian Learner
%
% 09-08-2013
% MPI.
%
% Adapted by Erie Boorman.
% Also see the original .cpp code by Tim Behrens.
clearvars -except reward schedules schedulesParams scores

% clc, close all

%% 

% load Subject1.txt;
%load Subject2.txt;
% reward = Subject1(:,5); % 5th column => A correct?
%reward = Subject1(:,6);  % 6th column => Advice correct?

if (exist('reward')==0) 
    reward = [rand(1,120)<0.75 rand(1,40)<0.2 rand(1,40)<0.8 rand(1,40)<0.2 rand(1,40)<0.8]; 
end
tic;

%%
str = 'no jump';
jumpModel = 0; 
% p-axis (reward rate)
pvec = .01:.02:.99;
pSize = length(pvec);

% Ilog axis 
% Ilog = log(0.01):(log(0.3)-log(0.01))/20:log(0.3); 
% Ilog = log(1/10000):0.2:log(1/2); % vlog
Ilog = log(2):0.2:log(10000); % in Erie's code because I = 1/v ;
% Ilog = log(1/300):.2:log(0.3);
Isize = length(Ilog);

% klog axis
% klog = log(1/3000):.2:log(0.5);
klog = log(5/10000):.2:log(20);
kSize = length(klog);

tmp = ones(pSize, Isize)*(1/(pSize*Isize*kSize));
% load tmp into pIk for each k
% pIk: 3D distribution that we're always updating p(p_i, I_i, k | every
% observation to date).
pIk = zeros(pSize, Isize, kSize);
for k = 1 : kSize
    pIk(:,:,k) = tmp;
end


% precompute p(Ii+1|Ii,k) for every I_{i+1}, I_i, k
Ip1gIk = zeros(Isize, Isize, kSize); % 3d array for I transition.
tmpI = zeros(Isize, Isize);  % to keep N(I_{i+1},k2) % Alon: actually N(I_{i},k2)
var_all = zeros(1,kSize);
for k = 1 : kSize
    for I = 1:Isize
        for Ip1 = 1:Isize
            % N(I_{i+1},k2) % Alon: actually N(I_{i},k2)
            var = exp(klog(k)*2); % k is stdev
            var_all(k)=var;
            tmpI(Ip1, I) = (exp(-power((Ilog(I) - Ilog(Ip1)),2)/(2*var))) / (sqrt(2*pi*var));
        end
        % normalise so p(i_i+1|I_i,k) sums to 1;
        tmpI(:,I) = tmpI(:,I)./sum(tmpI(:,I));
    end
    Ip1gIk(:,:,k) = tmpI; % place tmpI in it.
end


% precompute p(pi+1|pi,Ii+1) for every p_{i+1}, p_i, I_{i+1}
if jumpModel
    pp1gpIp1 = zeros(pSize, pSize, Isize);
    for Ip1=1:Isize
        noJump = eye(pSize).*(1-exp(-Ilog(Ip1)));
        jump = (ones(pSize)./pSize).*(exp(-Ilog(Ip1)));
        pp1gpIp1(:,:,Ip1) = jump + noJump;
    end
    
        

    
%     
%     
%     pp1gpIp1 = (reshape(repmat(eye(pSize),1,Isize),pSize,pSize,Isize)... % no jump occurred
%         .* permute(reshape(repmat(1-exp(-Ilog),pSize,pSize),Isize,pSize,pSize),[2 3 1]))...
%         + ((ones(pSize,pSize,Isize)./pSize)... % jump occurred
%         .* permute(reshape(repmat(exp(-Ilog),pSize,pSize),Isize,pSize,pSize),[2 3 1]));
%             
else
    
    pp1gpIp1 = zeros(pSize, pSize, Isize);
    tmpp = zeros(pSize, pSize);
    for Ip1 = 1:Isize
        for p = 1:pSize
            for pp1 = 1:pSize
                % Beta(x; a,b)
                a = 1 + (exp(Ilog(Ip1)))*pvec(p);
                b = 1 + (exp(Ilog(Ip1)))*(1-pvec(p));
                x = pvec(pp1);
                if ~(x==0)||~(x==1)
                    logkerna = (a-1)*log(x);
                    logkernb = (b-1)*log(1-x);
                    betaln_ab = gammaln(a) + gammaln(b) - gammaln(a+b);
                    tmpp(pp1, p) = exp(logkerna + logkernb - betaln_ab);
                else
                    tmpp(pp1, p) = 0;
                end
                
            end
            tmpp(:,p) = tmpp(:,p)./sum(tmpp(:,p));
        end
        pp1gpIp1(:,:,Ip1) = tmpp; % place tmpp in it.
    end
end    
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
    IDist(trial,:) = sum(pI); % column-wise sum.
    IEst(trial,:)  = sum(IDist(trial,:).*Ilog);
    vEst(trial,:)  = sum(IDist(trial,:).*(-Ilog));    
    IEstExp(trial,:)  = sum(IDist(trial,:).*exp(Ilog));
    vEstExp(trial,:)  = sum(IDist(trial,:).*exp(-Ilog));
    
    
    % p(k)
    Ik = squeeze(sum(pIk,1));
    kDist(trial,:) = sum(Ik); % column-wise sum.
    kEst(trial,:)  = sum(kDist(trial,:).*klog);
    kEstExp(trial,:) = sum(kDist(trial,:).*exp(klog));
    
    
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
    
    %

 end



t=toc;


% make some plots
X = 1:length(reward);

figure
suptitle(str)
subplotNum = 6;
subplot(subplotNum,1,1)
plot (X(reward==0),IEst(reward==0),'.r',X(reward==1),IEst(reward==1),'.b') , title('IEst')
hold on
subplot(subplotNum,1,2)
plot(X(reward==0),1./IEst(reward==0),'.r',X(reward==1),1./IEst(reward==1),'.b')  , title('1./IEst')  % volatility
subplot(subplotNum,1,3)
plot(X(reward==0),IEstExp(reward==0),'.r',X(reward==1),IEstExp(reward==1),'.b')  , title('IEstExp')
subplot(subplotNum,1,4)
plot(X(reward==0),1./IEst(reward==0),'.r',X(reward==1),1./IEst(reward==1),'.b') , title('1./IEstExp')
subplot(subplotNum,1,5)
plot(X(reward==0),vEst(reward==0),'.r',X(reward==1),vEst(reward==1),'.b') , title('vEst')
subplot(subplotNum,1,6)
plot(X(reward==0),vEstExp(reward==0),'.r',X(reward==1),vEstExp(reward==1),'.b') , title('vEstExp (the important one)')

figure
suptitle(str)
hold on
subplot(2,1,1)
plot (kEst) , title('kEst')
subplot(2,1,2)
plot(kEstExp), title('kEstExp')


probSched = [ones(1,120)'*.75; ones(1,40)'*.20; ones(1,40)'*.80; ones(1,40)'*.20; ones(1,40)'*.80]; % experimental setup
%probSchedAdvice = [ones(1,30)'*.75; ones(1,10)'*.20; ones(1,10)'*.80; ones(1,10)'*.20; ones(1,10)'*.80; ones(1,50)'*.15];
figure, plot(X(reward==0),pEst(reward==0),'.r',X(reward==1),pEst(reward==1),'.b') % reward rate
suptitle(str)
hold on, plot(probSched,'r')
ylim([0 1])




