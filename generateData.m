function [cues, outcomes, vol, probs ] = generateData (stateBlockDuration, totalTrials, trialsUntilFlipVol, trialsUntilFlipStab,highProb)


% inputs: 
% stateBlockDuration - # trials in a block of correlation (either AB or AC)
% totalTrials - # trials in all the session  
% trialsUntilFlipVol - # trials after which the probabilities of reward 
%                   will flip for cues with a current volatile environment
% trialsUntilFlipStab - same for cues with a current stable environment.
                          
cues = randi(3,1,totalTrials) - 1; % 0<=>A, 1<=>B, 2<=>C
outcomes = nan(1,totalTrials);
vol = nan(1,totalTrials);
probs = nan(1,totalTrials);

% Hardcoded: AC are correlated in the first block. In the first half of the
% block AC are stable (B is volatile), and in the second half vica versa.
% In the second block AB are correlated with AB being volatile in the first
% half of the block.

blockHalf=0;
counterVol = 0;
counterStab = 0;
lowProb = 1-highProb;
probVol = highProb; % probability of reward for cues currently with a volatile schedule.
probStab = lowProb; % probability of reward for cues currently with a stable schedule.

for i = 1:totalTrials
    lastBlockHalf=blockHalf;
	blockNumber = ceil(i/stateBlockDuration); % starts at 1
    blockHalf = round((i - ((blockNumber-1)*stateBlockDuration)) / stateBlockDuration) + 1 ; % 1 - first half, 2 - 2nd half
    if blockHalf ~= lastBlockHalf
        counterVol = 0;
        counterStab = 0;
        % set initial probabilities for each half of a block
        if (blockNumber==1 && blockHalf==1) 
            probVol = highProb; % probability of reward for cues currently with a volatile schedule.
            probStab = lowProb; % probability of reward for cues currently with a stable schedule.        
        elseif(blockNumber==1 && blockHalf==2)
            probVol = lowProb;
            probStab = lowProb;
        elseif(blockNumber==2 && blockHalf==1)
            probVol = highProb;
            probStab = lowProb;
        elseif(blockNumber==2 && blockHalf==2)
            probVol = highProb;
            probStab = lowProb;
        end
    end 
    
    if counterVol == trialsUntilFlipVol 
    probVol = 1 - probVol;
    counterVol = 0;
    end
    if counterStab == trialsUntilFlipStab
        probStab = 1 - probStab;
        counterStab = 0;
    end
    counterVol = counterVol + 1;
    counterStab = counterStab + 1; 
    
    if (blockNumber==1 && blockHalf==1) % AC corr & stab
        if (cues(i)==0 || cues(i) ==2)
            probs(i) = probStab;
            vol(i) = 1/trialsUntilFlipStab;
        else probs(i) = probVol;
            vol(i) = 1/trialsUntilFlipVol;
        end
    elseif (blockNumber==1 && blockHalf==2) % AC corr & vol
        if (cues(i)==0 || cues(i) ==2)
            probs(i) = probVol;
            vol(i) = 1/trialsUntilFlipVol;
        else probs(i) = probStab;
             vol(i) = 1/trialsUntilFlipStab; 
        end         
    elseif (blockNumber==2 && blockHalf==1) % AB corr & vol
        if (cues(i)== 0 || cues(i) == 1)
            probs(i) = probVol;
            vol(i) = 1/trialsUntilFlipVol;
        else probs(i) = probStab;
            vol(i) = 1/trialsUntilFlipStab;
        end 
    elseif (blockNumber==2 && blockHalf==2) % AB corr & stab
        if (cues(i)== 0 || cues(i) == 1)
            probs(i) = probStab;
            vol(i) = 1/trialsUntilFlipStab;
        else probs(i) = probVol;
             vol(i) = 1/trialsUntilFlipVol;
        end 
    end
    outcomes(i) = (rand()<probs(i));
end
