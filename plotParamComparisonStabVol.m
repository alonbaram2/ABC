function plotParamComparisonStabVol (evids, schedules, schedulesParams,scores, highProbToPlot)

i=find(schedulesParams.highProbCand==highProbToPlot);
figure('position',[50,50,1200,700])
suptitle(sprintf('mean (joint - separate) model evidence for high Prob of reward = %g',schedulesParams.highProbCand(i)));
subplot (2,3,1)
hold on
title('AC 1st half (corr.)')
imagesc(schedulesParams.flipStabCand,schedulesParams.flipVolCand,scores.meanDiffEvidAC1(:,:,i))
colorbar
ylabel('flipVol')

subplot (2,3,2)
hold on
title('AC 2nd half(uncorr)')
imagesc(schedulesParams.flipStabCand,schedulesParams.flipVolCand,scores.meanDiffEvidAC2(:,:,i))
colorbar
subplot (2,3,3)
hold on
title('AC: 1st - 2nd (corr-uncorr)' )
imagesc(schedulesParams.flipStabCand,schedulesParams.flipVolCand,scores.meanDiffEvidAC1(:,:,i) - scores.meanDiffEvidAC2(:,:,i))
colorbar
subplot (2,3,4)
hold on
title('AB 1st half (uncorr)')
imagesc(schedulesParams.flipStabCand,schedulesParams.flipVolCand,scores.meanDiffEvidAB1(:,:,i))
colorbar
ylabel('flipVol')
xlabel('flipStab')
subplot (2,3,5)
hold on
title('AB 2nd half (corr)')
imagesc(schedulesParams.flipStabCand,schedulesParams.flipVolCand,scores.meanDiffEvidAB2(:,:,i))
colorbar
xlabel('flipStab')
subplot (2,3,6)
hold on
title('AB: 2nd-1st (corr-uncorr)')
imagesc(schedulesParams.flipStabCand,schedulesParams.flipVolCand,scores.meanDiffEvidAB2(:,:,i)-scores.meanDiffEvidAB1(:,:,i))
colorbar
xlabel('flipStab')

