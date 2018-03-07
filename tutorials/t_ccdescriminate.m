function t_ccdescriminate
%% cc = contrast classification
%%
clc, clear all, close all;
%%
ieInit
%%
numFrameStmls = 100;
numFrameNoStmls = 100;
numFrameTotal = numFrameStmls + numFrameNoStmls;
nTrials = 100;
%%
stmlType = {'Yes', 'No'};
%%
frameStmls = squeeze(ccframeStmlsCreate(numFrameStmls, stmlType{1}, nTrials));
[trials, row, col, numFrameStmls] = size(frameStmls);
meanFrameStmls = mean(frameStmls, 4);
frameStmlsReshp = reshape(meanFrameStmls, [row * col, trials]);
frameStmlsReshp = frameStmlsReshp';
%%
frameNoStmls = squeeze(ccframeStmlsCreate(numFrameNoStmls, stmlType{2}, nTrials));
[trials, row, col, numFrameNoStmls] = size(frameNoStmls);
meanFrameNoStmls = mean(frameNoStmls, 4);
frameNoStmlsReshp = reshape(meanFrameNoStmls, [row * col, trials]);
frameNoStmlsReshp = frameNoStmlsReshp';
%% plot mean results
meanStmlsPlot(meanFrameStmls);

hold on;

meanStmlsPlot(meanFrameNoStmls);

%%
dataStmls = [frameStmlsReshp;frameNoStmlsReshp];
classStmls = cell(numFrameTotal,1);
for i = 1 : numFrameStmls
    classStmls{i} = stmlType{1};
    classStmls{i + numFrameStmls} = stmlType{2};
end
%classStmls{1:numFrameStmls} = stmlType{1};
%classStmls{numFrameStmls + 1 : end} = stmlType{2};

%%
kernelFunction = 'rbf';
standardizeSVMpredictors = false;

kFold = 10;


%svm = fitcsvm(dataStmls, classStmls, 'KernelFunction',kernelFunction,'KernelScale','auto', 'standardize', standardizeSVMpredictors);
svm = fitcsvm(dataStmls, classStmls);
CVSVM = crossval(svm,'KFold',kFold);
percentCorrect = 1 - kfoldLoss(CVSVM,'lossfun','classiferror','mode','individual');
stdErr = std(percentCorrect)/sqrt(kFold);
percentCorrect = mean(percentCorrect);

%% Plot result
sv = svm.SupportVectors;
figure
gscatter(dataStmls(:,1),dataStmls(:,2),classStmls)
hold on
plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
legend('Harmonic Stimulus','No Harmonic Stimulus','Support Vector')
hold off
%%














end