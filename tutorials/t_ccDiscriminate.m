function t_ccDiscriminate
%% cc = contrast classification
%
%   Create Gabor stimuli with different contrasts.  Run an SVM to see
%   if we can tell them apart.
%   
%   Currently gives 99.75% correctness.
%
% ZL, SCIEN STANFORD, 2018
clc, close all, clear all;
%%
ieInit

%%
nFrames   = 100;   % Number of temporal samples
numFrameNoStmls = 100;

nTrials = 50;
stmlType = {'Yes', 'No'};   % Stimulus or no stimulus

%% Calculate the total number of absorptions

% <trials,row,col,time>
absorptions = ccAbsorptions(stmlType{1}, nTrials);
%{
% Look at the movie from a trial
thisTrial = 1;
trialData = squeeze(absorptions(thisTrial,:,:,:));
ieMovie(trialData);
%}

%% Create a vector of the mean absorptions on each trial
%
% We place the data into one big matrix with <space, trials>
meanAbsorptions = mean(absorptions, 4);
%{
trialData = squeeze(meanAbsorptions(1,:,:));
vcNewGraphWin;
imagesc(trialData); colormap(gray); colorbar; axis image
%}

% We want to put each vector of the cone absorptions in a trial into a
% row of the stimulus matrix. Later we will add the no stimulus trials with their label.
frameStmlsReshp = permute(meanAbsorptions,[2 3 1]);
frameStmlsReshp = RGB2XWFormat(frameStmlsReshp)';
% size(frameStmlsReshp)

%%  ZL

absorptionsNoStml = ccAbsorptions(stmlType{2}, nTrials);

%%
meanAbsorptionsNoStmls = mean(absorptionsNoStml, 4);


frameNoStmlsReshp = permute(meanAbsorptionsNoStmls, [2 3 1]);
frameNoStmlsReshp = RGB2XWFormat(frameNoStmlsReshp)';
%% plot mean results
meanStmlsPlot(meanAbsorptions);

hold on;

meanStmlsPlot(meanAbsorptionsNoStmls);

%%
dataStmls = [frameStmlsReshp;frameNoStmlsReshp];
classStmls = cell(2 * nTrials,1);
for i = 1 : nTrials
    classStmls{i} = stmlType{1};
    classStmls{i + nTrials} = stmlType{2};
end


%% parameter optimization
%{
%svm = fitcsvm(dataStmls, classStmls, 'OptimizeHyperparameters', 'all', 'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'));
%}


%% optimized parameter try
% best parameter for now
boxConstraint = 0.0017183;
kernelScale = 23.443;
kernelFunction = 'gaussian';
standarize = true;

kFold = 10;
svmOptimize = fitcsvm(dataStmls, classStmls, 'BoxConstraint', boxConstraint, 'KernelScale', kernelScale, 'KernelFunction', kernelFunction,'Standardize', standarize);
CVSVMOptimize = crossval(svmOptimize,'KFold',kFold);
percentCorrect = 1 - kfoldLoss(CVSVMOptimize,'lossfun','classiferror','mode','individual');
stdErr = std(percentCorrect)/sqrt(kFold);
meanPercent = mean(percentCorrect)

end