function singleMeanCorrect = t_ccDiscriminate(spatialF, sContrast, fov)
%% Calculate the probability correct detection for a Gabor 
%
% Create a Gabor stimulus with the specified contrast.  Run an SVM to see
% if it can be discriminated from a zero contrast stimulus with matched eye
% movements
%
%
% Input parameters
%
% Optional key/value pairs
%
% Return
%
% ZL, SCIEN STANFORD, 2018

%{

%}

%%
ieInit
% clc, close all, clear all;

%% Make the stimulus

%%  Generate oisequence

clear hparams


% Make the time varying part
hparams(2) = harmonicP;
hparams(2).freq      = spatialF;     % Cycles per field of view
hparams(2).GaborFlag = 0.2;
hparams(2).contrast  = sContrast;

% Make the constant part
hparams(1) = hparams(2);
hparams(1).contrast = 0;
sparams.fov = fov;
fprintf('Cycles per degree %.1f\n',hparams(1).freq/sparams.fov);

% These are the scalar over time for the oi sequence
nTimeSteps = 100;
tSD = 30;
stimWeights = ieScale(fspecial('gaussian',[1,nTimeSteps],tSD),0,1);

% Build the sequence
fprintf('Stimulus sequence (freq, contrast): %.1f, %.2f\n',spatialF, sContrast);
ois = oisCreate('harmonic','blend',stimWeights, ...
    'testParameters',hparams,'sceneParameters',sparams);

%{
 ois.visualize('movie illuminance');
%}
    
nTrials = 100;      % Debug with a small number of trials
stmlType = {'Yes', 'No'};   % Stimulus or no stimulus

%% Calculate the total number of absorptions

% absorptionsStim:  <trials,row,col,time>
[absorptionsStim, cm] = ccAbsorptions(ois, nTrials);
% cm.window;

%{
% Look at the movie from a trial
thisTrial = 10;
trialData = squeeze(absorptionsStim(thisTrial,:,:,:));
ieMovie(trialData);
%}

%% Create the zero contrast version

fprintf('Creating noise stimulus\n');
hparams(2).contrast = 0;
ois = oisCreate('harmonic','blend',stimWeights, ...
    'testParameters',hparams,'sceneParameters',sparams);
absorptionsNostim = ccAbsorptions(ois, nTrials);

%% Create a vector of the mean absorptions on each trial
%

% This could be a function that converts absorptions into the right
% format for classification

% We place the data from each temporal trial into one temporal matrix
% meanAbsorptionsStim:  <trials, row, col> 
meanAbsorptionsStim = mean(absorptionsStim, 4);
%{
trialData = squeeze(meanAbsorptionsStim(1,:,:));
vcNewGraphWin;
imagesc(trialData); colormap(gray); colorbar; axis image
%}

% We want to put each vector of the cone absorptions in a trial into a row
% of the stimulus matrix. Later we will combine this with the no stimulus
% trials with their label.  The original ordering is trial,row,col.  We
% shift to row,col,trial and then reshape
frameStmlsReshp = permute(meanAbsorptionsStim,[2 3 1]);
% ieMovie(frameStmlsReshp);

% <trials, row*col>
frameStmlsReshp = RGB2XWFormat(frameStmlsReshp)';
% size(frameStmlsReshp)

% Do it again for the no contrast case
frameNoStmlsReshp = permute(mean(absorptionsNostim, 4),[2 3 1]);
frameNoStmlsReshp = RGB2XWFormat(frameNoStmlsReshp)';
% size(frameNoStmlsReshp)

%% Combine the stimuli and the labels

dataStmls = [frameStmlsReshp; frameNoStmlsReshp];
% size(dataStmls)
classStmls = cell(2 * nTrials,1);
for i = 1 : nTrials
    classStmls{i} = stmlType{1};
    classStmls{i + nTrials} = stmlType{2};
end


%% Simple svm fit

svm = fitcsvm(dataStmls, classStmls);
kFold = 10;
CVSVMOptimize = crossval(svm,'KFold',kFold);
probabilityCorrect = 1 - kfoldLoss(CVSVMOptimize,'lossfun','classiferror','mode','individual');

fprintf('Mean probability correct %.2f\n',mean(probabilityCorrect));
singleMeanCorrect = mean(probabilityCorrect);

%% Curves to make

% 1. Frequency on x and probability correct, for a fixed contrast
% level.
% 2.  Sweep out many contrast levels for each frequency.  Show the
% contrast needed to achieve a fixed probability correct (75% or 80%)
% iso-performance level.
% 

% Have fun.

%%

%{
% Signal known approximately idea
% And using the whole time series.

We aren't using the whole time series.  We are using only the mean of
the time series.  Can we set this up to run with the time series?
One way is to take a 100% contrast example of the stimulus, and make a
time series.  S(x,y,t).  Then compute the principal components of the
images, S(x,y,:).  This would be

   [ ..,  vector(S(x,y,ii),... ]

Take the SVD of this matrix.  THe first 10 biggest singular values
would be the basis functions.  We replace the  images in the SVM with
the weights on those basis functions.  BW thinks it will be reasonably
accurate to use less than 5 basis functions for a known harmonic.

%}






% stdErr = std(probabilityCorrect)/sqrt(kFold);

% 
% 
% %% parameter optimization
% %{
% %svm = fitcsvm(dataStmls, classStmls, 'OptimizeHyperparameters', 'all', 'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
%     'expected-improvement-plus'));
% %}
% sigma = optimizableVariable('sigma',[1e-5,1e5],'Transform','log');
% box = optimizableVariable('box',[1e-5,1e5],'Transform','log');
% 
% %% optimized parameter try
% % best parameter for now
% boxConstraint = 0.0017183;
% kernelScale = 23.443;
% kernelFunction = 'gaussian';
% standarize = true;
% 
% kFold = 10;
% svm = fitcsvm(dataStmls, classStmls);
% %svmOptimize = fitcsvm(dataStmls, classStmls, 'BoxConstraint', boxConstraint, 'KernelScale', kernelScale, 'KernelFunction', kernelFunction,'Standardize', standarize, 'Cost', [0 100;100 0]);
% CVSVMOptimize = crossval(svm,'KFold',kFold);
% percentCorrect = 1 - kfoldLoss(CVSVMOptimize,'lossfun','classiferror','mode','individual');
% stdErr = std(percentCorrect)/sqrt(kFold);
% meanPercent = mean(percentCorrect)
% %%
% %% Bayesian
% c = cvpartition(2 * nTrials,'KFold',10);
% minfn = @(z)kfoldLoss(fitcsvm(dataStmls,classStmls,'CVPartition',c,...
%     'KernelFunction','rbf','BoxConstraint',z.box,...
%     'KernelScale',z.sigma));
% %%
% results = bayesopt(minfn,[sigma,box],'IsObjectiveDeterministic',true,...
%     'AcquisitionFunctionName','expected-improvement-plus')
% %%
% z(1) = results.XAtMinObjective.sigma;
% z(2) = results.XAtMinObjective.box;
% SVMModel = fitcsvm(dataStmls,classStmls,'KernelFunction','rbf',...
%     'KernelScale',z(1),'BoxConstraint',z(2));
% %%
% kFold = 10;
% CVSVMTest = crossval(SVMModel,'KFold',kFold);
% percentCorrect = 1 - kfoldLoss(CVSVMTest,'lossfun','classiferror','mode','individual');
% stdErr = std(percentCorrect)/sqrt(kFold);
% meanPercent = mean(percentCorrect)
end