function t_ccDiscriminate
%% cc = contrast classification
%
%   Create Gabor stimuli with different contrasts.  Run an SVM to see
%   if we can tell them apart.
%
% ZL, SCIEN STANFORD, 2018

%%
ieInit

%%
nFrames   = 100;   % Number of temporal samples
numFrameNoStmls = 100;
numFrameTotal = nFrames + numFrameNoStmls;
nTrials = 100;
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
classStmls = cell(numFrameTotal,1);
for i = 1 : nFrames
    classStmls{i} = stmlType{1};
    classStmls{i + nFrames} = stmlType{2};
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