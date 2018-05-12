function [absorptionsStim, absorptionsNoise] = csfAbsorptions(varargin)
%
%   We should probably call csfStimuli() to get the ois
%   Then we call csfAbsorptions(ois) to get the absorptions
%
% 
% ZL/BW

%% Create 100 samples of a large signal plus noise trials
nTimeSteps = 20;
integrationTime = 0.005;
sampleTimes = ((1:nTimeSteps) - 1) * integrationTime;   % Five ms integration time
nTrials    = 100;

% Make the harmonic with some contrast optical image
clear hparams
hparams(2)           = harmonicP;
hparams(2).freq      = sFreq;       % Cycles per field of view
hparams(2).GaborFlag = 0;
hparams(2).contrast  = sContrast;   % High contrast

% Make the constant part
hparams(1)           = hparams(2);
hparams(1).contrast  = 0;

% Make side parame6ters
sparams.fov          = fov;

% Temporal weights for the stimulus
stimWeights = ones(1, nTimeSteps);

% tSD = 30;
% stimWeights = ieScale(fspecial('gaussian',[1,nTimeSteps],tSD),0,1);

% Summarize
fprintf('Cycles per degree %.1f\n',hparams(1).freq/sparams.fov);

%%
ois = oisCreate('harmonic','blend',stimWeights, ...
    'testParameters',hparams,'sceneParameters',sparams, ...
    'sampleTimes',sampleTimes);

% nTrials, row, col, time
[absorptionsStim, cmStim, emPath] = ccAbsorptions(ois, nTrials);
% meanStim = mean(absorptionsStim,4);
% meanStim = permute(meanStim,[2 3 1]);

%{
vcNewGraphWin; imagesc(squeeze(meanStim(:,:,5)))
ieMovie(meanStim);
%}

%% Create 100 samples of noise only trials

hparams(2).contrast  = 0;   % High contrast
ois = oisCreate('harmonic','blend',stimWeights, ...
    'testParameters',hparams,'sceneParameters',sparams, ...
    'sampleTimes',sampleTimes);

% nTrials, row, col, time
[absorptionsNoise, cmNoise, emPath] = ccAbsorptions(ois, nTrials);
% meanNoise = mean(absorptionsNoise,4);
% meanNoise = permute(meanNoise,[2 3 1]);

%{
vcNewGraphWin; imagesc(squeeze(meanNoise(:,:,5)))
ieMovie(meanNoise);
%}
