function [PC, samples] = csfPC(sFreq,sContrast, fov, nPCs)
% Calculate the principal components from signal and noise trials
%
%  We pre-compute the likely stimuli,  The logic is to use a high
%  contrast stimulus and make many samples of the mean absoprtion
%  pattern.  Then use a noise stimulus and make many samples.  Then
%  put all the means (noise and signal) into a big matrix.  Compute
%  the SVD to get the top N principal components.  When there are only
%  two classes and no eye movements we probably only need two
%  components. If we have eye movements, then we will probably need a
%  few more principal components.
%
%  We use these components to do CSF calculation.
%  
%  After checking the PCs with different contrast level, we decide to
%  choose to use the PCs with the contrast = 1 (for now), and calculate the
%  weight 
%
%  When there are no eye movements, 
%
%% ZL/BW

%{

sFreq = 4;
fov  = 1;
nPCs = 200;

% sFreq, sContrast, fov, nPCs

PC1 = csfPC(sFreq,1.0,fov,nPCs);
vcNewGraphWin; imagesc(reshape(PC1(:,1),cmStim.rows,cmStim.cols))
PC2 = csfPC(7,0.01,fov,nPCs);
vcNewGraphWin; imagesc(reshape(PC2(:,1),cmStim.rows,cmStim.cols))

vcNewGraphWin; plot(PC1(:,1),PC2(:,1),'.'); 
grid on; axis equal; identityLine;

PC1(:,1)'*PC2(:,1)

%}

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
meanStim = mean(absorptionsStim,4);
meanStim = permute(meanStim,[2 3 1]);

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
meanNoise = mean(absorptionsNoise,4);
meanNoise = permute(meanNoise,[2 3 1]);

%{
vcNewGraphWin; imagesc(squeeze(meanNoise(:,:,5)))
ieMovie(meanNoise);
%}


%% Put the mean absorptions from all trials into one big matrix and calculate the PCs

samples = [RGB2XWFormat(meanStim), RGB2XWFormat(meanNoise)];

[U, ~, ~] = svd(samples,'econ');
    %{
      nPCs = 50;
      [A, B, C] = svd(samples, 'econ');
      wgts = A' * samples;
      wgts_two = wgts(1:nPCs,:);
      
    %}
PC = U;

    %{
        sample_rec = PC * wgts_two;
        nframe = 110;
        vcNewGraphWin;imagesc(reshape(sample_rec(:,nframe),cmStim.rows,cmStim.cols))
        vcNewGraphWin;imagesc(reshape(samples(:,nframe),cmStim.rows,cmStim.cols))
    %}
% sValues = diag(S);
% vcNewGraphWin; semilogy(sValues(1:4));
% vcNewGraphWin; imagesc(reshape(PC(:,1),cmStim.rows,cmStim.cols))
% vcNewGraphWin; imagesc(reshape(PC(:,2),cmStim.rows,cmStim.cols))

end