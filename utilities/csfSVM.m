% Explore the svm classifier on the PC components for harmonic
% stimuli.
% The steps comes
% 
%
% ZL

%% Parameter initialization
sFreq         = 16; % logspace(0, 1.5, 1);
nPCs          = 2;
fov           = 1;
sContrast     = 1;

% Set up the stimulus parameters
clear hparams
hparams(2)           = harmonicP;
hparams(2).freq      = sFreq;
hparams(2).contrast  = sContrast;
hparams(1)           = hparams(2);
hparams(1).contrast  = 0;

sparams.fov = 1;

nTimeSteps = 20;
stimWeights = ones(1, nTimeSteps);

%% Create the OIS for the PC calculation

% Stimulus for the PC
ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);
% ois.visualize('movie illuminance');

%% Set up cone mosaic parameters

integrationTime = 0.005;
sampleTimes = ((1:nTimeSteps) - 1) * integrationTime;   % Five ms integration time
nTrials    = 100;

cm = coneMosaic;
cm.integrationTime = ois.timeStep;

% Make the cm smaller than the oi size, but never smaller than 0.2 deg
fovDegs = max(oiGet(ois.oiFixed,'fov') - 0.2, 0.2);  % Degrees
cm.setSizeToFOV(fovDegs);

% EM path is set to be zero (meaning no eyemovement for now). Will 
% implement the fixEM in the future.
empath = zeros(nTrials, nTimeSteps, 2);

%% Calculate the absorption template for the high contrast example of the stimulus

cm.noiseFlag = 'none';
templateHighContrast     = mean(squeeze(cm.compute(ois)), 3);
% vcNewGraphWin; imagesc(template); colormap(gray);

%% Calculate the absorption template for zero contrast of the stimulus
hparams(2).contrast  = 0.0;
ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);
templateZeroContrast     = mean(squeeze(cm.compute(ois)), 3);

%% Calculate the PCs (the whole sets of principal components)
PCs = csfPC(templateHighContrast, templateZeroContrast, nTrials);

%% Create the test stimulus at a lower contrast level
cm.noiseFlag = 'random';

hparams(2).contrast  = 0.001;
ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);
% ois.visualize('movie illuminance');

% Eye movements?

absorptions = cm.compute(ois, 'empath', empath);
meanabsorptionsStimulus = mean(absorptions, 4);
%{
    thisTrial = 10;
    thisFrame = 5;
    vcNewGraphWin; imagesc(squeeze(absorptions(thisTrial, :, :, thisFrame))); colormap(gray);

%}

%% Set up for PCs parameters
nPC = 2;
%% Calculate the wgts for the absorptions

wgtsStimulus = csfWgts(meanabsorptionsStimulus, PCs, nPC);

%% Create a "blank" pattern without stimulus
hparams(2).contrast  = 0.0;
ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);
% ois.visualize('movie illuminance');

% No eyemovements for now.
absorptionsNoise = cm.compute(ois, 'empath', empath);
meanabsorptionsNoise = mean(absorptionsNoise, 4);

%% Calculate the wgts for the Noise
wgtsNoise = csfWgts(meanabsorptionsNoise, PCs, nPC);

% Let's plot the weights for the stimulus and no stimulus conditions

%% Process SVM

% correctness   = zeros(numel(Contrast), numel(Freq)); will be implemented
% in the future.

classStmls = [ones(size(wgtsStimulus, 2), 1); zeros(size(wgtsNoise, 2), 1)];

dataStmls = [wgtsStimulus'; wgtsNoise'];
meanCorrect = svmProcess(dataStmls, classStmls);


%%