% Explore the svm classifier on the PC components for harmonic
% stimuli.
% The steps comes
% 
%
% ZL

%% Parameter initialization
digits(10);
sFreq         = 6; % logspace(0, 1.5, 1);
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
%% Calculate the PCs for the high contrast example of the stimulus

cm.noiseFlag = 'none';
template     = squeeze(cm.compute(ois));
% vcNewGraphWin; imagesc(template); colormap(gray);

%Calculate the PCs (the whole sets of principal components)
PCs = csfPC(template, nTrials);

%% Create the test stimulus at a lower contrast level

hparams(2).contrast  = 0.1;
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

%%
hparams(2).contrast  = 0.0;
ois = oisCreate('harmonic', 'blend', stimWeights, ...
    'testParameters', hparams, 'sceneParameters', sparams);
% ois.visualize('movie illuminance');

% Eye movements?
absorptionsNoise = cm.compute(ois, 'empath', empath);
meanabsorptionsNoise = mean(absorptionsNoise, 4);

%% Calculate the wgts for the Noise
wgtsNoise = csfWgts(meanabsorptionsNoise, PCs, nPC);

% Let's plot the weights for the stimulus and no stimulus conditions

%% Process SVM

correctness   = zeros(numel(Contrast), numel(Freq));

stmlType = {'Yes', 'No'};
classStmls = cell(size(wgts{1, 1}, 2),1);

% Put a comment here
for i = 1 : size(wgts{1, 1}, 2) / 2
    classStmls{i} = stmlType{1};
    classStmls{i + size(wgts{1, 1}, 2) / 2} = stmlType{2};
end

for c = 1 : size(Contrast, 2)
    for f = 1 : size(Freq, 2)
        dataStmls = wgts{c, f}';
        meanCorrect = svmProcess(dataStmls, classStmls);
        correctness(c, f) = meanCorrect;
    end
end

%%