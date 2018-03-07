function t_fixationalEMConeSamp
% Editted version of t_fixationalEyeMovementsTypes
% Log
%   02/16/2018 Edit
%
%

%%
ieInit
clc, clear all, close all;
%% Generate eye movements

% Generate eye movement data for 2048 1-second long trials, with a
% sample time of 1 millisecond.
emDuration = 1.0; sampTime = 1/1000; nTrials = 2;

% Do not compute velocity of eye movements
compVelocity = false;
usePar = false;

% Random seed to be used in all eye movement compute() calls
ranSeed = 1;

% Initialize object
fixEMobj = fixationalEM();

% First case: No micro-saccades, only drift
% Set all params to their default value
fixEMobj.setDefaultParams();
fixEMobj.microSaccadeType = 'none';
fixEMobj.randomSeed = ranSeed;

fixEMobj.compute(emDuration, sampTime, nTrials, ...
    compVelocity, 'useParfor', usePar);

%{
visualizedTrial = 1;
colors = [1 0 0; 0 0 1; 0 0 0; 0.4 0.4 0.4];
legends = {};
intTimeIndex = 1;

vcNewGraphWin;
plot(fixEMobj.timeAxis*1000, squeeze(fixEMobj.emPosArcMin(visualizedTrial,:,1)), 's-', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 4, 'MarkerFaceColor', [0.8 0.8 0.8], ...
            'Color', squeeze(colors(intTimeIndex,:)));
    
legend(legends);
xlabel('time (ms)'); ylabel('x-position (arc min)');
grid on; set(gca, 'FontSize', 14);
%}

%%  What is the time step here?
clear hparams
hparams(2) = harmonicP; hparams(2).freq = 6; hparams(2).GaborFlag = .2;
hparams(1) = hparams(2); hparams(1).contrast = 0;
sparams.fov = 1;
stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
ois = oisCreate('harmonic','blend',stimWeights, ...
    'testParameters',hparams,'sceneParameters',sparams);
ois.visualize('movie illuminance');

%%  Generate oisequence
%{
% Two scenes, for oiFixed and oiModulated
scene = cell(1,2);
% oiModulated harmonic parameters
tparams(2) = harmonicP;
tparams(2).freq = 4;
tparams(2).GaborFlag = 0.2;

% oiFixed has the same parameters, but zero contrast
tparams(1) = tparams(2);
tparams(1).contrast = 0;

% Create the harmonic scenes
for ii=1:2
    scene{ii} = sceneCreate('harmonic',tparams(ii)); %% KEY TO KNOW HOW sceneCreate() function works.
end


% Compute optical images from the scene
OIs = cell(1, 2);
oi = oiCreate;
for ii = 1:2
    OIs{ii} = oiCompute(oi,scene{ii});
end

integrateTimeOI = 0.001;

modulation = ieScale(fspecial('gaussian',[1,100],10),0,.5);
sampleTimes = ((1:length(modulation))-1)*integrateTimeOI;
ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
    'composition', 'blend');
ois.visualize('movie illuminance');

%}
%%  Rectangular mosaic.
fovDegs = 1;
TimeIntegrat = 1/1000;

% Instantiate a hex mosaic with a specific resampling factor
cm = coneMosaic('fovDegs', fovDegs);
cm.integrationTime = TimeIntegrat;

%% Compute the number of eye movements 
%
% Accounts for this integration time, oiSequence and cone mosaic
eyeMovementsPerTrial = ois.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
fixEMobj.computeForConeMosaic(cm, eyeMovementsPerTrial, ...
    'nTrials', nTrials, ...
    'computeVelocity', true, ...
    'rSeed', 1);
cm.emPositions = fixEMobj.emPos;
cm.compute(ois);
cm.window;

%% Hex compute

% cone mosaic parms
fovDegs = 1;
TimeIntegrat = 1/1000;
resamplingFactors = 1; %[1 3 6 13];

% Instantiate a hex mosaic with a specific resampling factor
cm = coneMosaicHex(resamplingFactors, 'fovDegs', fovDegs);
cm.integrationTime = TimeIntegrat;
% legends{numel(legends)+1} = sprintf('resampling: %2.0f ms', resamplingFactors(iSampleIndex));

eyeMovementsPerTrial = ois.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
fixEMobj.computeForConeMosaic(cm, eyeMovementsPerTrial, ...
    'nTrials', nTrials, ...
    'computeVelocity', true, ...
    'rSeed', 1);

cm.emPositions = fixEMobj.emPos;
cm.compute(ois);
cm.window;


