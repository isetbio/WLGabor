function absorptions = ccAbsorptions(stmlType, nTrials)
% Calculate the total cone absorptions for multiple trials  
%
% Description:
%  Create a harmonic stimulus, a cone mosaic, and calculate the total
%  number of cone absorptions from the stimulus.
%
%
% Inputs:
%
% Optional key/value pairs:
%   emDuration         - total simulation time (s)
%   sampTime           - integrate time unit (s)
%   nTrials            - number of trials
%   moduleLength       - change according to sampTime
%   fovDegs            - Set accordingly
%   resamplingFactors  - Cone mosaic for hex case
%
% Outputs:
%   absorptions - matrix of numbers for the cone mosaic absorptions
%
% ZL, SCIEN, 2018
%
% Based on t_fixationalEyeMovementsTypes
%
% See also:


%%  ZL to implement the parameters

%p = inputParser;

%% Generate eye movements

% Generate eye movement data for 2 1-second long trials, with a
% sample time of 1 millisecond.
emDuration = 1; sampTime = 1/100;

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
fixEMobj.microSaccadeType = 'stats based';
fixEMobj.randomSeed = ranSeed;

fixEMobj.compute(emDuration, sampTime, nTrials, ...
    compVelocity, 'useParfor', usePar);

%%  Generate oisequence
% Two scenes, for oiFixed and oiModulated
scene = cell(1,2);
% oiModulated harmonic parameters
tparams(2) = harmonicP;
tparams(2).freq = 4;
tparams(2).GaborFlag = 0.2;

% oiFixed has the same parameters, but zero contrast
tparams(1) = tparams(2);
tparams(1).contrast = 0;

if (strcmp(stmlType, 'No'))   % No stimulus
    tparams(2).contrast = 0;
end
% Create the harmonic scenes
for ii=1:2
    scene{ii} = sceneCreate('harmonic',tparams(ii));
end


% Compute optical images from the scene
OIs = cell(1, 2);
oi = oiCreate;
for ii = 1:2
    OIs{ii} = oiCompute(oi,scene{ii});
end

integrateTimeOI = 0.01;
moduleLength = emDuration / integrateTimeOI;
modulation = ieScale(fspecial('gaussian',[1,moduleLength],10),0,.5);
sampleTimes = ((1: length(modulation))-1)*integrateTimeOI;
ois = oiSequence(OIs{1}, OIs{2}, sampleTimes, modulation, ...
    'composition', 'blend');
%ois.visualize('movie illuminance');
%%
fovDegs = 1;
TimeIntegrat = 1/100;

cm = coneMosaic('fovDegs', fovDegs);
cm.integrationTime = TimeIntegrat;

%% Generate cone mosaic
% cone mosaic parms
% fovDegs = 1;
% TimeIntegrat = 1/100;
% resamplingFactors = 3; %[1 3 6 13];
%
% % % Instantiate a hex mosaic with a specific resampling factor
% cm = coneMosaicHex(resamplingFactors, 'fovDegs', fovDegs);
% cm.integrationTime = TimeIntegrat;
%legends{numel(legends)+1} = sprintf('resampling: %2.0f ms', resamplingFactors(iSampleIndex));

%% Compute the number of eye movements for this integration time and oiSequence
eyeMovementsPerTrial = ois.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
fixEMobj.computeForConeMosaic(cm, eyeMovementsPerTrial, ...
    'nTrials', nTrials, ...
    'computeVelocity', true, ...
    'rSeed', 1);
%%
cm.emPositions = fixEMobj.emPos;
[absorptions, ~, ~, ~ ] = cm.compute(ois);
cm.emPositions = fixEMobj.emPos;

%% Might need to change the window function
%cm.window;
%% Extract  based on number of frame
absorptions = absorptions;

end