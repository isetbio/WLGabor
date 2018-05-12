function [absorptions, cm, emPath] = ccAbsorptions(ois, nTrials, varargin)
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


%%  Parse inputs

p = inputParser;
p.addRequired('ois', @(x)isa(x, 'oiSequence'));
p.addParameter('nTrials', nTrials, @isnumeric);

p.parse(ois,varargin{:})

ois     = p.Results.ois;
nTrials = p.Results.nTrials;

%% Initialize the cone mosaic
%% Make the cone mosaic

% cone mosaic params
% fovDegs = 1;
% TimeIntegrat = 1/100;
% resamplingFactors = 3; %[1 3 6 13];
%
% % Instantiate a hex mosaic with a specific resampling factor
%
% cm = coneMosaicHex(resamplingFactors, 'fovDegs', fovDegs);
% cm.integrationTime = TimeIntegrat;
% legends{numel(legends)+1} = sprintf('resampling: %2.0f ms', resamplingFactors(iSampleIndex));

% Create the object
cm = coneMosaic;

% Match the oisequence timing
cm.integrationTime = ois.timeStep;

% Make the cm smaller than the oi size, but never smaller than 0.2 deg
fovDegs = max(oiGet(ois.oiFixed,'fov') - 0.2, 0.2);  % Degrees
cm.setSizeToFOV(fovDegs);

%% Deal with eye movements
if (numel(varargin) >= 1)
    emPath = varargin{1};
else
    
    % Generate eye movement data for 2 1-second long trials, with a
    % sample time of 1 millisecond.  We want to be able to control the
    % size of the eye movements, making sure that they aren't too big so
    % the stimulus is not visible on the mosaic.
    emDuration = ois.timeStep*numel(ois.timeAxis);   % This should be get/set type stuff
    sampTime   = ois.timeStep;
    
    % Do not compute velocity of eye movements
    compVelocity = false;
    usePar       = false;
    
    % Random seed to be used in all eye movement compute() calls
    ranSeed = 1;
    
    % Initialize object
    fixEMobj = fixationalEM();
    
    % First case: micro-saccades, and drift
    % Set all params to their default value
    fixEMobj.setDefaultParams();
    fixEMobj.microSaccadeType = 'stats based';
    fixEMobj.randomSeed = ranSeed;
    
    fixEMobj.compute(emDuration, sampTime, nTrials, ...
        compVelocity, 'useParfor', usePar);
    %% Compute the number of eye movements for this integration time and oiSequence
    
    eyeMovementsPerTrial = ois.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
    % We needed the cone mosaic to calculate an example eye movement path
    fixEMobj.computeForConeMosaic(cm, eyeMovementsPerTrial, ...
        'nTrials', nTrials, ...
        'computeVelocity', true, ...
        'rSeed', 1);
    %%  Invoke cone compute with eye movements
    
    emPath  = fixEMobj.emPos;
end

% Compute the absorptions for all trials
absorptions = cm.compute(ois, 'empath', emPath);

end