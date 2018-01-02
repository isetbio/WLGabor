%% Cone mosaic - Gabor stimuli response
%

%% A monochrome Gabor patch
clear hparams
% The modulating harmonic parameters.  The possibilities are defined in the
% sceneCreate('harmonic', ...) function, and the complete list of
% parameters is returned by harmonicP.
hparams(2) = harmonicP; 
hparams(2).freq = 6; hparams(2).GaborFlag = .2;

% The matched, zero contrast, harmonic parameters
hparams(1) = hparams(2); hparams(1).contrast = 0;

% The general scene properties can also be set.  Here we make the scene 1
% deg of visual angle
sparams.fov = 1;

% And then we make a Gaussian temporal modulation that brings the stimulus
% on and off
stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
integrationTime = 0.050;
sampleTimes = (1:length(stimWeights))*integrationTime;

% The two harmonics are 'blended', which means at each moment in time we
% have a weighted sum of the two where the weights sum to 1.
ois = oisCreate('harmonic','blend',stimWeights, ...
    'sampleTimes', sampleTimes,...
    'testParameters',hparams,...
    'sceneParameters',sparams);
% ois.visualize('movie illuminance');

%%
fov = 1;
cm = coneMosaic;
cm.integrationTime = integrationTime;
cm.setSizeToFOV(fov);
cm.emGenSequence(ois.length);
cm.compute(ois);
cm.window;

%% A color Gabor patch
dsp        = displayCreate('LCD-Apple.mat');
wave       = displayGet(dsp,'wave');
backSPD    = displayGet(dsp,'spd primaries')*0.5*ones(3,1);
backSPD    = Energy2Quanta(wave,backSPD);
[~,modSPD] = humanConeIsolating(dsp);
modSPD     = Energy2Quanta(wave,modSPD);

clear hparams
hparams(2)           = harmonicP;
hparams(2).freq      = 4; 
hparams(2).GaborFlag = 0.2;
hparams(2).ang       = pi/6;

hparams(2).backSPD   = backSPD;
hparams(2).modSPD    = modSPD*[3 3 1]';  % Could be a mixture of the modSPDs
hparams(2).wave      = wave;
hparams(1)           = hparams(2);
hparams(1).contrast  = 0;

% Fifty time steps.  Each will be one integration time.
stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
integrationTime = 0.005;
sampleTimes = (1:length(stimWeights))*integrationTime;
ois = oisCreate('harmonic','blend',stimWeights, ...
    'sampleTimes', sampleTimes,...
    'testParameters',hparams,...
    'sceneParameters',sparams);

% ois.visualize('movie rgb');
% ois.visualize('weights')
% humanConeContrast(hparams(2).modSPD,hparams(2).backSPD,hparams(2).wave,'energy')

% vname = fullfile(isetbioRootPath,'local','oisVideo2.mp4');
% ois.visualize('movie rgb','vname',vname);


%%
fov = 1;
cm = coneMosaic;
cm.integrationTime = integrationTime;
cm.setSizeToFOV(fov);
cm.emGenSequence(ois.length);
cm.compute(ois);
cm.window;
