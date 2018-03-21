function t_principleComponent
% The function is a short demo to give the idea that how to generate the 
% principleComponent for the cone absorptions of each trial. It'll first
% generate the absorptions, vectorize every single frame, and calculate the
% Singular value decomposition


%% Temporary parameters

freq = 10;
contrast = 0.8;
%%  Generate oisequence

clear hparams

% Make the time varying part
hparams(2) = harmonicP;
hparams(2).freq      = freq;     % Cycles per field of view
hparams(2).GaborFlag = 0.2;
hparams(2).contrast  = contrast;

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
ois = oisCreate('harmonic','blend',stimWeights, ...
    'testParameters',hparams,'sceneParameters',sparams);
%% Generate absorption for single trail

nTrials = 1;
[absorptions, cm] = ccAbsorptions(ois, nTrials);

%% Vectorize the absorption
absorptionVec = permute(squeeze(absorptions),[1 2 3]);
absorptionVec = RGB2XWFormat(absorptionVec)';

%% Calculate the svd
[U, S, V] = svd(absorptionVec);
end