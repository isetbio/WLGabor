%% Test eye movement and cone mosaic calculation

%% Initialize scene, oi, cone mosaic
scene = sceneCreate('rings rays');
scene = sceneSet(scene,'fov',2);
% ieAddObject(scene); sceneWindow;

oi = oiCreate;
oi = oiCompute(oi,scene);
% ieAddObject(oi); oiWindow;

%% Make the cone mosaic and one trial worth of eye movements

cmosaic = coneMosaic;
cmosaic.emGenSequence(10);
cmosaic.compute(oi);

%{
cmosaic.window;
%}

%% Now try multiple trials

cmosaic = coneMosaic;

% Generate the eye movements
nTrials = 100; nEyeMovements = 50;
emPaths = zeros(nTrials, nEyeMovements, 2);
% Generate the eye movements for each trial
for iTrial = 1:nTrials
    theEMPaths(iTrial , :,:) = cMosaic.emGenSequence(nEyeMovements);
end

% Build an oi sequence
hparams(2) = harmonicP; hparams(2).freq = 6; hparams(2).GaborFlag = .2;
hparams(1) = hparams(2); hparams(1).contrast = 0;
sparams.fov = 1;
stimWeights = ieScale(fspecial('gaussian',[1,50],15),0,1);
ois = oisCreate('harmonic','blend',stimWeights, ...
    'testParameters',hparams,'sceneParameters',sparams);

% Check that we get the absorptions back
absorptions = cMosaic.compute(ois, 'emPaths', theEMPaths);

%%