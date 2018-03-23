function t_principalComponent
%
% Generate the principal components for the cone absorptions across trials.
% 
%   Generate the absorptions
%   Vectorize every time step
%   Calculate the Singular value decomposition of these
%
%
% ZL, SCIEN Stanford, 2018

%% Temporary parameters

spatialF = 2;
spatialC = 0.8;
fov = 0.4;

%%  Generate oisequence

clear hparams

% Make the time varying part
hparams(2) = harmonicP;
hparams(2).freq      = spatialF;     % Cycles per field of view
hparams(2).GaborFlag = 0.2;
hparams(2).contrast  = spatialC;

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

%{
% Show the mean
meanAbs = mean(absorptionVec);
reshape(meanAbs,cm.rows,cm.cols);
imagesc(tmp); colormap(gray); axis image
%}
absorptionVec0 = absorptionVec - meanAbs;
%{
z = mean(absorptionVec0);
max(z(:))
%}

%% Calculate the svd.  absorptionVec = U*S*V' 
[~, S, V] = svd(absorptionVec0,'econ');
vcNewGraphWin; plot(diag(S),'o-');
title('Singular values');

%{
% Confirm that we did the multiplication properly.  
% The absorptions are integers.  So, if the
% largest difference is much less than 1, we are probably OK.
tmp = U*S*V' ;
max(absorptionVec0(:) - tmp(:))
%}

%% In this formulation, the weights are U*S and the PCs are the rows of V

% Convert the PCs to images (matrices) if we did not subtract the mean
allPC = XW2RGBFormat(V,cm.rows,cm.cols);

% Convert the PCs to images if we do subtract the mean
% allPC = XW2RGBFormat(absorptionVec0',cm.rows,cm.cols);

%{
vcNewGraphWin;
% It seems to me that the first principal component accounts for the fact
% that the S-cones differ from the others
% The other terms also account for the stimulus position change.
for ii=1:3
    imagesc(allPC(:,:,ii)); axis image; colormap(gray)
    title(sprintf('PC %d',ii));
    pause(1);
end
%}

%% Find the PC weights from the stimulus absorptions

nPC = 3;
thesePC = allPC(:,:,1:nPC);
% thesePC:  <X, nPC>
thesePC = RGB2XWFormat(thesePC);
% size(thesePC)

wgts = absorptionVec*thesePC;
%{
vcNewGraphWin;
plot3(wgts(:,1),wgts(:,2),wgts(:,3),'o');
axis equal
grid on;

%}
%{
% Reconstruct

%}



end