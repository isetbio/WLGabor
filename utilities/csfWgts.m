% Calculate the PC weights for each (sContrast, sFreq). Used to 
% confirm the PC calculate from (sContrast = 1, sFreq = 4) can be used
% for all (sContrast, sFreq) pair
% 

%% Initiation parameters
nPCs = 200;
%% Calculate the PC for sContrast = 1, sFreq = 4
sFreq = 4;
sContrast = 1;
fov  = 1;
PC = csfPC(sFreq, sContrast, fov, nPCs);
%% Calculate the wgts for another (more in the future) (C, F) with the PC above.
sFreq = 4;
sContrast = 1;
fov  = 1;

[~, samples] = csfPC(sFreq, sContrast, fov, nPCs);

wgts = PC' * samples;
wgtsnPC = wgts(1:nPCs, :);
%{
    rows = size(samples, 1) ^ 0.5;
    cols = rows;
    sample_rec = PC(:, 1:nPCs) * wgtsnPC;
    nframe = 150;
    vcNewGraphWin;imagesc(reshape(sample_rec(:,nframe),rows, cols))
    vcNewGraphWin;imagesc(reshape(samples(:,nframe),rows, cols))
%}
