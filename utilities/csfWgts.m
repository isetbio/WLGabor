function wgts = csfWgts(absorptions, PCs, nPC)
% Calculate the PC weights for each (sContrast, sFreq). Used to 
% confirm the PC calculate from (sContrast = 1, sFreq = 4) can be used
% for all (sContrast, sFreq) pair
%
%  
%
% Examples:
%    freq = 4;
%    wgts = csfWgts(freq,'contrast',0.5,'n pcs',3);
%

%% Parse the inputs
p = inputParser;

p.addRequired('absorptions',@isnumeric);
p.addRequired('PCs',@isnumeric);
p.addRequired('nPC',@isnumeric);

p.parse(absorptions, PCs, nPC);
absorptions = p.Results.absorptions;
PCs = p.Results.PCs;
nPC     = p.Results.nPC;
%% Initialte parameters
nTrials = size(absorptions, 1);
vecLength = size(PCs, 1);
%% Vecterize the absorptions
absorptions = permute(absorptions, [2 3 1]);
absorptionsVec = RGB2XWFormat(absorptions);

%% Calculate the wgts
wgtsWhole = PCs' * absorptionsVec;
wgts = wgtsWhole(1:nPC, :);

%{
    % Check if the PCs and thisPC can be considered as the "Same"
    [thisPC, ~, ~] = svd(absorptionsVec, 'econ');
    PCs(:,1)' * thisPC(:,1)
%}
%{
    sample =  PCs(:, 1:nPC) * wgtsWhole(1:nPC, :);
    thisTrial = 1;
    row = vecLength^0.5;
    col = row;
    vcNewGraphWin; imagesc(squeeze(reshape(sample(:,thisTrial), row,...
    col)));
    vcNewGraphWin; imagesc(squeeze(reshape(absorptionsVec(:,thisTrial), row,...
    col)));
%}




end
