function wgtsWhole = csfWgts(absorptions, PCs)
% Calculate the PC weights for each (sContrast, sFreq). Used to 
% confirm the PC calculate from (sContrast = 1, sFreq = 4) can be used
% for all (sContrast, sFreq) pair
% 
% Description: 
%   The weights are calculated simply by multiplying the transpose of the
%   PCs with the absorpitons matrix.
% 
%
% 
%

%% Parse the inputs
p = inputParser;

p.addRequired('absorptions',@isnumeric);
p.addRequired('PCs',@isnumeric);

p.parse(absorptions, PCs);
absorptions = p.Results.absorptions;
PCs = p.Results.PCs;

%% Vecterize the absorptions
absorptions = permute(absorptions, [2 3 1]);
absorptionsVec = RGB2XWFormat(absorptions);

%% Calculate the wgts
wgtsWhole = PCs' * absorptionsVec;
%{
    plot(abs(wgtsWhole(:,1)), '-o')
%}

%{
    % Check if the PCs and thisPC can be considered as the "Same"
    [thisPC, ~, ~] = svd(absorptionsVec, 'econ');
    PCs(:,1)' * thisPC(:,1)
%}
%{
    nPC = 2;
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
