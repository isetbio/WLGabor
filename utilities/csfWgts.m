function wgts = csfWgts(freq, varargin)
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

% Removes spaces and forces to lower case
varargin = ieParamFormat(varargin);

p.addRequired('freq',@isscalar);
p.addParameter('npcs',2,@isscalar);
p.addParameter('contrast',1.0,@isscalar);

p.parse(freq,varargin{:});
contrast = p.Results.contrast;
nPCs     = p.Results.npcs;

%% Initiation parameters

wgts = cell(numel(contrast), numel(freq));

%% Calculate the wgts for another (more in the future) (C, F) with the PC above.
for c = 1 : size(contrast, 2)
    for f = 1 : size(freq, 2)
        thisFreq = freq(f);
        thisContrast = contrast(c);
        
        [~, samples] = csfPC(thisFreq, thisContrast, fov, nPCs);
        wgts{c, f} = PC' * samples;
        wgtsnPC = wgts{c, f}(1:nPCs, :);
        wgts{c, f} = wgtsnPC;
        
        
        %{
            rows = size(samples, 1) ^ 0.5;
            cols = rows;
            sample_rec = PC(:, 1:nPCs) * wgtsnPC;
            nframe = 50;
            vcNewGraphWin;imagesc(reshape(sample_rec(:,nframe),rows, cols))
            vcNewGraphWin;imagesc(reshape(samples(:,nframe),rows, cols))
        %}
    end
end




