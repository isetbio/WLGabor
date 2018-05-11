function wgts = csfWgts(contrast, freq)
% Calculate the PC weights for each (sContrast, sFreq). Used to 
% confirm the PC calculate from (sContrast = 1, sFreq = 4) can be used
% for all (sContrast, sFreq) pair
% 

%% Initiation parameters
nPCs = 2;

wgts = cell(size(contrast,2), size(freq, 2));
%% Calculate the PC for sContrast = 1, sFreq = 4
sFreq = 4;
sContrast = 1;
fov  = 1;
PC = csfPC(sFreq, sContrast, fov, nPCs);
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




