function wgts = wgtsGenerate(absorptions, cm, varargin)
% Calculate the weights of the spatial principal components of the
% absorptions
%
% Description
%  For each trial, we calculate the principal components of the
%  spatial cone absorption pattern, with the PCs derived from the
%  muultiple time steps. We then calculate the PC weights at each time
%  step, and return the mean of these weights.
%  
% Inputs
%  absorptions - 4D vector of cone [trial,row,col,time] absorptions
%  cm - cone mosaic
%
% Optional key/val parameters
%   N/A
%
% Returns:
%   wgts - Returned temporal series of weights for the absorptions
%
% Zheng Liu

%%
p = inputParser;
p.addRequired('absorptions',@isnumeric);
p.addRequired('cm',@(x)(isa(x,'coneMosaic')));

p.addParameter('npc',3,@isinteger);

p.parse(absorptions,cm,varargin{:});

nPC     = p.Results.npc;
nTrials = size(absorptions,1);

% wgtsAbsorption = [];
wgts = zeros(nTrials,nPC);
% wgtsAbsorption = zeros(nTrials,nPC);

for ii = 1 : nTrials
    % Vectorize the absorption for this trial.  It becomes (r,c,time)
    absorptionVec = permute(squeeze(absorptions(ii,:,:,:)),[1 2 3]);
    
    % Each row is the spatial array of cone responses at one sample
    % time.
    absorptionVec = RGB2XWFormat(absorptionVec)';
    
    %% Calculate the svd.  absorptionVec = U*S*V'
    [~, S, V] = svd(absorptionVec,'econ');
    % vcNewGraphWin; plot(diag(S),'o-');
    % title('Singular values');
    
    %% In this formulation, the weights are U*S and the PCs are the rows of V
    
    % Convert the PCs to images (matrices) if we did not subtract the mean
    allPC = XW2RGBFormat(V,cm.rows,cm.cols);
    
    %% Find the PC weights from the stimulus absorptions
    
    thesePC = allPC(:,:,1:nPC);
    % thesePC:  <X, nPC>
    thesePC = RGB2XWFormat(thesePC);
    % size(thesePC)
    
    % wgts has the size of [steps * nPC]
    wgtsTemp           = absorptionVec * thesePC;
    wgts(ii,:) = mean(wgtsTemp,1);
    
    %{
          Check the correctness of weight calculation with

          absorptionTemp = wgtsTemp * thesePC';
          nFrame = 30;
          approx = reshape(absorptionTemp(nFrame,:),cm.rows,cm.cols);
          vcNewGraphWin; imagesc(approx); colormap(gray); axis image;
          actual = reshape(absorptionVec(nFrame,:),cm.rows,cm.cols);
          vcNewGraphWin; imagesc(approx); colormap(gray); axis image;
    %}
end
end