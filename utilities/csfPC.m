function PCs = csfPC(templateHC, templateZC, nTrials)
% Calculate the principal components from signal and noise trials
% 
%  Description:
%  We pre-compute the likely stimuli,  The logic is to use a high
%  contrast stimulus and make many samples of the mean absoprtion
%  pattern.  Then use a noise stimulus and make many samples.  Then
%  put all the means (noise and signal) into a big matrix.  Compute
%  the SVD to get the top N principal components.  When there are only
%  two classes and no eye movements we probably only need two
%  components. If we have eye movements, then we will probably need a
%  few more principal components.
%
%  We use these components to do CSF calculation.
%  
%  After checking the PCs with different contrast level, we decide to
%  choose to use the PCs with the contrast = 1 (for now), and calculate the
%  weight 
%
%  
%  
%
%% ZL/BW

%{

sFreq = 4;
fov  = 1;
nPCs = 20;

% sFreq, sContrast, fov, nPCs

PC1 = csfPC(sFreq,1.0,fov,nPCs);
vcNewGraphWin; imagesc(reshape(PC1(:,1),cmStim.rows,cmStim.cols))
PC2 = csfPC(7,0.01,fov,nPCs);
vcNewGraphWin; imagesc(reshape(PC2(:,1),cmStim.rows,cmStim.cols))

vcNewGraphWin; plot(PC1(:,1),PC2(:,1),'.'); 
grid on; axis equal; identityLine;

PC1(:,1)'*PC2(:,1)

%}
%% Parse the inputs

p = inputParser;

p.addRequired('templateHC',@isnumeric);
p.addRequired('templateZC',@isnumeric);
p.addRequired('nTrials', @isnumeric);
p.parse(templateHC, templateZC, nTrials);

templateHC = p.Results.templateHC;
templateZC = p.Results.templateZC;
nTrials    = p.Results.nTrials;

%% Vectorize the template

templateHCRep = repmat(templateHC, [1, 1, nTrials / 2]);
templateZCRep = repmat(templateZC, [1, 1, nTrials / 2]);
templateHCVec = RGB2XWFormat(templateHCRep);
templateZCVec = RGB2XWFormat(templateZCRep);
templateRep = [templateHCVec templateZCVec];
%% Use the template to calculate the PCs

[PCs, ~, ~] = svd(templateRep,'econ'); 
%{
    wgts = PCs' * templateRep;
    sample = PCs(:,1:2) * wgts(1:2,:);
    anycopy = 30;
    [row, col] = size(templateHC);
    vcNewGraphWin; imagesc(reshape(sample(:,anycopy), row, col));
%}
end