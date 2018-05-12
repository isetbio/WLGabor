function [PC, samples] = csfPC(sFreq, fov)
% Calculate the principal components from signal and noise trials
%
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
%  When there are no eye movements, 
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


%% Put the mean absorptions from all trials into one big matrix and calculate the PCs

samples = [RGB2XWFormat(meanStim), RGB2XWFormat(meanNoise)];

[U, ~, ~] = svd(samples,'econ');
    %{
      nPCs = 50;
      [A, B, C] = svd(samples, 'econ');
      wgts = A' * samples;
      wgts_two = wgts(1:nPCs,:);
      
    %}
PC = U;

    %{
        sample_rec = PC(:, 1:nPCs) * wgts_two;
        nframe = 120;
        vcNewGraphWin;imagesc(reshape(sample_rec(:,nframe),cmStim.rows,cmStim.cols))
        vcNewGraphWin;imagesc(reshape(samples(:,nframe),cmStim.rows,cmStim.cols))
    %}
% sValues = diag(S);
% vcNewGraphWin; semilogy(sValues(1:4));
% vcNewGraphWin; imagesc(reshape(PC(:,1),cmStim.rows,cmStim.cols))
% vcNewGraphWin; imagesc(reshape(PC(:,2),cmStim.rows,cmStim.cols))

end