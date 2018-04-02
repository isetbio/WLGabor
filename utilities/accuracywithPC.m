function correctness = accuracywithPC(Contrast, Freq, fov)
%
% Calculate the correctness with a series of contrast and frequency with 
% the principal component by using SVM

%   Input:
%        - Contrast
%        - Frequency
        

%%
correctness = zeros(numel(Contrast), numel(Freq));

%%
for c = 1 : numel(Contrast)
    for f = 1 : numel(Freq)
    %%  Generate oisequence
    clear hparams

    % Make the time varying part
    hparams(2) = harmonicP;
    hparams(2).freq      = Freq(f);     % Cycles per field of view
    hparams(2).GaborFlag = 0.2;
    hparams(2).contrast  = Contrast(c);

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
    
    %% Generate absorption for single trail WITH stimulus (multiple trails to be implemented)

    nTrials = 1;
    [absorptionsStim, cmStim, emPath] = ccAbsorptions(ois, nTrials);
    %% Generate absorption without stimulus
    hparams(2).contrast = 0;
    ois = oisCreate('harmonic','blend',stimWeights, ...
        'testParameters',hparams,'sceneParameters',sparams);
    [absorptionsNoStim, cmNoStim] = ccAbsorptions(ois, nTrials, emPath);
%%
    %{
        Check absorptions w/ and w/o stimulus
        nFrame = 20;
        FramStim = reshape(squeeze(absorptionsStim(:,:,:,nFrame)),cmStim.rows,cmStim.cols);
        vcNewGraphWin; imagesc(FramStim); colormap(gray); axis image;
        FramNoStim = reshape(squeeze(absorptionsNoStim(:,:,:,nFrame)),cmStim.rows,cmStim.cols);
        vcNewGraphWin; imagesc(FramNoStim); colormap(gray); axis image;
    %}
    
    
    %% Calculate weights for PC
    wgtsStim   = wgtsGenerate(absorptionsStim, cmStim);
    wgtsNoStim = wgtsGenerate(absorptionsNoStim, cmNoStim);
    

    %% Combine the stimuli and the labels
    stmlType = {'Yes', 'No'};
    
    dataStmls = [wgtsStim; wgtsNoStim];
    % size(dataStmls)
    classStmls = cell(2 * size(wgtsStim, 1),1);
    for i = 1 : size(wgtsStim, 1)
        classStmls{i} = stmlType{1};
        classStmls{i + size(wgtsStim, 1)} = stmlType{2};
    end
    
    meanCorrect = svmProcess(dataStmls, classStmls);
    correctness(c, f) = meanCorrect;
    end
end

end