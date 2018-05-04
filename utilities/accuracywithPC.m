function [correctness, ois, emPathStim, emPathNoStim,... 
                absorptionsStim, absorptionsNoStim, wgtsAbsorptionStim, wgtsAbsorptionNoStim] = accuracywithPC(Contrast, Freq, fov)
%
% Calculate the correctness with a series of contrast and frequency with 
% the principal component by using SVM

%   Input:
%        - Contrast
%        - Frequency
        

%%
correctness  = zeros(numel(Contrast), numel(Freq));
ois          = cell(numel(Contrast), numel(Freq));
emPathStim   = cell(numel(Contrast), numel(Freq));
emPathNoStim = cell(numel(Contrast), numel(Freq));
absorptionsStim = cell(numel(Contrast), numel(Freq));
absorptionsNoStim = cell(numel(Contrast), numel(Freq));
wgtsAbsorptionStim = cell(numel(Contrast), numel(Freq));
wgtsAbsorptionNoStim = cell(numel(Contrast), numel(Freq));
%%
for c = 1 : numel(Contrast)
    for f = 1 : numel(Freq) 
    %%  Generate oisequence
    clear hparams

    % Make the time varying part
    hparams(2) = harmonicP;
    hparams(2).freq      = Freq(f);     % Cycles per field of view
    hparams(2).GaborFlag = 0;
    hparams(2).contrast  = Contrast(c);

    % Make the constant part
    hparams(1) = hparams(2);
    hparams(1).contrast = 0;
    sparams.fov = fov;
    fprintf('Current Constrast: %.2f; current frequency: %.1f\n', Contrast(c), Freq(f));
    fprintf('Cycles per degree %.1f\n',hparams(1).freq/sparams.fov);

    % These are the scalar over time for the oi sequence
    nTimeSteps = 100;
    tSD = 30;
    %stimWeights = ieScale(fspecial('gaussian',[1,nTimeSteps],tSD),0,1);
    stimWeights = ones(1, nTimeSteps);
    % Build the sequence
    ois{c, f} = oisCreate('harmonic','blend',stimWeights, ...
        'testParameters',hparams,'sceneParameters',sparams);
    %{
        oisTmp = ois{c, f}
        oisTmp.visualize('movie illuminance');
    %}
    %% Generate absorption for single trail WITH stimulus (multiple trails to be implemented)

    nTrials = 100;
    [absorptionsStim{c, f}, cmStim, emPath] = ccAbsorptions(ois{c, f}, nTrials);
    emPathStim{c, f} = emPath;
    %% Generate absorption without stimulus
    hparams(2).contrast = 0;
    oisNoStim = oisCreate('harmonic','blend',stimWeights, ...
        'testParameters',hparams,'sceneParameters',sparams);
    [absorptionsNoStim{c, f}, cmNoStim, emPath] = ccAbsorptions(oisNoStim, nTrials, emPath);
    emPathNoStim{c, f} = emPath;
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
    [wgtsStim, wgtsAbsorptionStim{c, f}]      = wgtsGenerate(absorptionsStim{c, f}, cmStim, nTrials);
    [wgtsNoStim, wgtsAbsorptionNoStim{c, f}]  = wgtsGenerate(absorptionsNoStim{c, f}, cmNoStim, nTrials);
    

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