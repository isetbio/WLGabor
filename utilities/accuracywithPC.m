function [correctness, ois, emPathStim, emPathNoStim,... 
                absorptionsStim, absorptionsNoStim, ...
                wgtsAbsorptionStim, wgtsAbsorptionNoStim] = ...
                accuracywithPC(Contrast, Freq, fov)
% Calculate the probability correct for a series of contrast and
% frequencies. Using the principal component method and SVM
%
%   Input:
%        - Contrast
%        - Frequency      
%
% Zheng Liu 

%%  Allocate space
correctness          = zeros(numel(Contrast), numel(Freq));
ois                  = cell(numel(Contrast), numel(Freq));
emPathStim           = cell(numel(Contrast), numel(Freq));
emPathNoStim         = cell(numel(Contrast), numel(Freq));
absorptionsStim      = cell(numel(Contrast), numel(Freq));
absorptionsNoStim    = cell(numel(Contrast), numel(Freq));
wgtsAbsorptionStim   = cell(numel(Contrast), numel(Freq));
wgtsAbsorptionNoStim = cell(numel(Contrast), numel(Freq));

%% Initialize basic parameters

nTimeSteps = 20;
integrationTime = 0.005;
sampleTimes = ((1:nTimeSteps) - 1) * integrationTime;   % Five ms integration time
nTrials    = 100;

% Make the harmonic with some contrast optical image
clear hparams
hparams(2)           = harmonicP;
hparams(2).freq      = Freq(1);     % Cycles per field of view
hparams(2).GaborFlag = 0;
hparams(2).contrast  = Contrast(1);

% Make the constant part
hparams(1)           = hparams(2);
hparams(1).contrast  = 0;

% Make side parame6ters
sparams.fov          = fov;

% Temporal weights for the stimulus
stimWeights = ones(1, nTimeSteps);

% tSD = 30;
% stimWeights = ieScale(fspecial('gaussian',[1,nTimeSteps],tSD),0,1);

% Summarize
fprintf('Cycles per degree %.1f\n',hparams(1).freq/sparams.fov);

%% Loop for each contrast and frequency
for c = 1 : numel(Contrast)
    hparams(2).contrast  = Contrast(c);
    for f = 1 : numel(Freq)
        hparams(2).freq      = Freq(f);     % Cycles per field of view
        fprintf('Current Constrast: %.2f; current frequency: %.1f\n', Contrast(c), Freq(f));

        % Build the sequence
        ois{c, f} = oisCreate('harmonic','blend',stimWeights, ...
            'testParameters',hparams,'sceneParameters',sparams, ...
            'sampleTimes',sampleTimes);
        % oisTmp = ois{c, f}
        % oisTmp.visualize('movie illuminance');
        
        
        %% Generate absorption for single trail WITH stimulus (multiple trails to be implemented)
        
        % trials, row, col, frames
        [absorptionsStim{c, f}, cmStim, emPath] = ccAbsorptions(ois{c, f}, nTrials);
        if sum(abs(emPath(:))) > 0, disp('Eye movements present'); end
        % whichTrial = 50;
        % thisTrial = mean(squeeze(absorptionsStim{c,f}(whichTrial,:,:,:)),3);
        % colSum = mean(thisTrial);
        % vcNewGraphWin; plot(colSum)
        % vcNewGraphWin; plot(abs(fft(colSum)));
        % emPathStim{c, f} = emPath;
        
        % Generate absorption with zero contrast for the harmonic
        hparams(2).contrast = 0;
        oisNoStim = oisCreate('harmonic','blend',stimWeights, ...
            'testParameters',hparams,'sceneParameters',sparams,...
            'sampleTimes',sampleTimes);
        [absorptionsNoStim{c, f}, cmNoStim, emPath] = ccAbsorptions(oisNoStim, nTrials, emPath);
        % whichTrial = 50;
        % thisTrial = mean(squeeze(absorptionsNoStim{c,f}(whichTrial,:,:,:)),3);
        % colSum = mean(thisTrial);
        % vcNewGraphWin; plot(colSum);
        % vcNewGraphWin; plot(abs(fft(colSum)));
        % emPathNoStim{c, f} = emPath;
        % Check absorptions w/ and w/o stimulus
        %{
            nFrame = 20;
            FramStim = reshape(squeeze(absorptionsStim(:,:,:,nFrame)),cmStim.rows,cmStim.cols);
            vcNewGraphWin; imagesc(FramStim); colormap(gray); axis image;
            FramNoStim = reshape(squeeze(absorptionsNoStim(:,:,:,nFrame)),cmStim.rows,cmStim.cols);
            vcNewGraphWin; imagesc(FramNoStim); colormap(gray); axis image;
        %}
            
            
        %% Calculate weights for PC
        wgtsStim   = wgtsGenerate(absorptionsStim{c, f}, cmStim);
        wgtsNoStim = wgtsGenerate(absorptionsNoStim{c, f}, cmNoStim);
        
            
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