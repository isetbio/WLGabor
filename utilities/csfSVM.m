% Explore the svm classifier on the PC components for harmonic
% stimuli.
% 
% Description:
% First we calculate the principal components with the set of absorptions 
% that have both high contrast lavel and zero contrast level. Then 
% we change the spatial frequency and the contrast level of the stimulus.
% By using the principal components we get the weights for each 
% (freq, contrast) pair. The two biggest weights are used to represent
% the absorptions. Finally we applied SVM on the weights from each trails
% for every (freq, contrast) pair.
% 
%
% ZL/BW, 2018
%%
ieInit;
%% Parameter initialization
sFreq         = 4; 
nPCs          = 2;
fov           = 1;
sContrast     = 1;

scanFreq      = 4 %logspace(0, 1.5, 5);
scanContrast  = 1 %logspace(-3.5, 0, 5);

accuracy = zeros(numel(scanFreq), numel(scanContrast));
%% Set up the stimulus parameters
clear hparams
hparams(2)           = harmonicP;
hparams(2).freq      = sFreq;  % Set the Frequency
hparams(2).contrast  = sContrast;
hparams(1)           = hparams(2);
hparams(1).contrast  = 0;

sparams.fov = 1;

nTimeSteps = 20;
stimWeights = ones(1, nTimeSteps);
%% Set up cone mosaic parameters

integrationTime = 0.005;
sampleTimes = ((1:nTimeSteps) - 1) * integrationTime;   % Five ms integration time
nTrials    = 100;

% EM path is set to be ZERO (meaning no eyemovement for now). Will 
% implement the fixEM in the future.
empath = zeros(nTrials, nTimeSteps, 2);
%% Loop over each frequency and contrast level
for f = 1 : numel(scanFreq)
    sprintf('Current frequency is %.2f', scanFreq(f))
    
    %% Change the frequency for the stimulus
    hparams(2).freq      = scanFreq(f);  % Set the Frequency
    
    %% Create the OIS for the PC calculation

    % Stimulus for the PC
    ois = oisCreate('harmonic', 'blend', stimWeights, ...
        'testParameters', hparams, 'sceneParameters', sparams);
    % ois.visualize('movie illuminance');
    
    %% Set the coneMosaic parameters according to the OI
    cm = coneMosaic;
    cm.integrationTime = ois.timeStep;

    % Make the cm smaller than the oi size, but never smaller than 0.2 deg
    fovDegs = max(oiGet(ois.oiFixed,'fov') - 0.2, 0.2);  % Degrees
    cm.setSizeToFOV(fovDegs);

    
    %% Calculate the absorption template for the high contrast example of the stimulus

    cm.noiseFlag = 'none';
    templateHighContrast     = mean(squeeze(cm.compute(ois)), 3);
    % vcNewGraphWin; imagesc(templateHighContrast); colormap(gray);
    
    
    %% Calculate the absorption template for zero contrast of the stimulus
    hparams(2).contrast  = 0.0;
    ois = oisCreate('harmonic', 'blend', stimWeights, ...
        'testParameters', hparams, 'sceneParameters', sparams);
    templateZeroContrast     = mean(squeeze(cm.compute(ois)), 3);

    %% Calculate the PCs (the whole sets of principal components)
    sprintf('Calculating the whole sets of PCs')
    
    PCs = csfPC(templateHighContrast, templateZeroContrast, nTrials);
    
    %% Create a "blank" pattern without stimulus
    hparams(2).contrast  = 0.0;
    ois = oisCreate('harmonic', 'blend', stimWeights, ...
        'testParameters', hparams, 'sceneParameters', sparams);
    % ois.visualize('movie illuminance');

    % No eyemovements for now.
    cm.noiseFlag = 'random';
    absorptionsNoise = cm.compute(ois, 'empath', empath);
    meanabsorptionsNoise = mean(absorptionsNoise, 4);

    %% Calculate the wgts for the Noise
    sprintf('Calculating the Noise PC weights')
    
    wgtsNoiseWhole = csfWgts(meanabsorptionsNoise, PCs);
    wgtsNoise      = wgtsNoiseWhole(1:nPCs, :);

    %% The second for loop to scan contrast level
    
    for c = 1 : numel(scanContrast)
        sprintf('Current contrast level is %.2f', scanContrast(c))
        %% Create the test stimulus at a lower contrast level
        
        hparams(2).contrast  = scanContrast(c);
        ois = oisCreate('harmonic', 'blend', stimWeights, ...
            'testParameters', hparams, 'sceneParameters', sparams);
        % ois.visualize('movie illuminance');

        %Eye movement is set to zero
        absorptions = cm.compute(ois, 'empath', empath);
        meanabsorptionsStimulus = mean(absorptions, 4);
        %{
            thisTrial = 10;
            thisFrame = 5;
            vcNewGraphWin; imagesc(squeeze(absorptions(thisTrial, :, :, thisFrame))); colormap(gray);
        %}
        
        %% Calculate the wgts for the absorptions
        sprintf('Calculating the stimulus PC weights')
        
        wgtsStimulusWhole = csfWgts(meanabsorptionsStimulus, PCs);
        wgtsStimulus      = wgtsStimulusWhole(1:nPCs, :);
        %{
            figure;
            hold all;
            plot(abs(wgtsStimulusWhole(:,1)),'-o'); plot(abs(wgtsNoiseWhole(:,1)), '-o');
            legend('Stimulus', 'NoStimulus')
        %}
        %% Process SVM

        % correctness   = zeros(numel(Contrast), numel(Freq)); will be implemented
        % in the future.
        sprintf('Processing the SVM')
        classStmls = [ones(size(wgtsStimulus, 2), 1); zeros(size(wgtsNoise, 2), 1)];

        dataStmls = [wgtsStimulus'; wgtsNoise'];
        meanCorrect = svmProcess(dataStmls, classStmls);
        accuracy(f, c) = meanCorrect;
    end
end

%% Plot the results

figure;
hold all;

pp = cell(1,numel(scanFreq));

for i = 1 : numel(scanFreq)
	plot(log10(scanContrast), accuracy(i,:),'-o')
	pp{i} = sprintf('Cycles per degree = %.2f',scanFreq(i));
end
legend(pp);
xlabel('Log Contrast')
ylabel('Mean Correctness')
