% Generate data for the CNN
%
% Required: isetbio
%
%%
ieInit;

%% Parameter initialization
sFreq         = 4;
nPCs          = 2;
fov           = 1;
sContrast     = 1;

% Save the generated data?
saveFlag = true;
saveName = 'freq_8_contrast_00005';
% accuracy = zeros(numel(scanFreq), numel(scanContrast));

%% Set up the stimulus parameters
clear hparams
hparams(2)           = harmonicP;
hparams(2).freq      = sFreq;  % Set the Frequency
hparams(2).contrast  = sContrast;
hparams(1)           = hparams(2);
hparams(1).contrast  = 0;

sparams.fov = 1.5;

nTimeSteps = 20;
stimWeights = ones(1, nTimeSteps);

%% Set up cone mosaic parameters

integrationTime = 0.005;
sampleTimes = ((1:nTimeSteps) - 1) * integrationTime;   % Five ms integration time
nTrials    = 100;

%% Randomly generate parameters:
% scanFreq
% scanContrast
% aberrations
rng(1); % Set the random seed

% Scan for test set
scanFreq = 8; % logspace(0.1, 1.5441, 5);
scanContrast = 0.005;% logspace(-3.5, -1, 10);
numSamples = 1000;

nImages = length(scanFreq)*length(scanContrast)*numSamples;

% Aberrations
measPupilMM = 4.5; % 4.5 mm pupil size (Thibos data)
calcPupilMM = 3; % 3 mm pupil size (what we want to calculate with)
zCoeffs = VirtualEyes(nImages,measPupilMM);

%% Generate images

% testTemp = zeros(249,249,length(scanFreq), 2*nImages);
% With stimulus and noise
imgNoiseStimulus = zeros(249,249, nImages);
imgNoiseNoStimulus = zeros(249,249, nImages);

% With stimulus and with no noise
imgNoNoiseStimulus = zeros(249,249,nImages);
imgNoNoiseNoStimulus = zeros(249,249,nImages);

% Labels
imgLabelsStimulus = zeros(nImages,1);
imgLabelsNoStimulus = zeros(nImages,1);

% Contrast
imgContrasts = zeros(nImages,1);

% Frequencies
imgFreqs = zeros(nImages,1);

k = 1;
ii = 1;
for ff = 1 : length(scanFreq)
    for cc = 1:length(scanContrast)
        for nn = 1:numSamples
            
            fprintf('Generating image: %i \n',ii)
            %% Change the frequency and contrast for the stimulus
            hparams(2).freq      = scanFreq(ff);  % Set the Frequency
            hparams(2).contrast  = scanContrast(cc);
            
            %% Create the oi with aberrations
            
            %z = zeros(65,1);
%             z(1:13) = zCoeffs(ii,1:13);
             ii = ii + 1;
%                         
%             % Create the example subject
%             sbjWvf = wvfCreate;                                     % Initialize
%             sbjWvf = wvfSet(sbjWvf,'zcoeffs',z);                    % Zernike
%             sbjWvf = wvfSet(sbjWvf,'measured pupil',measPupilMM);   % Data
%             sbjWvf = wvfSet(sbjWvf,'calculated pupil',calcPupilMM); % What we calculate
%             sbjWvf = wvfSet(sbjWvf,'measured wavelength',550);
%             sbjWvf = wvfSet(sbjWvf,'calc wave',[450:10:650]');            % Must be a column vector
%             sbjWvf = wvfComputePSF(sbjWvf);
%             oi = wvf2oi(sbjWvf);
            
            % Default
            oi = oiCreate('wvf human');
            
            %% Create the OIS
            
            ois = oisCreate('harmonic', 'blend', stimWeights, ...
                'testParameters', hparams, 'sceneParameters', sparams,...
                'oi',oi);
            % ois.visualize('movie illuminance');
            
            %% Set the coneMosaic parameters according to the OI
            cm = coneMosaic;
            cm.integrationTime = ois.timeStep;
            
            % Make the cm smaller than the oi size, but never smaller than 0.2 deg
            fovDegs = max(oiGet(ois.oiFixed,'fov') - 0.2, 0.2);  % Degrees
            cm.setSizeToFOV(fovDegs);
            
            %% Calculate the absorption template for the high contrast example of the stimulus
            
            cm.noiseFlag = 'random';
            imgNoiseStimulus(:,:,k) = mean(squeeze(cm.compute(ois)), 3);
            imgLabelsStimulus(k) = 1;
            imgContrasts(k) = scanContrast(cc);
            imgFreqs(k) = scanFreq(ff);
            
            % Calculate without noise
            cm.noiseFlag = 'none';
            imgNoNoiseStimulus(:,:,k) = mean(squeeze(cm.compute(ois)), 3);
            
            
            
            %% Create a "blank" pattern without stimulus
            hparams(2).contrast  = 0.0;
            ois = oisCreate('harmonic', 'blend', stimWeights, ...
                'testParameters', hparams, 'sceneParameters', sparams,...
                'oi',oi);
            % ois.visualize('movie illuminance');
            
            cm.noiseFlag = 'random';
            imgNoiseNoStimulus(:,:,k) = mean(squeeze(cm.compute(ois)), 3);
            imgLabelsNoStimulus(k) = 0;
            imgContrasts(k) = scanContrast(cc);
            imgFreqs(k) = scanFreq(ff);
            
            % Calculate without noise
            cm.noiseFlag = 'none';
            imgNoNoiseNoStimulus(:,:,k) = mean(squeeze(cm.compute(ois)), 3);
            
            k = k+1;
            

        end
    end
end


%% Exam the mean cone absorption
meanImgNoiseStimulus = mean(imgNoiseStimulus, 3);
meanImgNoNoiseStimulus = mean(imgNoNoiseStimulus, 3);
meanImgNoiseNoStimulus = mean(imgNoiseNoStimulus, 3);
meanImgNoNoiseNoStimulus = mean(imgNoNoiseNoStimulus, 3);

figure(1); imagesc(meanImgNoiseStimulus, [6 8]); colorbar; title('meanImgNoiseStimulus');
figure(2); imagesc(meanImgNoNoiseStimulus, [6 8]); colorbar; title('meanImgNoNoiseStimulus');

figure(3); imagesc(meanImgNoiseNoStimulus, [6 8]); colorbar; title('meanImgNoiseNoStimulus');
figure(4); imagesc(meanImgNoNoiseNoStimulus, [6 8]); colorbar; title('meanImgNoNoiseNoStimulus');

figure(5); imagesc(meanImgNoiseStimulus - meanImgNoiseNoStimulus, [-5 5]); colorbar; title('mean Noise sub');
figure(6); imagesc(meanImgNoNoiseStimulus - meanImgNoNoiseNoStimulus, [-0.02 0.02]);colorbar; title('mean No Noise sub');
%% Crop
imgNoiseStimulus = imgNoiseStimulus(:,:,:);
imgNoiseNoStimulus = imgNoiseNoStimulus(:,:,:);
imgNoNoiseStimulus = imgNoNoiseStimulus(:, :,:);
imgNoNoiseNoStimulus = imgNoNoiseNoStimulus(:, :,:);
%% Save everything

if(saveFlag)
    currDate = datestr(now,'mm-dd-yy_HH_MM');
    save(sprintf('%s_%s.mat',saveName,currDate),...
        'imgNoiseStimulus',...
        'imgNoiseNoStimulus',...
        'imgNoNoiseStimulus',...
        'imgNoNoiseNoStimulus',...
        'imgLabelsStimulus',...
        'imgLabelsNoStimulus',...
        'imgContrasts',...
        'imgFreqs');
end

%% Need to be fixed Below

%% Display a random sampling of images

n = 16;
figure(1);
[dataSamp,idx] = datasample(Images,n,3);
labelSamp = imgLabels(idx);
contrastSamp = imgContrasts(idx);
freqSamp = imgFreqs(idx);
k = 1;
for ii = 1:n
    subplot(4,4,k);
    imagesc(dataSamp(:,:,ii));
    axis image; axis off;
    title(sprintf('label = %i \n c = %0.2f | f = %0.2f',...
        labelSamp(ii),contrastSamp(ii),freqSamp(ii)))
    k = k+1;
end

%% Generate range

unique_c = unique(imgContrasts);
unique_c = unique_c(2:2:end);
unique_f = unique(imgFreqs);
N = length(unique_c);
M = length(unique_f);

figure;
k = 1;
for cc = 1:length(unique_c)
    for ff = 1:length(unique_f)
        
        curr_c = unique_c(cc);
        curr_f = unique_f(ff);
        
        subplot(N,M,k)
        curr_i = ((imgContrasts == curr_c) & ...
            (imgLabels == 1) & ...
            (imgFreqs == curr_f));
        curr_i = find(curr_i == 1); % Choose the first image
        curr_i = curr_i(1);
        curr_image = Images(:,:,curr_i);
        imagesc(curr_image);
        axis image; axis off;
        title(sprintf('label = %i \n c = %0.2f | f = %0.2f',...
            imgLabels(curr_i),imgContrasts(curr_i),imgFreqs(curr_i)))
        k = k +1;
    end
end
