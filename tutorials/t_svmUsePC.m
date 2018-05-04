function t_svmUsePC
%% Here is a curve we want

%{
   %For a series of contrast levels, say from 1 percent to 50%
   % maybe 5 or 7 different levels  
   sContrast = logspace(-2,-0.3,5)
   sFreq     = logspace(0,1.5, 8)
   calculate the probability correct for 3 principal component approx
   make the image or matrix that has rows of contrast, columns of
   spatial frequency, entries of probability correct detection
   Show it as an image or a mesh or a set of curves

   Make sure you do some checks that when it is 0 contrast, you are at
   chance.

   Make sure that the Gabor stimulus with contrast and without have
   the same eye movements!

   Have a great time.
%}

%% Here is the curve that I created

sContrast   = 0.02;%logspace(-2.5, -1e-5, 4);
sFreq       = 4; %logspace(0, 1.5, 4);
fov         = 1;
[probCorrect, ois, emPathStim, emPathNoStim,...
  absorptionsStim, absorptionsNoStim, wgtsAbsorptionStim, wgtsAbsorptionNoStim]= accuracywithPC(sContrast, sFreq, fov);

%% Sanity check for emPath
%{
    c = 1; 
    f = 1;
    nTrials = size(emPathStim{c, f}, 1);
    for i = 1 : nTrials
        subplot(ceil(nTrials/3), 3, i);
        xPosStim = emPathStim{c,f}(i,:,1);
        yPosStim = emPathStim{c,f}(i,:,2);
        plot(xPosStim, yPosStim,'-ko');
        hold all;
        xPosNoStim = emPathNoStim{c,f}(i,:,1);
        yPosNoStim = emPathNoStim{c,f}(i,:,2);
        plot(xPosNoStim, yPosNoStim,'-go');

        plot(xPosStim - xPosNoStim, yPosStim - yPosNoStim, '-bo');
    end
%}

%% Sanity check for original absorption and SVD reconstruction
%{
    c = 1;
    f = 1;
    nTrials = size(absorptionsStim{c, f}, 1);
    for i = 1 : 15
        subplot(ceil(15/3), 3, i);
        meanabsorption = mean(absorptionsStim{c, f}, 4);
        imagesc(squeeze(meanabsorption(i,:,:)));colormap(gray); axis image;
    end
%}
%{
    Check the correctness of weight calculation with 
    c = 1;
    f = 1;
    [Trials, rows, cols, Frames] = size(absorptionsStim{c, f});
    nFrame = 100;
    nTrial = 55;
    wgtsTemp = wgtsAbsorptionStim{c, f};
    absorptionsStimTemp = absorptionsStim{c, f};
    approx = reshape(wgtsTemp((nTrial - 1) * 100 + nFrame,:,:), rows, cols);
    vcNewGraphWin; imagesc(approx); colormap(gray); axis image;
    actual = reshape(absorptionsStimTemp(nTrial,:, :, nFrame), rows, cols);
    vcNewGraphWin; imagesc(approx); colormap(gray); axis image;
%}

%% plot Figure
figure;
stem3(sFreq, sContrast ,probCorrect,'linestyle','none')
figure;
hold all;

pp = cell(1,numel(sFreq));

for i = 1 : numel(sFreq)
    plot(log10(sContrast), probCorrect(:,i),'-o')
    pp{i} = sprintf('Cycles per degree = %.2f',sFreq(i));
end
legend(pp);
xlabel('Log Contrast')
ylabel('Mean Correctness')
end