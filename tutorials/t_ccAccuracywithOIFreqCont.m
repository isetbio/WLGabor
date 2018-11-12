function meanProbCorrect = t_ccAccuracywithOIFreqCont
%%
% Calculate the accuracy of SVM discrimination on whether there is stimulus or
% not given certain contrast and frequency range.
%
%
% Output:
%       meanProbCorrect     - Mean probability of correctness
%%
ieInit;
%%  call t_ccDiscriminate
% Set the certain range of contrast, frequency and FOV for the OI sequence, calculate the mean
% probability of correctness by call t_ccDiscriminate
    contrast = [0.1 : 0.3 : 0.7];
    frequency = [1 : 1 : 15];
    meanProbCorrect = zeros(numel(frequency), numel(contrast));
    fov = 0.6;
    
    for f = 1 : numel(frequency)
        for c = 1 : numel(contrast)
            meanProbCorrect(f, c) = t_ccDiscriminate(frequency(f), ...
                                                        contrast(c), fov);
        end
        
    end
    
    
    
%%  Plot the result
% This function is used to plot the mean probability of correctness with
% different frequency and contrast.
% Todo: fit the curve and find the freqeuency under different contrast that
% has 80% (75%) correctness.
    plotMeanProb(meanProbCorrect, frequency, contrast);
   
end