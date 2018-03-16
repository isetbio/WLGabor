function meanProbCorrect = t_ccAccuracywithOIFreqCont
%%
% Calculate the accuracy of SVM discrimination on whether there is stimulus or
% not given certain contrast and frequency range.
%
%
% Output:
%       meanProbCorrect     - Mean probability of correctness
%%  call t_ccDiscriminate
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
    plotMeanProb(meanProbCorrect, frequency, contrast);
   
end