function plotMeanProb(meanProbCorrect, frequency, contrast)
%%  plot the Probability of correctness for fixed contrast given frequency  
    figure;
    for c = 1 : numel(contrast)
        hold on;
        plot(frequency, meanProbCorrect(:, c),'o-');
    end
    xlabel('Spatial frequency')
    ylabel('Probability of correctness')
end