function meanStmlsPlot(frame)
        [trials, row, col] = size(frame);
        figure;
        nPlot = 10;
        for i = 1 : trials
            subplot(10, 10, i);
            imshow(squeeze(frame(i,:,:) / max(max(frame(i,:,:)))));
        end
end