function wgts = wgtsGenerate(absorptions, cm, nTrials)
    wgts = [];
    for i = 1 : nTrials
        %% Vectorize the absorption
        absorptionVec = permute(squeeze(absorptions(i,:,:,:)),[1 2 3]);
        absorptionVec = RGB2XWFormat(absorptionVec)';

        %% Calculate the svd.  absorptionVec = U*S*V' 

        [~, S, V] = svd(absorptionVec,'econ');
        %vcNewGraphWin; plot(diag(S),'o-');
        %title('Singular values');

        %% In this formulation, the weights are U*S and the PCs are the rows of V


        % Convert the PCs to images (matrices) if we did not subtract the mean
        allPC = XW2RGBFormat(V,cm.rows,cm.cols);

        %% Find the PC weights from the stimulus absorptions

        nPC = 3;
        thesePC = allPC(:,:,1:nPC);
        % thesePC:  <X, nPC>
        thesePC = RGB2XWFormat(thesePC);
        % size(thesePC)

        % wgts has the size of [steps * nPC]
        wgtsTemp = absorptionVec*thesePC;
        wgtsMean = mean(wgtsTemp,1);
        wgts = [wgts;wgtsMean];
        %{
          Check the correctness of weight calculation with 

          absorptionTemp = wgts * thesePC';
          nFrame = 30;
          approx = reshape(absorptionTemp(nFrame,:),cm.rows,cm.cols);
          vcNewGraphWin; imagesc(approx); colormap(gray); axis image;
          actual = reshape(absorptionVec(nFrame,:),cm.rows,cm.cols);
          vcNewGraphWin; imagesc(approx); colormap(gray); axis image;
        %}
    end
end