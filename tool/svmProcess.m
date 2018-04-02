function meanCorrect = svmProcess(dataStmls, classStmls)
    svm = fitcsvm(dataStmls, classStmls);
    kFold = 10;
    CVSVMOptimize = crossval(svm,'KFold',kFold);
    probabilityCorrect = 1 - kfoldLoss(CVSVMOptimize,'lossfun','classiferror','mode','individual');

    fprintf('Mean probability correct %.2f\n',mean(probabilityCorrect));
    meanCorrect = mean(probabilityCorrect);
end