% Try the svm classifier on the PC components.

%% Parameter initialization
Contrast             = 0.05; %logspace(-2, 0, 5);
Freq                 = 6; % logspace(0, 1.5, 1);
correctness          = zeros(numel(Contrast), numel(Freq));

%% Calculate the wgts
wgts = csfWgts(Contrast, Freq);

%% Process SVM
stmlType = {'Yes', 'No'};
classStmls = cell(size(wgts{1, 1}, 2),1);

for i = 1 : size(wgts{1, 1}, 2) / 2
    classStmls{i} = stmlType{1};
    classStmls{i + size(wgts{1, 1}, 2) / 2} = stmlType{2};
end

for c = 1 : size(Contrast, 2)
    for f = 1 : size(Freq, 2)
        dataStmls = wgts{c, f}';
        meanCorrect = svmProcess(dataStmls, classStmls);
        correctness(c, f) = meanCorrect;
    end
end