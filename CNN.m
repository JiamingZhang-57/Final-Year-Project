clear all;clc;
%% Data input
load('Trainnoise.mat');
load('Testset.mat')
D=[1,0,0,1,0,1,0,1,1,0];
D=D';
D=categorical(D);
%% Define Layers
layers = [
    %imageInputLayer([28 28 1])
    imageInputLayer([120 120 1])
    
    %convolution2dLayer(3,8,'Padding','same')
    convolution2dLayer(3,16,'Padding','same')
    %batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    %convolution2dLayer(3,16,'Padding','same')
   convolution2dLayer(3,32,'Padding','same')
    %batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    %convolution2dLayer(3,32,'Padding','same')
    convolution2dLayer(3,64,'Padding','same')
    %batchNormalizationLayer
    reluLayer

    maxPooling2dLayer(2,'Stride',2)

    %convolution2dLayer(3,16,'Padding','same')
    convolution2dLayer(3,64,'Padding','same')
    %batchNormalizationLayer
    reluLayer
    
    
    %fullyConnectedLayer(10)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm', ...
    'MaxEpochs',20,...
    'InitialLearnRate',1e-4, ...
    'Verbose',false, ...
    'Plots','training-progress');
net = trainNetwork(Trainnoise,D,layers,options);
%% Classify objects
YPred = classify(net,Testset);
YValidation = transpose([0,1,0,1,0,1]);
YValidation = categorical(YValidation);
accuracy = sum(YPred == YValidation)/numel(YValidation);
