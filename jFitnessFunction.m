function f=jFitnessFunction(dataset,X,func_num)
if func_num==1
  classif='KNN';
elseif func_num==2
  classif='LDA';
elseif func_num==3
  classif='RegTree';
elseif func_num==4
  classif= 'NB';
elseif func_num==5
  classif='SVM';
end

FeatIndexKNN = find(X==1); %Feature Index
if ~any(X(:))
    FeatIndexKNN(1) = 1;        %Check that at least one feature is selected
end 

feat = dataset(:,1:end-1);
data_tr = feat(:,[FeatIndexKNN]);
data_tr=[data_tr dataset(:,end)];

[Ninst, NF] = size(data_tr);
NF = NF - 1;                

CVF = 3;
indices = crossvalind('Kfold',Ninst,CVF);
fac = 0.00000000000000000000001;
f = Fit_classifyCV(data_tr(:,1:NF), data_tr(:,end), indices, classif) + fac;
end

function [test_Predicted_labels] = knn(Xtrain,Ltrain,Xtest, K)
%% This is the matlab implemenatation for K nearest neighbors
% Input Arguments: 
% Xtrain : training data set
% Ltrain : Labels of training samples
% Xtest : test data set
% K : number of nearest neighbors 
% Output Arguments :
% TestLabel : Predicted labels of the output data set
% default value of K if not given by user is 8.
if(nargin < 3)
error('Incorrect number of inputs.');
end
if(nargin < 4)
   K = 8; 
end
[N , ~] = size(Xtrain);
[Ntest,~] = size(Xtest);
distance = zeros(N,Ntest);
%descendingDistances = zeros(N,Nt);
%Ltest = repmat(Xtest(1,:),N,1);

% calculating the euclidean distance of the test samples from training
% samples
for i = 1: Ntest
     for j = 1: N 
distance(j,i) = norm(Xtest(i,:)-Xtrain(j,:));
     end
end

% ascendingdistances stores all the distances of the test samples
% from the all training samples in cloumns
% Index will have indices of the corresponding training sample
[~,Index]= sort(distance,'ascend');

% consider only top K nearest neighbors to predict the label for test
% sample
Ltest = zeros(K,Ntest);
for i = 1:Ntest
    for j=1:K
    Ltest(j,i) = Ltrain(Index(j,i));
    end
end

test_Predicted_labels = mode(Ltest);
test_Predicted_labels = test_Predicted_labels';

end


function[Error] =  Misclassification_error(test_labels, test_Predicted_labels)
%% this is the matlab Implementation for KNN misclassification error
% if you have labels for your test data set
% use KNN.m to find the predicted labels for test samples.
% then you can use this function to find misclassification error for your 
% test data set
% Input arguments: 
% test_labels : Known labels of your test data
% test_predicted_labels : predicted test labels

if(nargin < 2)
error('Incorrect number of inputs.');
end

% check if the input labels are in row vector
% if yes make it column vectors
if(isrow(test_labels)==1)
    test_labels = test_labels';
end

if(isrow(test_Predicted_labels)==1)
    test_Predicted_labels = test_Predicted_labels';
end

% matrix dimensions should be equal
[Ltest,~] = size(test_labels);
if (Ltest ~= size(test_Predicted_labels,1))
   Error('matrix dimensions are not consistent');
end

count=0;

for i=1:Ltest   
if(test_labels(i,1)==test_Predicted_labels(i,1))
    count = count+1;
end    
end

missPredictions = count;
Error = missPredictions / Ltest;

end


function Fit = Fit_classifyCV(data, dataLab, indices, classif)  
CVF = max(indices);  
Fit = zeros(CVF,1);  
for k=1:CVF
    testn = (indices == k); 
    trainn = ~testn;       
    %NTest = sum(testn);   
    switch classif
        case 'KNN'
            test_Predicted_labels = knn(data(trainn,:),dataLab(trainn),data(testn,:),3);
            Fit(k) = Misclassification_error(dataLab(testn),test_Predicted_labels);      
    end
end
Fit =1 - mean(Fit);  
end


 
