function f = jFitnessFunction1(dataset,X,func_num)
% subx=find(X==1);
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


function Fit = Fit_classifyCV(data, dataLab, indices, classif)  
CVF = max(indices);
Fit = zeros(CVF,1);  
for k=1:CVF
    testn = (indices == k); 
    trainn = ~testn;       
    NTest = sum(testn);   
%     predicted_labels = zeros(size(testn, 1), 1);
    switch classif
        case 'KNN'
            mdl= fitcknn(data(trainn,:),dataLab(trainn),'NumNeighbors',3);
            predicted_labels = predict(mdl,data(testn,:));
        case 'LDA'
            mdl = fitcdiscr(data(trainn,:),dataLab(trainn,:));
            predicted_labels = predict(mdl,data(testn,:));
        case 'RegTree'
            mdl = fitctree(data(trainn),dataLab(trainn,:));
            predicted_labels = predict(mdl,data(testn,:));
        case 'NB'
            mdl = fitcnb(data(trainn),dataLab(trainn,:));
            predicted_labels = predict(mdl,data(testn,:));
        case 'SVM'
            mdl = fitcsvm(data(trainn),dataLab(trainn,:));
            predicted_labels = predict(mdl,data(testn,:));
    end
    if size(predicted_labels,1) == sum(testn)
        Fit(k) = sum(predicted_labels~=dataLab(testn))/NTest; 
    else
        Fit(k) = 1; 
    end
%     p = dataLab(testn,:) == 1;
%     tp = sum(predicted_labels == 1 & p);
%     fp = sum(predicted_labels == 1 & ~p);
%     fn = sum(predicted_labels == 0 & p);
%     precision = tp / (tp + fp);
%     recall = tp / (tp + fn);
%     Fit(k) = 2 * (precision * recall) / (precision + recall);
%     if isnan(Fit(k))
%         Fit(k) = 0;
%     end
end    
Fit = mean(Fit);  





