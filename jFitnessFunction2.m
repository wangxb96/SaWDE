function f=jFitnessFunction2(dataset_train,dataset_test,X,func_num)
% subx=find(X==1);
if func_num==1
classif='KNN';
end

FeatIndexKNN = find(X==1); %Feature Index
if ~any(X(:))
    FeatIndexKNN(1) = 1;        %Check that at least one feature is selected
end 

feat = dataset_train(:,1:end-1);
data_tr = feat(:,[FeatIndexKNN]);
% data_tr=dataset(:,subx);  
data_tr=[data_tr dataset_train(:,end)];

[Ninst, NF] = size(data_tr);
NF = NF - 1;                

CVF = 3;
indices = crossvalind('Kfold',Ninst,CVF); 
fac = 0.00000000000000000000001;
f = Fit_classifyCV(data_tr(:,1:NF), data_tr(:,end), indices, classif,dataset_test,FeatIndexKNN) + fac;


function Fit = Fit_classifyCV(data, dataLab, indices, classif,dataset_test,FeatIndexKNN)  
CVF = max(indices);
for k=1:CVF
    testn = (indices == k); 
    trainn = ~testn;        
    switch classif
        case 'KNN'
            mdl= fitcknn(data(trainn,:),dataLab(trainn),'NumNeighbors',3);
    end
end

feat = dataset_test(:,1:end-1);
data_tr = feat(:,[FeatIndexKNN]);

Ac2=predict(mdl,data_tr);
matchKNN = find(Ac2 ~= dataset_test(:,end));
Fit = size(matchKNN,1)/size(data_tr,1);

Fit = 1 - mean(Fit);  




