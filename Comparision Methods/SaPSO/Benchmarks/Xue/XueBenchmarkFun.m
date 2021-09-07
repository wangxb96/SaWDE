function f=XueBenchmarkFun(x,dataset,func_num)

subx = find(x==1);
if func_num==1
classif='KNN';
end


data_tr=dataset(:,subx);  
data_tr=[data_tr dataset(:,end)];

[Ninst, NF] = size(data_tr);
NF = NF - 1;                

CVF = 3;      
indices = crossvalind('Kfold',data_tr(:,end),CVF);

fac = 0.00000000000000000000001;
f = Fit_classifyCV(data_tr(:,1:NF), data_tr(:,end), indices, classif) + fac; 

function Fit = Fit_classifyCV(data, dataLab, indices, classif)  
CVF = max(indices);  
Fit = zeros(CVF,1);  
for k=1:CVF,   
    testn = (indices == k); 
    trainn = ~testn;       
    NTest = sum(testn);   


    switch classif
        case 'KNN'

            mdl=ClassificationKNN.fit(data(trainn,:),dataLab(trainn),'numneighbors',3);  %used for windows
            Ac1=predict(mdl,data(testn,:));
%             Ac1=knnclassify(data(testn,:),data(trainn,:),dataLab(trainn),3);
%             Used for Linux
            if size(Ac1,1) == sum(testn),
                Fit(k) = sum(Ac1~=dataLab(testn))/NTest; 

            else
                Fit(k) = 1; 
            end

    end
end
Fit = mean(Fit);  



 
