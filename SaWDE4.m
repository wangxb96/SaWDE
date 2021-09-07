%**************************************************************************************************
% Reference: 
% G. Wu, R. Mallipeddi, P. N. Suganthan, R. Wang, and H. Chen, "Differential evolution with multi-population based ensemble of mutation strategies," Information Sciences, vol. 329, pp. 329-345, 2016.
%**************************************************************************************************

clc;
clear all;

format long;
format compact;

'SaWDE'

%% Problem
%Problem = {'Alizadeh-2000-v1','Alizadeh-2000-v2','Bittner-2000','Garber-2001','West-2001'};
Problem = {'MadelonValid.txt','OpticalRecognitionofHandwritten.txt'};

%% Objective function
fun=@jFitnessFunction1;
fun2 = @jFitnessFunction2;
% fun=@fmeasure;

strParameterDescription = 'xxpopsize = 100 ';
prepopsize = 100;
xxpopsize = 100;%popsize
xxleastSelectionPro = 0.2;
result = zeros(12,3);%column 1 for training acc, 2 for feat, 3 for test acc
selectionResult = zeros(12,1000);
% selectionResults = zeros(12,1000);
AccProcession = zeros(100, 12);
SizeProcession = zeros(100, 12);

% Additional data
for n = 1 : length(Problem)   
   tic;
p_name = Problem{n};
results.p_name = p_name;   
dataset = load(['C:\Users\c\Desktop\SaWDE\train\',p_name]);
% dataset = dataset.Train;
feat=dataset(:,1:end-1); labels=dataset(:,end);
D = size(feat,2);
D1 = D;
selectD = floor(0.5 * D); % 50% features 
if selectD < 10
   selectD = 10;
end
countp = 0;
count = 0;
% lu: define the upper and lower bounds of the variables
lu = [-1* ones(1, D); 1 * ones(1, D)];
weight1 = zeros(1,D);
weight = zeros(1,D);
time = 1;
threshold = 0.6;
% The total number of runs
totalTime = 1;
 while time <= totalTime     
   %% the values and indices of the best solutions
    FES = 0;
    leastSelectionPro = xxleastSelectionPro;
    arrayGbestChange = [1,1,1,1,1,1,1,1,1,1];
    arrayGbestChangeRate = [0,0,0,0,0,0,0,0,0,0];
    genForChange = 20;
    MaxFESpre = 1000000;
    mixPopSizePre = prepopsize;
    mixPopSize = xxpopsize;
    mixPopSizeNext = D;
    mixPopNext = zeros(mixPopSizeNext,D);
    mixValNext = zeros(mixPopSizeNext,1);
    nfeatValNext = zeros(mixPopSizeNext,1);
    fitnessValNext = zeros(mixPopSizeNext,1);
    MaxGen = round(MaxFESpre/mixPopSize);
    strategyNum = [0,0,0,0,0,0,0,0,0,0];
    strategySelect = [0,0,0,0,0,0,0,0,0,0];
    strategyRate = [0,0,0,0,0,0,0,0,0,0];
    mixPopMax = zeros(1,D);%target 1
    nfeatMax = 0;
    mixValMax = 0;
    mixPopfitMin = zeros(1,D);%target 2
    nfeatfitMin = 0;
    mixValfitMin = 0;
    fitnessMin = 0;
    %fitness
    premixVal = zeros(mixPopSizePre,1);
%   mixPop1 = repmat(lu(1, :), mixPopSizePre, 1) + rand(mixPopSizePre, D) .* (repmat(lu(2, :) - lu(1, :), mixPopSizePre, 1));
    mixPop1 = rand(mixPopSizePre, D); 
    nfeat1 = zeros(mixPopSizePre,1);
    fitness1 = zeros(mixPopSizePre,1);
    mixPop0 = zeros(mixPopSize,D);   
    mixVal0 = zeros(mixPopSize,1);
    nfeat0 = zeros(mixPopSize,1);
    fitness0 = zeros(mixPopSize,1);
    fitness = zeros(mixPopSize,1);
    for i = 1 : floor(mixPopSizePre/3)
          IndexPre = floor(rand()*D) + 1;
          for j = 1 : IndexPre  
              Index = floor(rand()*IndexPre) + 1;
              mixPop1(i,Index) = 1;
          end
          for j = 1 : D
             if mixPop1(i,j) ~= 1
                 mixPop1(i,j) = 0;
             end
          end              
    end
    for i = floor(mixPopSizePre/3)+1 : 2 * floor(mixPopSizePre/3)
       for j = 1 : D     
           if mixPop1(i,j) > 0.6
              mixPop1(i,j) = 1;
           else
              mixPop1(i,j) = 0;
           end
       end
    end
    for i = 2 * floor(mixPopSizePre/3) + 1 : mixPopSizePre
        for j = 1 : 0.50 * D
            Index = floor(rand()*D) + 1;
            mixPop1(i,Index) = 1;
        end
        for j = 1 : D
           if mixPop1(i,j) ~= 1
              mixPop1(i,j) = 0;
           end
        end 
    end
    for i = 1 : mixPopSizePre
       premixVal(i) = fun(dataset,mixPop1(i,:),1);
       nfeat1(i) = size(find(mixPop1(i,:) == 1 ),2);
%        fitness1(i) = (0.1*(nfeat1(i)/D))-0.9*premixVal(i);
       fitness1(i) = 0.1 - 0.9 * premixVal(i);
    end
    FES = FES + mixPopSizePre;
    [~,I3]=sort(fitness1, 'ascend');
    for i = 1 : mixPopSize
       mixPop0(i,:) = mixPop1(I3(i),:);
       mixVal0(i) = premixVal(I3(i)); 
       nfeat0(i) = nfeat1(I3(i));
       fitness0(i) = fitness1(I3(i));
       if nfeat0(i) <= selectD
          if mixVal0(i) > mixValMax
             mixPopMax = mixPop0(i,:);
             mixValMax = mixVal0(i);
             nfeatMax = nfeat0(i);
%              fitness(i) = (0.1 * (nfeat0(i) / D)) - 0.9 * mixVal0(i);
             fitness1(i) = 0.1 - 0.9 * mixVal0(i);
          end
       end 
    end     
    [~,II]=sort(fitness0, 'ascend');
    mixPop = mixPop0;
    mixVal = mixVal0;
    fitness = fitness0;
    nfeat = nfeat0;
    [overallBestVal,ibest] = min(fitness0);
    mixPopfitMin = mixPop(ibest,:);
    mixValfitMin = mixVal(ibest);
    nfeatfitMin = nfeat(ibest);
    if overallBestVal < fitnessMin  
       fitnessMin = overallBestVal;
    end
    c = 1/10;
    pj = 0.1;
    Afactor = 1;
    archive.NP = mixPopSize; % the maximum size of the archive
    archive.pop = zeros(0, D); % the solutions stored in te archive
    archive.funvalues = zeros(0, 1); % the function value of the archived solutions
    FESj = 0;
    consumedFES = [1,1,1,1,1,1,1,1,1,1];
    %% F & CR Initialization
    % the five control parameter settings
    F = [0.5 1 0.6 0.9 0.5 0.9 0.6 1];
    CR = [0.1 0.2 0.9 0.8 0.9 0.1 0.8 0.2];
    goodCR=0;goodF=0;CRm1=0;Fm1=0;goodCR2=0;goodF2=0;CRm2=0;Fm2=0;goodCR4=0;goodF4=0;CRm4=0;Fm4=0;
    goodCR5=0;goodF5=0;CRm5=0;Fm5=0;goodCR6=0;goodF6=0;CRm6=0;Fm6=0;
%% Choose the best results from good offspring to evolve(CBFG mechanism)
    gen = 0;
    MaxFES = MaxFESpre;
    while  FES < MaxFES  
        gen = gen + 1;  
        permutation = randperm(mixPopSize);  
        % subpopulation = 4
        arrayFirst = permutation(1:25);
        arraySecond = permutation(26:50);
        arrayThird = permutation(51:75);
        arrayFourth = permutation(76:100);
        
%         arrayThird= permutation(1:leastSelectionPro*mixPopSize);
%         arraySecond = permutation(leastSelectionPro*mixPopSize+1: 2*leastSelectionPro*mixPopSize);
%         arrayFirst = permutation(2*leastSelectionPro*mixPopSize+1: 2*leastSelectionPro*mixPopSize + leastSelectionPro*mixPopSize);
%         arrayFourth = permutation(2*leastSelectionPro*mixPopSize + leastSelectionPro*mixPopSize+1: 2*leastSelectionPro*mixPopSize + 2*leastSelectionPro*mixPopSize);
%         arrayFifth = permutation(2*leastSelectionPro*mixPopSize + 2*leastSelectionPro*mixPopSize+1:end);            
%         if mixPopSize<20
%            arrayFirst =  permutation;    
%            arraySecond = [];
%            arrayThird = [];
%            arrayFourth = [];
%            arrayFifth = [];
%         end   
    %% Muti-population mechanism 
        a1 = arrayFirst;
        a2 = arraySecond;
        a3 = arrayThird;  
        a4 = arrayFourth;
%         a5 = arrayFifth;
        paraIndex = floor(rand(1, 8) * length(F)) + 1;        
        for k = 1 : 4
    %% a selection for CoDE
          if k == 1
             array = a1;
          elseif k == 2
             array = a2;
          elseif k == 3
             array = a3;          
          elseif k == 4
             array = a4;            
          elseif k == 5    
             array = a5;  
          end
          if isempty(array) 
              continue;
          end
          if FES <= 500000
             strategy = floor(rand() * 10) + 1;
             strategyNum(strategy) = strategyNum(strategy) + 1;
          else
              strategyRate = strategySelect ./ strategyNum;
              [~,r] = sort(strategyRate,'descend');
              t = floor(rand() * 5) + 1;%select from the top 5 strategies
              strategy = r(t);
              strategyNum(strategy) = strategyNum(strategy) + 1;
          end
          if mod(gen,genForChange) == 0
             arrayGbestChangeRate(1) = arrayGbestChange(1)/consumedFES(1);
             arrayGbestChangeRate(2) = arrayGbestChange(2)/consumedFES(2);
             arrayGbestChangeRate(3) = arrayGbestChange(3)/consumedFES(3);
             arrayGbestChangeRate(4) = arrayGbestChange(4)/consumedFES(4);
             arrayGbestChangeRate(5) = arrayGbestChange(5)/consumedFES(5);
             arrayGbestChangeRate(6) = arrayGbestChange(6)/consumedFES(6);
             arrayGbestChangeRate(7) = arrayGbestChange(7)/consumedFES(7);
             arrayGbestChangeRate(8) = arrayGbestChange(8)/consumedFES(8);
             arrayGbestChangeRate(9) = arrayGbestChange(9)/consumedFES(9);
             arrayGbestChangeRate(10) = arrayGbestChange(10)/consumedFES(10);            
             [~,indexBestLN]=max(arrayGbestChangeRate);
             if sum(arrayGbestChangeRate == arrayGbestChangeRate(1)) == 10
                 indexBestLN = randi([1,10],1);
                 sprintf('changed');
             end
             strategySelect(indexBestLN) = strategySelect(indexBestLN) + 1;
             arrayGbestChange = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
             arrayGbestChangeRate =  [0,0,0,0,0,0,0,0,0,0];
             consumedFES = [1,1,1,1,1,1,1,1,1,1];    
             if indexBestLN == 1               
                strategy = 1;
             elseif indexBestLN == 2
                strategy = 2;
             elseif indexBestLN == 3
                strategy = 3;
             elseif indexBestLN == 4
                strategy = 4;
             elseif indexBestLN == 5
                strategy = 5;
             elseif indexBestLN == 6
                strategy = 6;            
             elseif indexBestLN == 7
                strategy = 7;            
             elseif indexBestLN == 8
                strategy = 8;            
             elseif indexBestLN == 9
                strategy = 9;            
             elseif indexBestLN == 10
                strategy = 10;            
             end
          end          
          switch strategy
              case 1
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR2,goodF2,CRm2,Fm2,goodCR4,goodF4,CRm4,Fm4,weight] = ...
                  CoDE1(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR2,goodF2,CRm2,Fm2,goodCR4,goodF4,CRm4,Fm4,paraIndex,D,D1,dataset,gen,array,weight); 
                  [~,co1] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co1(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end
              case 2
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR2,goodF2,CRm2,Fm2,goodCR5,goodF5,CRm5,Fm5,weight] = ...
                  CoDE2(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR2,goodF2,CRm2,Fm2,goodCR5,goodF5,CRm5,Fm5,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co2] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co2(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end              
              case 3
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR2,goodF2,CRm2,Fm2,goodCR6,goodF6,CRm6,Fm6,weight] = ...
                  CoDE3(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR2,goodF2,CRm2,Fm2,goodCR6,goodF6,CRm6,Fm6,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co3] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co3(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end              
              case 4
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR4,goodF4,CRm4,Fm4,goodCR5,goodF5,CRm5,Fm5,weight] = ...
                  CoDE4(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR4,goodF4,CRm4,Fm4,goodCR5,goodF5,CRm5,Fm5,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co4] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co4(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end              
              case 5
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR4,goodF4,CRm4,Fm4,goodCR6,goodF6,CRm6,Fm6,weight] = ...
                  CoDE5(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR4,goodF4,CRm4,Fm4,goodCR6,goodF6,CRm6,Fm6,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co5] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co5(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end              
              case 6
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR5,goodF5,CRm5,Fm5,goodCR6,goodF6,CRm6,Fm6,weight] = ...
                  CoDE6(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR5,goodF5,CRm5,Fm5,goodCR6,goodF6,CRm6,Fm6,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co6] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co6(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end              
              case 7
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR2,goodF2,CRm2,Fm2,goodCR4,goodF4,CRm4,Fm4,goodCR5,goodF5,CRm5,Fm5,weight] = ...
                  CoDE7(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR2,goodF2,CRm2,Fm2,goodCR4,goodF4,CRm4,Fm4,goodCR5,goodF5,CRm5,Fm5,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co7] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co7(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end             
              case 8
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR2,goodF2,CRm2,Fm2,goodCR4,goodF4,CRm4,Fm4,goodCR6,goodF6,CRm6,Fm6,weight] = ...
                  CoDE8(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR2,goodF2,CRm2,Fm2,goodCR4,goodF4,CRm4,Fm4,goodCR6,goodF6,CRm6,Fm6,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co8] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co8(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end              
              case 9
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR2,goodF2,CRm2,Fm2,goodCR5,goodF5,CRm5,Fm5,goodCR6,goodF6,CRm6,Fm6,weight] = ...
                  CoDE9(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR2,goodF2,CRm2,Fm2,goodCR5,goodF5,CRm5,Fm5,goodCR6,goodF6,CRm6,Fm6,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co9] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co9(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end              
              case 10
                  [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR4,goodF4,CRm4,Fm4,goodCR5,goodF5,CRm5,Fm5,goodCR6,goodF6,CRm6,Fm6,weight] = ...
                  CoDE10(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR4,goodF4,CRm4,Fm4,goodCR5,goodF5,CRm5,Fm5,goodCR6,goodF6,CRm6,Fm6,paraIndex,D,D1,dataset,gen,array,weight);
                  [~,co10] = sort(mixVal(array), 'descend');
                  for i = 1 : floor(0.2 * length(array))
                     A = find(mixPop0(co10(i),:) == 1 );
                     num = size(A,2);
                     for j = 1 : num
                        weight1(1,A(j)) = weight1(1,A(j)) + 1;
                     end
                  end
          end
        end
%         if mod(gen,100) == 0 %disp the results per 9k fes
%             count = count + 1;
%             disp(count);
%             disp('ibestVal');
%             disp(nfeatfitMin);
%             disp(mixValfitMin);
%             disp('ValMax');
%             disp(nfeatMax);
%             disp(mixValMax);
%         end
      if mod(gen,20) == 0 %change 20% popsize per 3w fes
           [~,change] = sort(mixVal,'ascend');%mixPopSize
           [~,select] = sort(weight,'descend');%D
           [~,select2] = sort(weight1,'descend');           
           for i = 1 :  floor(mixPopSizeNext/2)
               mixPopNext(i,:) = zeros(1,D);  
               for j = 1 : i
                   mixPopNext(i,select(j)) = 1;
               end
               mixValNext(i) = fun(dataset,mixPopNext(i,:),1);
               nfeatValNext(i) = size(find(mixPopNext(i,:) == 1 ),2);                  
               if nfeatValNext(i) <= selectD
                  if mixValNext(i) > mixValMax
                     mixPopMax = mixPopNext(i,:);
                     mixValMax = mixValNext(i);
                     nfeatMax = nfeatValNext(i);
                  end              
               end
               if FES >= 500000
                  fitnessValNext(i) = (0.1*(nfeatValNext(i)/D))-0.9*mixValNext(i);         
               else
                  fitnessValNext(i) = 0.1 - 0.9 * mixValNext(i);
               end
           end
           for i = floor(mixPopSizeNext / 2) + 1 :  mixPopSizeNext               
                mixPopNext(i,:) = zeros(1,D);
               for j = 1 : i - floor(mixPopSizeNext / 2)
                   mixPopNext(i,select2(j)) = 1;
               end
               mixValNext(i) = fun(dataset,mixPopNext(i,:),1);
               nfeatValNext(i) = size(find(mixPopNext(i,:) == 1 ),2);                  
               if nfeatValNext(i) <= selectD
                  if mixValNext(i) > mixValMax
                     mixPopMax = mixPopNext(i,:);
                     mixValMax = mixValNext(i);
                     nfeatMax = nfeatValNext(i);
                  end
               end
               if FES >= 500000
                  fitnessValNext(i) = (0.1*(nfeatValNext(i)/D))-0.9*mixValNext(i);         
               else
                  fitnessValNext(i) = 0.1 - 0.9 * mixValNext(i);
               end               
           end
           FES = FES + mixPopSizeNext;
           [~,changeNext] = sort(mixValNext,'descend');
%            for k = 1 : 1
%                if mixValNext(changeNext(k)) > mixVal(change(k))
%                    mixVal(change(k)) = mixValNext(changeNext(k));
%                    mixPop(change(k)) = mixPopNext(changeNext(k));
%                    nfeat(change(k)) = nfeatValNext(changeNext(k));
%                    fitness(change(k)) = fitnessValNext(changeNext(k));
%                    if fitness(change(k)) < fitnessMin
%                       mixPopfitMin = mixPop(change(k));
%                       mixValfitMin = mixVal(change(k));
%                       nfeatfitMin = nfeat(change(k));
%                       fitnessMin = fitness(change(k));
%                    end
%                end
%            end
%         end
          k = 1; j = mixPopSize;
           while k <= mixPopSizeNext
               if mixValNext(changeNext(k)) > mixVal(change(j))
                   mixVal(change(j)) = mixValNext(changeNext(k));
                   mixPop(change(j)) = mixPopNext(changeNext(k));
                   nfeat(change(j)) = nfeatValNext(changeNext(k));
                   fitness(change(j)) = fitnessValNext(changeNext(k));
                   if fitness(change(j)) < fitnessMin
                      mixPopfitMin = mixPop(change(j));
                      mixValfitMin = mixVal(change(j));
                      nfeatfitMin = nfeat(change(j));
                      fitnessMin = fitness(change(j));
                   end
                   k = k + 1;
                   j = j + 1;
               else
                   j = j + 1;
               end
               if j > mixPopSize
                   break;
               end
           end
        [outcome, ibest] = min(fitness);
        if mixVal(ibest) > mixValfitMin
            mixPopfitMin = mixPop(ibest,:);
            mixValfitMin = mixVal(ibest);
            nfeatfitMin = nfeat(ibest);
            fitnessMin = outcome;
        end
%         if(mixValfitMin == 1 && nfeatfitMin == 1)
%             break;
%         end
%         weight = zeros(1,D);
%         weight1 = zeros(1,D);
      end
      if mod(gen, 100) == 0 %record it per 100 gen
         countp = countp + 1;
         if mixValMax > mixValfitMin %case 1
                AccProcession(countp, n) =  mixValMax;
                SizeProcession(countp, n) = nfeatMax;
         elseif mixValMax == mixValfitMin%case 2
                if nfeatMax > nfeatfitMin
                    AccProcession(countp, n) =  mixValfitMin;
                    SizeProcession(countp, n) = nfeatfitMin;
                else
                    AccProcession(countp, n) =  mixValMax;
                    SizeProcession(countp, n) = nfeatMax;
                end
         else  %case 3
                AccProcession(countp, n) =  mixValfitMin;
                SizeProcession(countp, n) = nfeatfitMin;    
         end      
      end
    end
    time = time + 1;     
  end
if mixValMax > mixValfitMin %case 1
        disp(n);
        disp(nfeatMax);
        disp(mixValMax);
%         selectionResult(n,i) = mixPopMax(1,i);
        for i = 1 : D
            if mixPopMax(1,i) == 1
                selectionResult(n,i) = 1; 
            end    
        end 
        selectionResults = mixPopMax(1,:);
        result(n, 1) = mixValMax;
        result(n, 2) = nfeatMax;
 elseif mixValMax == mixValfitMin%case 2
        if nfeatMax > nfeatfitMin
            disp(n);
            disp(nfeatfitMin);
            disp(mixValfitMin);
            for i = 1 : D
                if mixPopfitMin(1, :) == 1
                    selectionResult(n,i) = 1; 
                end    
            end 
            selectionResults = mixPopfitMin(1, :);
            result(n, 1) = mixValfitMin;
            result(n, 2) = nfeatfitMin;
        else
            disp(n);
            disp(nfeatMax);
            disp(mixValMax);
            for i = 1 : D
                if mixPopMax(1,i) == 1
                    selectionResult(n,i) = 1; 
                end    
            end
            selectionResults = mixPopMax(1,:);
            result(n, 1) = mixValMax;
            result(n, 2) = nfeatMax;
        end
 else  %case 3
        disp(n);
        disp(nfeatfitMin);
        disp(mixValfitMin); 
        for i = 1 : D
            if mixPopfitMin(1,i) == 1
                selectionResult(n,i) = 1; 
            end    
        end 
        selectionResults = mixPopfitMin(1, :);
        result(n, 1) = mixValfitMin;
        result(n, 2) = nfeatfitMin;      
end 
results.trainacc = mixValfitMin;
results.selectedfeatures = nfeatfitMin;
testdataset = load(['C:\Users\c\Desktop\SaWDE\test\',p_name]);
% testdataset = testdataset.Test;
fit = fun2(dataset, testdataset, mixPopfitMin ,1);
results.testacc = fit;
toc;
time = num2str(toc);
disp(time);
results.time = time;
saveResults(results);
end

