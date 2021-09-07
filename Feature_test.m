clc;
clear all;
tic;

format long;
format compact;

'SaCoDE'

% Objective function
fun=@jFitnessFunction2;
load('outGeEvaNumGbest1_3.mat');
Feature = seletionResult;
for n = 3 :3 
   if n == 1 
        %Data 1
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\grammatical_facial_expression01.txt');
        dataset_train=grammatical_facial_expression01;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\grammatical_facial_expression01.txt');
        dataset_test=grammatical_facial_expression01;
        
   elseif n == 2  
        %Data 2
         load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\SemeionHandwrittenDigit.txt');
        dataset_train=SemeionHandwrittenDigit;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\SemeionHandwrittenDigit.txt');
        dataset_test=SemeionHandwrittenDigit;
        %load('SemeionHandwrittenDigit.txt');
        %dataset=SemeionHandwrittenDigit;
   elseif n == 3  
        %Data 3
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\isolet5.txt');
        dataset_train=isolet5;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\isolet5.txt');
        dataset_test=isolet5;       
   elseif n == 4  
        %Data 4
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\MultipleFeaturesDigit.txt');
        dataset_train=MultipleFeaturesDigit;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\MultipleFeaturesDigit.txt');
        dataset_test=MultipleFeaturesDigit;
   elseif n == 5   
        %Data 5
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\HAPTDataSet.txt');
        dataset_train=HAPTDataSet;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\HAPTDataSet.txt');
        dataset_test=HAPTDataSet;
   elseif n == 6
        %Data 6
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\har.txt');
        dataset_train=har;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\har.txt');
        dataset_test=har;
   elseif n == 7  
        %Data 7
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\UJIIndoorLoc.txt');
        dataset_train=UJIIndoorLoc;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\UJIIndoorLoc.txt');
        dataset_test=UJIIndoorLoc;
   elseif n == 8  
        %Data 8
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\MadelonValid.txt');
        dataset_train=MadelonValid;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\MadelonValid.txt');
        dataset_test=MadelonValid;
   elseif n == 9 
        %Data 9
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\OpticalRecognitionofHandwritten.txt');
        dataset_train=OpticalRecognitionofHandwritten;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\OpticalRecognitionofHandwritten.txt');
        dataset_test=OpticalRecognitionofHandwritten;
   elseif n == 10  
        %Data 10
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\ConnectionistBenchData.txt');
        dataset_train=ConnectionistBenchData;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\ConnectionistBenchData.txt');
        dataset_test=ConnectionistBenchData;
   elseif n == 11  
        %Data 11
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\wdbc.txt');
        dataset_train=wdbc;
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\wdbc.txt');
        dataset_test=wdbc;
   elseif n == 12  
        %Data 12
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\LungCancer.txt');
        dataset_train=LungCancer; 
        
        load('C:\Users\Administrator\Desktop\Íø¿Î\SaCoDE\test\LungCancer.txt');
        dataset_test=LungCancer;
   end
%    feat=dataset(:,1:end-1); labels=dataset(:,end);
%    X_train = feat;
%    Y_train = labels;
%    % number of features
%    numF = size(X_train,2);
%    % Infinite Latent Feature Selection - ICCV 2017
%    [ranking, weights] = ILFS(X_train, Y_train , 6, 0 );
%    k = round(0.5 * numF); 
%    data_tr = X_train(:,ranking(1:k));
%    dataset = [data_tr labels];
   fit = zeros(1,100);        

   for i = 1 :100
       fit(i) = fun(dataset_train,dataset_test,Feature(n,:),1);
   end
   disp(mean(fit));

end  