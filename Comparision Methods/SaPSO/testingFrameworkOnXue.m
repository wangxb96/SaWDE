% **********************************************************
% We are more than happy to share the related Matlab codes with you. If you do some further studies based on this research work,  please cite this study as a reference in your paper(s) as follows:
% 
%  Yu Xue, Bing Xue, and Mengjie Zhang. 2019. Self-adaptive Particle Swarm Optimization for Large-scale Feature Selection in Classification. ACM Transactions on Knowledge Discovery from Data. 1, 1, Article 1 (January 2019), 26 pages. https://doi.org/10.1145/3340848
% 
% Please notice that the reference format and the codes are not the final versions, we will continuously update them in the following some days. Please feel free to contact me if you have any questions or comments. My Email address is: xueyu@nuist.edu.cn
% 
% I will also update the reference and codes in my ResearchGate website at 
% https://www.researchgate.net/profile/Yu_Xue19?ev=hdr_xprf
% 
% (1) File "testingFrameworkOnXue" is the main procedure
% (2) The maximum number of fitness evaluation can be modified in the file "\LSFS\algLib\Benchmarks\Xue\getFunParamOnXue.m"
% (3) The codes for the proposed method are included in the file of "\LSFS\algLib\algorithms\LSFS"

beforeTesting;

classiNum1=1:10;
for classiNum2=1:1  
    classiNum=classiNum1(classiNum2);
    pro =classiNum;
  

    fname=@XueBenchmarkFun;


algorithmPool = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44];
%% Problem
% Problem = {'Alizadeh-2000-v1','Alizadeh-2000-v2','Bittner-2000','Garber-2001','West-2001','Nutt-2003-v2'};
Problem = {'Nutt-2003-v3','Pomeroy-2002-v1','Pomeroy-2002-v2','Shipp-2002-v1','Armstrong-2002-v1','Dyrskjot-2003','Liang-2005'};
 for curAlgorithm=44:44
    algori=curAlgorithm ;
    algorithm = selectAlgorithm(algorithmPool(curAlgorithm));  


%       dataNum1=1:15;
%       for dataNum2=14:14
%           
%           dataNum=dataNum1(dataNum2)
%         
%           switch dataNum
%                 case 1
%                     Particle_Number=100; D = 5000;
%                     load('./DataSet/gisette_valid.txt');       
%                     dataset=gisette_valid;
%                  case 2
%                     Particle_Number=100; D = 1300;
%                     load('./DataSet/MicroMass.txt');          
%                     dataset=MicroMass;     
%                 case 3
%                     Particle_Number=100; D = 856;
%                     load('./DataSet/CNAE.txt');          
%                     dataset=CNAE;
%                 case 4
%                     Particle_Number=100; D = 301;
%                     load('./DataSet/grammatical_facial_expression01.txt');         
%                     dataset=grammatical_facial_expression01;    
%                  case 5
%                     Particle_Number=100; D = 256;
%                     load('./DataSet/SemeionHandwrittenDigit.txt');         
%                     dataset=SemeionHandwrittenDigit;                     
%                  case 6
%                     Particle_Number=100; D = 617;  
%                     load('./DataSet/isolet5.txt');      
%                     dataset=isolet5;
%                  case 7
%                     Particle_Number=100; D = 649;
%                     load('./DataSet/MultipleFeaturesDigit.txt');         
%                     dataset=MultipleFeaturesDigit;                             
%                  case 8
%                     Particle_Number=100; D = 561;
%                     load('./DataSet/HAPTDataSet.txt');        
%                     dataset=HAPTDataSet;   
%                  case 9
%                     Particle_Number=100; D =561;
%                     load('./DataSet/har.txt');          
%                     dataset=har;                       
%                  case 10
%                     Particle_Number=100; D = 522;
%                     load('./DataSet/UJIIndoorLoc.txt');         
%                     dataset=UJIIndoorLoc;                      
%                  case 11
%                     Particle_Number=100; D = 500;
%                     load('./DataSet/MadelonValid.txt');         
%                     dataset=MadelonValid;                     
%                  case 12
%                     Particle_Number=100; D = 64;
%                     load('./DataSet/OpticalRecognitionofHandwritten.txt');          
%                     dataset=OpticalRecognitionofHandwritten;                    
%                  case 13
%                     Particle_Number=100; D = 60;
%                     load('./DataSet/ConnectionistBenchData.txt');         
%                     dataset=ConnectionistBenchData;                     
%                  case 14
%                     Particle_Number=100; D = 30;
%                     load('./DataSet/wdbc.txt');           
%                     dataset=wdbc;                      
%                  case 15
%                     Particle_Number=100; D =56;
%                     load('./DataSet/LungCancer.txt');          
%                     dataset=LungCancer;         
%                    
%           end
      for n = 1 : length(Problem)   
        tic;  
        p_name = Problem{n};
        results.p_name = p_name;   
        dataset = load(['C:\Users\c\Desktop\SaWDE\train\',p_name]);
        dataset = dataset.Train;
        feat=dataset(:,1:end-1); labels=dataset(:,end);
        D = size(feat,2);
        Particle_Number=100;
        dataNum = 1;
        runtimes2=1:30;        
        for runtimes1 = 1:1     
            runtimies=runtimes2(runtimes1)
            [  XRmin, XRmax, Lbound, Ubound, oo, maxfeval ] = getFunParamOnXue( D, pro );

            Max_Gen = ceil(maxfeval/Particle_Number); 
            

            algHdl = str2func(algorithm);

            [onef, gbestval, outGeEvaNumGbest,outStrSelPRecord]= ...
                feval(algHdl,fname,Max_Gen,maxfeval,Particle_Number,D,classiNum,dataset,pro);

        path2=sprintf('%02d%02d%02d%02d',classiNum,algori,dataNum,runtimies);
        ppath2=strcat(path2,'outGeEvaNumGbest.mat');  
        save(ppath2,'outGeEvaNumGbest');          
        outGeEvaNumGbest(end,4:5)

            
        end
        results.trainacc = 1 - outGeEvaNumGbest{end,5};
        results.selectedfeatures = outGeEvaNumGbest{end,4};
        testdataset = load(['C:\Users\c\Desktop\SaWDE\test\',p_name]);
        testdataset = testdataset.Test;
        results.testacc = 1 - Fun(outGeEvaNumGbest{end,3},testdataset,1);
        toc;
        time = num2str(toc);
        disp(time);
        results.time = time;
        saveResults(results);
       end 
      end %dataNum  
              
     end %algori 
    
 %end % classiNum

% pause
% afterTesting;

