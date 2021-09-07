%**************************************************************************************************
% Developed by WANG, Xubin School of Artificial Intelligence, Jilin University,
% Changchun, China
% wangxb19@mails.jlu.edu.cn

% Reference: 
% G. Wu, R. Mallipeddi, P. N. Suganthan, R. Wang, and H. Chen, "Differential evolution with multi-population based ensemble of mutation strategies," Information Sciences, vol. 329, pp. 329-345, 2016.
%**************************************************************************************************

clc;
clear all;
tic;

format long;
format compact;

'SaWDE'
  
% Objective function
fun=@jFitnessFunction1;
% fun=@fmeasure;  

strParameterDescription = 'xxpopsize = 100 ';
prepopsize = 100;
xxpopsize = 100;%popsize
xxleastSelectionPro = 0.2;
selectedFeature = zeros(12, 1000);
AccResults = zeros(12, 2);
AccProcessionCSGS = zeros(100, 12);
SizeProcessionCSGS = zeros(100, 12);
for CSGS = 8 : 8
%load data
  for n = 6 : 6 
   l = 1;
   if n == 1 
        %Data 1
        load('grammatical_facial_expression01.txt');
        dataset=grammatical_facial_expression01;
   elseif n == 2  
        %Data 2
        load('SemeionHandwrittenDigit.txt');
        dataset=SemeionHandwrittenDigit;
   elseif n == 3  
        %Data 3
        load('isolet5.txt');
        dataset=isolet5;
   elseif n == 4  
        %Data 4
        load('MultipleFeaturesDigit.txt');
        dataset=MultipleFeaturesDigit;
   elseif n == 5   
        %Data 5
        load('HAPTDataSet.txt');
        dataset=HAPTDataSet;
   elseif n == 6
        %Data 6
        load('har.txt');
        dataset=har;
   elseif n == 7  
        %Data 7
        load('UJIIndoorLoc.txt');
        dataset=UJIIndoorLoc;
   elseif n == 8  
        %Data 8
        load('MadelonValid.txt');
        dataset=MadelonValid;
   elseif n == 9 
        %Data 9
        load('OpticalRecognitionofHandwritten.txt');
        dataset=OpticalRecognitionofHandwritten;
   elseif n == 10  
        %Data 10
        load('ConnectionistBenchData.txt');
        dataset=ConnectionistBenchData;
   elseif n == 11  
        %Data 11
        load('wdbc.txt');
        dataset=wdbc;
   elseif n == 12  
        %Data 12
        load('LungCancer.txt');
        dataset=LungCancer;
   end 
feat=dataset(:,1:end-1); labels=dataset(:,end);
D = size(feat,2);
countp = 0;
% lu: define the upper and lower bounds of the variables
lu = [-1* ones(1, D); 1 * ones(1, D)];
%fitness
premixVal = zeros(prepopsize,1);
mixVal0 = zeros(xxpopsize,1);
time = 1;
threshold = 0.6;
% The total number of runs
totalTime = 1;
 while time <= totalTime     
   %% the values and indices of the best solutions
    FES = 0;
    leastSelectionPro = xxleastSelectionPro;
    MaxFESpre = 1000000;
    mixPopSizePre = prepopsize;
    mixPopSize = xxpopsize;
    MaxGen = round(MaxFESpre/mixPopSize);  

    %------Evaluate the best member after initialization------%
    mixPop0 = zeros(mixPopSize,D);
%     mixPop1 = repmat(lu(1, :), mixPopSizePre, 1) + rand(mixPopSizePre, D) .* (repmat(lu(2, :) - lu(1, :), mixPopSizePre, 1));
    mixPop1 = rand(mixPopSizePre, D);
    mixPop2 = mixPop1;
    p1 = 0.6;
    p2 = 0.6;
    for i = 1 : mixPopSizePre
      r1 = rand;
      if r1 <= p1
          for d =  1 : D
            if mixPop2(i,d) >= p2
                mixPop1(i,d) =1;
            else
                mixPop1(i,d) =0;
            end
          end
          if size(find(mixPop1(i,:)==1),2)==D   
               r3=randperm(D);
               r4=r3(1);
               mixPop1(i,r4)=1; 
               mixPop2(i,r4) = p2 + (1 - p2) * rand;
          end 
          if size(find(mixPop1(i,:)==0),2)==D   
               r3=randperm(D);
               r4=r3(1);
               mixPop1(i,r4)=0; 
               mixPop2(i,r4) = p2 * rand;
          end
       else 
          for d = 1:D
              if mixPop2(i,d) <= p2
                  mixPop1(i,d) = 1;
                  mixPop2(i,d) = p2 + (1 - p2) * rand;
              else
                  mixPop1(i,d) = 0;
                  mixPop2(i,d) = p2 * rand;
              end
          end
          if size(find(mixPop1(i,:)==0),2)==D
              r3 = randperm(D);
              r4 = r3(1);
              mixPop2(i,r4) = p2 + (1 - p2)*rand;
          end
          if size(find(mixPop1(i,:)==1),2)==D
              r3 = randperm(D);
              r4 = r3(1);
              mixPop1(i,r4) = 0;
              mixPop2(i,r4) = p2 *rand;
          end
      end  
    end
    parfor i = 1 : mixPopSizePre
       premixVal(i) = fun(dataset,mixPop1(i,:),1);
    end
    [~,I3]=sort(premixVal, 'descend');
    parfor i = 1 : mixPopSize
       mixPop0(i,:) = mixPop1(I3(i),:);
       mixVal0(i) = fun(dataset,mixPop0(i,:),1); 
    end  
    [~,II]=sort(mixVal0, 'descend');
    for i = 1:5
       disp(i);
       disp(mixVal0(II(i)));
    end    
    
    mixPop = mixPop0;
    mixVal = mixVal0;
    overallBestVal = max(mixVal0);
    c = 1/10;
    pj = 0.1;
    Afactor = 1;
    archive.NP = mixPopSize; % the maximum size of the archive
    archive.pop = zeros(0, D); % the solutions stored in te archive
    archive.funvalues = zeros(0, 1); % the function value of the archived solutions

    FESj = 0;
    %% F & CR Initialization
    % the five control parameter settings
    F = [0.5 1 0.6 0.9 0.5 0.9 0.6 1];
    CR = [0.1 0.2 0.9 0.8 0.9 0.1 0.8 0.2];
    paraIndex = floor(rand(1, 8) * length(F)) + 1;
%% Choose the best results from good offspring to evolve(CBFG mechanism)
    gen = 0;
    MaxFES = MaxFESpre;
    permutation = randperm(mixPopSize);
%     dlmwrite('genChange.txt',n,'\t')
    while  FES < MaxFES%%gen < MaxGen   
        gen = gen + 1;  
        if CSGS == 1
            arrayFirst = permutation(1:mixPopSize);
%%  %% ===========================mutation 1=====================================%%%%
           if ~isempty(arrayFirst) 
            pop1 = mixPop(arrayFirst,:); % the old population becomes the current population
            valParents1 = mixVal(arrayFirst);
            popsize = length(arrayFirst);
            [~,I1]=sort(mixVal, 'descend');%I1:the indices of xxpopsize,big->small
            [~,I2]=sort(valParents1, 'ascend');%I2:the indices of arrayFirst,small->big
            for r = 1 : 3
                pop1(I2(r),:) = mixPop(I1(r),:);%put the overall best into arrayFirst's smallest part
                valParents1(I2(r)) = mixVal(I1(r));%put the overall best value into arrayFirst's smallest part
            end
            prevalParents1 = valParents1;

            if FESj > 1 && ~isempty(goodCR) && sum(goodF) > 0 % If goodF and goodCR are empty, pause the update
                CRm1 = (1 - c) * CRm1 + c * mean(goodCR);
                Fm1 = (1 - c) * Fm1 + c * sum(goodF .^ 2) / sum(goodF); % Lehmer mean
            else
                CRm1 =  CR(paraIndex(1)); 
                Fm1 = F(paraIndex(1));               
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [Fj, CRj] = randFCR(popsize, CRm1, 0.1, Fm1, 0.1);
            r0 = [1 : popsize];
            popAll = [pop1; archive.pop];
            [r1, r2] = gnR1R2(popsize, size(popAll, 1), r0);
            % Find the p-best solutions
            [~, indBest] = sort(valParents1, 'descend');
            pNP = max(round(pj * popsize), 5); % choose at least two best solutions
            randindex = ceil(rand(1, popsize) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop1(indBest(randindex), :); % randomly choose one of the top 100p% solutions
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            %DE/current to best/1
            vi = pop1 + Fj(:, ones(1, D)) .* (pbest - pop1 + pop1(r1, :) - popAll(r2, :));
            % == == == == = Crossover == == == == =
            mask = rand(popsize, D) > CRj(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize)'; cols = floor(rand(popsize, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize D], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop1(mask);

            valOffspring1 = zeros(popsize,1);
           for i = 1 : popsize
             for j = 1 : D
                if ui(i,j) >= threshold
                   ui(i,j) = 1;
                else
                   ui(i,j) = 0;
                end   
                if size(find(ui(i,:)==0),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 1;
                end
                if size(find(ui(i,:)==1),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 0;
                end
              end  
            end
            parfor i = 1 : popsize
               valOffspring1(i) = fun(dataset,ui(i,:),1);
            end
            FESj = FESj + popsize;
            FES = FES + popsize;
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents1, I] = max([valParents1, valOffspring1], [], 2);
            popold1 = pop1;
            archive = updateArchive(archive, popold1(I == 2, :), valParents1(I == 2));
            popold1(I == 2, :) = ui(I == 2, :);
            valParents1(I == 2,:) = valOffspring1(I == 2, :);
            goodCR = CRj(I == 2);
            goodF = Fj(I == 2);
            if max(valParents1) > overallBestVal
                overallBestVal = max(valParents1);
            end
%             arrayGbestChange(1) = arrayGbestChange(1) + sum(prevalParents1- valParents1);
            if prevalParents1(I2(1)) == valParents1(I2(1)) %if nothing changed at last,restore it.
                for r = 1 : 3
                   popold1(I2(r),:) =  mixPop(arrayFirst(I2(r)),:);
                   valParents1(I2(r)) = mixVal(arrayFirst(I2(r)));
                end   
            end
            mixPop(arrayFirst,:) = popold1;
            mixVal(arrayFirst) = valParents1;
           end
        elseif CSGS == 2
            arraySecond = permutation(1:mixPopSize);
            %% ===========================mutation 2=====================================%%%%
           if ~isempty(arraySecond)
            pop2 = mixPop(arraySecond,:); % the old population becomes the current population 
            valParents2 = mixVal(arraySecond);
            popsize2 = length(arraySecond);
            [~,I1]=sort(mixVal, 'descend');
            [~,I2]=sort(valParents2, 'ascend');
            for r = 1 : 3
                pop2(I2(r),:) = mixPop(I1(r),:);
                valParents2(I2(r)) = mixVal(I1(r));
            end 
            prevalParents2 = valParents2;

            if gen > 1 && ~isempty(goodCR2) && sum(goodF2) > 0 % If goodF and goodCR are empty, pause the update
                CRm2 = (1 - c) * CRm2 + c * mean(goodCR2);
                Fm2 = (1 - c) * Fm2 + c * sum(goodF2 .^ 2) / sum(goodF2); % Lehmer mean
            else
                CRm2 =  CR(paraIndex(2)); 
                Fm2 = F(paraIndex(2));
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F2, CR2] = randFCR(popsize2, CRm2, 0.1, Fm2, 0.1);
            rot = (0:1:popsize2-1);
            ind = randperm(2);
            a1  = randperm(popsize2);             % shuffle locations of vectors
            rt = rem(rot+ind(1),popsize2);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),popsize2);
            a3  = a2(rt+1);
            pm1 = pop2(a1,:);             % shuffled population 1
            pm2 = pop2(a2,:);             % shuffled population 2
            pm3 = pop2(a3,:);             % shuffled population 3
           % ============================================mutation ===============================================
           %DE/current to rand/1
            vi =pop2 + repmat(rand(popsize2,1),1,D) .* (pm1 - pop2) + F2(:, ones(1, D)) .* (pm2 - pm3);
            ui = vi;
            valOffspring2 = zeros(popsize2,1);
            for i = 1 : popsize2
              for j = 1 : D
                if ui(i,j) >= threshold
                   ui(i,j) = 1;
                else
                   ui(i,j) = 0;
                end   
                if size(find(ui(i,:)==0),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 1;
                end
                if size(find(ui(i,:)==1),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 0;
                end
              end  
            end 
            parfor i = 1 : popsize2
               valOffspring2(i) = fun(dataset,ui(i,:),1);
            end
            FES = FES + popsize2;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents2, I] = max([valParents2, valOffspring2], [], 2);
            popold2 = pop2;
            popold2(I == 2, :) = ui(I == 2, :);
            valParents2(I == 2,:) = valOffspring2(I == 2, :);
            goodCR2 = CR2(I == 2);
            goodF2 = F2(I == 2);            
            if max(valParents2) > overallBestVal
                overallBestVal = max(valParents2);
            end 
%             arrayGbestChange(2) = arrayGbestChange(2) + sum(prevalParents2- valParents2);
            if prevalParents2(I2(1)) == valParents2(I2(1))
                for r = 1 : 3
                    popold2(I2(r),:) =  mixPop(arraySecond(I2(r)),:);
                    valParents2(I2(r)) = mixVal(arraySecond(I2(r)));
                end
            end
            mixPop(arraySecond,:) = popold2;
            mixVal(arraySecond) = valParents2;
           end
        elseif CSGS == 3
           arrayThird = permutation(1:mixPopSize);
           %% ===========================mutation 3 =====================================%%%%
           if ~isempty(arrayThird)
            pop3 = mixPop(arrayThird,:); % the old population becomes the current population 
            valParents3 = mixVal(arrayThird);
            popsize3 = length(arrayThird);
            [~,I1]=sort(mixVal, 'descend');
            [~,I2]=sort(valParents3, 'ascend');
            for r = 1 : 3
                pop3(I2(r),:) = mixPop(I1(r),:);
                valParents3(I2(r)) = mixVal(I1(r)); 
            end
            prevalParents3 = valParents3;
            if gen > 1 && ~isempty(goodCR3) && sum(goodF3) > 0 % If goodF and goodCR are empty, pause the update
                CRm3 = (1 - c) * CRm3 + c * mean(goodCR3);
                Fm3 = (1 - c) * Fm3 + c * sum(goodF3 .^ 2) / sum(goodF3); % Lehmer mean
            else
                CRm3 =  CR(paraIndex(3)); 
                Fm3 = F(paraIndex(3));
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F3, CR3] = randFCR(popsize3, CRm3, 0.1, Fm3, 0.1);
            rot = (0:1:popsize3-1);
            ind = randperm(2);
            a1  = randperm(popsize3);             % shuffle locations of vectors
            rt = rem(rot+ind(1),popsize3);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),popsize3);
            a3  = a2(rt+1);
            pm1 = pop3(a1,:);             % shuffled population 1
            pm2 = pop3(a2,:);             % shuffled population 2
            pm3 = pop3(a3,:);             % shuffled population 3
            %DE/rand/1
            vi =pm1 + F3(:, ones(1, D)) .* repmat(rand(popsize3,1),1,D) .* (pm1 - pop3 + pm2 - pop3 - pm3 - pop3); 

            mask = rand(popsize3, D) > CR3(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize3)'; cols = floor(rand(popsize3, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize3 D], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop3(mask);

            valOffspring3 = zeros(popsize3,1);
             for i = 1 : popsize3
              for j = 1 : D
                if ui(i,j) >= threshold
                   ui(i,j) = 1;
                else
                   ui(i,j) = 0;
                end   
                if size(find(ui(i,:)==0),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 1;
                end
                if size(find(ui(i,:)==1),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 0;
                end
              end  
            end
            parfor i = 1 : popsize3
               valOffspring3(i) = fun(dataset,ui(i,:),1);
            end
            FES = FES + popsize3;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents3, I] = max([valParents3, valOffspring3], [], 2);
            popold3 = pop3;
            popold3(I == 2, :) = ui(I == 2, :);
            valParents3(I == 2,:) = valOffspring3(I == 2, :);
            goodCR3 = CR3(I == 2);
            goodF3 = F3(I == 2);
            if max(valParents3) > overallBestVal
                overallBestVal = max(valParents3);
            end
%             arrayGbestChange(3) = arrayGbestChange(3) + sum(prevalParents3- valParents3);
            if prevalParents3(I2(1)) == valParents3(I2(1))
                for r = 1 : 3
                    popold3(I2(r),:) =  mixPop(arrayThird(I2(r)),:);
                    valParents3(I2(r)) = mixVal(arrayThird(I2(r)));
                end
            end 
            mixPop(arrayThird,:) = popold3;
            mixVal(arrayThird) = valParents3; 
           end 
        elseif CSGS == 4
           arrayFourth = permutation(1:mixPopSize);
          %% ===========================mutation 4 =====================================%%%%
           if ~isempty(arrayFourth)
            pop4 = mixPop(arrayFourth,:); % the old population becomes the current population 
            valParents4 = mixVal(arrayFourth);
            popsize4 = length(arrayFourth);
            [~,I1]=sort(mixVal, 'descend');
            [~,I2]=sort(valParents4, 'ascend');
            for r = 1 : 3
                pop4(I2(r),:) = mixPop(I1(r),:);
                valParents4(I2(r)) = mixVal(I1(r));  
            end
            prevalParents4 = valParents4;
            if gen > 1 && ~isempty(goodCR4) && sum(goodF4) > 0 % If goodF and goodCR are empty, pause the update
                CRm4 = (1 - c) * CRm4 + c * mean(goodCR4);
                Fm4 = (1 - c) * Fm4 + c * sum(goodF4 .^ 2) / sum(goodF4); % Lehmer mean
            else
                CRm4 =  CR(paraIndex(4)); 
                Fm4 = F(paraIndex(4));
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F4, CR4] = randFCR(popsize4, CRm4, 0.1, Fm4, 0.1);
            r0 = [1 : popsize4];
            popAll = [pop4; archive.pop];
            [r1, r2] = gnR1R2(popsize4, size(popAll, 1), r0);
            % Find the p-best solutions
            [~, indBest] = sort(valParents4, 'descend');
            pNP = max(round(pj * popsize4), 5); % choose at least two best solutions
            randindex = ceil(rand(1, popsize4) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop4(indBest(randindex), :); % randomly choose one of the top 100p% solutions
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            %-----DE/best/1
            vi = pbest + F4(:, ones(1, D)) .* (pop4(r1, :) - popAll(r2, :));          

            mask = rand(popsize4, D) > CR4(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize4)'; cols = floor(rand(popsize4, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize4 D], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop4(mask);

            valOffspring4 = zeros(popsize4,1);
             for i = 1 : popsize4
              for j = 1 : D
                if ui(i,j) >= threshold
                   ui(i,j) = 1;
                else
                   ui(i,j) = 0;
                end   
                if size(find(ui(i,:)==0),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 1;
                end
                if size(find(ui(i,:)==1),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 0;
                end
              end  
            end
            parfor i = 1 : popsize4
               valOffspring4(i) = fun(dataset,ui(i,:),1);
            end
            FES = FES + popsize4;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents4, I] = max([valParents4, valOffspring4], [], 2);
            popold4 = pop4;
            archive = updateArchive(archive, popold4(I == 2, :), valParents4(I == 2));
            popold4(I == 2, :) = ui(I == 2, :);
            valParents4(I == 2,:) = valOffspring4(I == 2, :);
            goodCR4 = CR4(I == 2);
            goodF4 = F4(I == 2);
            if max(valParents4) > overallBestVal
                overallBestVal = max(valParents4);
            end
%             arrayGbestChange(4) = arrayGbestChange(4) + sum(prevalParents4- valParents4);
            if prevalParents4(I2(1)) == valParents4(I2(1))
                for r = 1 : 3
                    popold4(I2(r),:) =  mixPop(arrayFourth(I2(r)),:);
                    valParents4(I2(r)) = mixVal(arrayFourth(I2(r)));
                end
            end 
            mixPop(arrayFourth,:) = popold4;
            mixVal(arrayFourth) = valParents4; 
           end
        elseif CSGS == 5
           arrayFifth = permutation(1:mixPopSize);
          %% ============================mutation 5 =====================================%%%%
           if ~isempty(arrayFifth)
            pop5 = mixPop(arrayFifth,:); % the old population becomes the current population 
            valParents5 = mixVal(arrayFifth);
            popsize5 = length(arrayFifth);
            [~,I1]=sort(mixVal, 'descend');
            [~,I2]=sort(valParents5, 'ascend');
            for r = 1 : 3
                pop5(I2(r),:) = mixPop(I1(r),:);
                valParents5(I2(r)) = mixVal(I1(r));   
            end
            prevalParents5 = valParents5;
            if gen > 1 && ~isempty(goodCR5) && sum(goodF5) > 0 % If goodF and goodCR are empty, pause the update
                CRm5 = (1 - c) * CRm5 + c * mean(goodCR5);
                Fm5 = (1 - c) * Fm5 + c * sum(goodF5 .^ 2) / sum(goodF5); % Lehmer mean
            else
                CRm5 =  CR(paraIndex(5)); 
                Fm5 = F(paraIndex(5));
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F5, CR5] = randFCR(popsize5, CRm5, 0.1, Fm5, 0.1);
            r0 = [1 : popsize5];
            popAll = [pop5; archive.pop];
            [r1, r2] = gnR1R2(popsize5, size(popAll, 1), r0);
            % Find the p-best solutions
            [~, indBest] = sort(valParents5, 'descend');
            pNP = max(round(pj * popsize5), 5); % choose at least two best solutions
            randindex = ceil(rand(1, popsize5) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop5(indBest(randindex), :); % randomly choose one of the top 100p% solutions
            rot = (0:1:popsize5-1);
            ind = randperm(2);
            a1  = randperm(popsize5);             % shuffle locations of vectors
            rt = rem(rot+ind(1),popsize5);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),popsize5);
            a3  = a2(rt+1);
            pm1 = pop5(a1,:);             % shuffled population 1
            pm2 = pop5(a2,:);             % shuffled population 2
            pm3 = pop5(a3,:);             % shuffled population 3
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            %DE/rand to best/1
            vi =repmat(rand(popsize5,1),1,D) .* (pm1 - pop5) + F5(:, ones(1, D)) .* pop5(r1, :) - popAll(r2, :) + F5(:, ones(1, D)) .* ( pbest - repmat(rand(popsize5,1),1,D)); 
            
            mask = rand(popsize5, D) > CR5(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize5)'; cols = floor(rand(popsize5, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize5 D], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop5(mask);

            valOffspring5 = zeros(popsize5,1);
             for i = 1 : popsize5
              for j = 1 : D
                if ui(i,j) >= threshold
                   ui(i,j) = 1;
                else
                   ui(i,j) = 0;
                end   
                if size(find(ui(i,:)==0),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 1;
                end
                if size(find(ui(i,:)==1),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 0;
                end
              end  
            end
            parfor i = 1 : popsize5
               valOffspring5(i) = fun(dataset,ui(i,:),1);
            end
            FES = FES + popsize5;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents5, I] = max([valParents5, valOffspring5], [], 2);
            popold5 = pop5;
            archive = updateArchive(archive, popold5(I == 2, :), valParents5(I == 2));
            popold5(I == 2, :) = ui(I == 2, :);
            valParents5(I == 2,:) = valOffspring5(I == 2, :);
            goodCR5 = CR5(I == 2);
            goodF5 = F5(I == 2);
            if max(valParents5) > overallBestVal
                overallBestVal = max(valParents5);
            end
%             arrayGbestChange(5) = arrayGbestChange(5) + sum(prevalParents5- valParents5);
            if prevalParents5(I2(1)) == valParents5(I2(1))
                for r = 1 : 3
                    popold5(I2(r),:) =  mixPop(arrayFifth(I2(r)),:);
                    valParents5(I2(r)) = mixVal(arrayFifth(I2(r)));
                end
            end 
            mixPop(arrayFifth,:) = popold5;
            mixVal(arrayFifth) = valParents5; 
           end 
        elseif CSGS == 6
           arraySixth = permutation(1:mixPopSize);
          %% ============================mutation 6 =====================================%%%%
           if ~isempty(arraySixth)
            pop6 = mixPop(arraySixth,:); % the old population becomes the current population 
            valParents6 = mixVal(arraySixth);
            popsize6 = length(arraySixth);
            [~,I1]=sort(mixVal, 'descend');
            [~,I2]=sort(valParents6, 'ascend');
            for r = 1 : 3
                pop6(I2(r),:) = mixPop(I1(r),:);
                valParents6(I2(r)) = mixVal(I1(r));   
            end
            prevalParents6 = valParents6;
            if gen > 1 && ~isempty(goodCR6) && sum(goodF6) > 0 % If goodF and goodCR are empty, pause the update
                CRm6 = (1 - c) * CRm6 + c * mean(goodCR6);
                Fm6 = (1 - c) * Fm6 + c * sum(goodF6 .^ 2) / sum(goodF6); % Lehmer mean
            else
                CRm6 =  CR(paraIndex(6)); 
                Fm6 = F(paraIndex(6));
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F6, CR6] = randFCR(popsize6, CRm6, 0.1, Fm6, 0.1);
            r0 = [1 : popsize6];
            popAll = [pop6; archive.pop];
            [r1, r2] = gnR1R2(popsize6, size(popAll, 1), r0);
            % Find the p-best solutions
            [~, indBest] = sort(valParents6, 'descend');
            pNP = max(round(pj * popsize6), 5); % choose at least two best solutions
            randindex = ceil(rand(1, popsize6) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop6(indBest(randindex), :); % randomly choose one of the top 100p% solutions
            rot = (0:1:popsize6-1);
            ind = randperm(2);
            a1  = randperm(popsize6);             % shuffle locations of vectors
            rt = rem(rot+ind(1),popsize6);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),popsize6);
            a3  = a2(rt+1);
            pm1 = pop6(a1,:);             % shuffled population 1
            pm2 = pop6(a2,:);             % shuffled population 2
            pm3 = pop6(a3,:);             % shuffled population 3
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            %DE/rand/2
            vi =pm1 + F6(:, ones(1, D)) .* repmat(rand(popsize6,1),1,D) .* ((pm2 - pm3) + (pop6(r1,:) - popAll(r2,:)));
             
            mask = rand(popsize6, D) > CR6(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize6)'; cols = floor(rand(popsize6, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize6 D], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop6(mask);

            valOffspring6 = zeros(popsize6,1);
             for i = 1 : popsize6
              for j = 1 : D
                if ui(i,j) >= threshold
                   ui(i,j) = 1;
                else
                   ui(i,j) = 0;
                end   
                if size(find(ui(i,:)==0),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 1;
                end
                if size(find(ui(i,:)==1),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 0;
                end
              end  
            end
            parfor i = 1 : popsize6
               valOffspring6(i) = fun(dataset,ui(i,:),1);
            end
            FES = FES + popsize6;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents6, I] = max([valParents6, valOffspring6], [], 2);
            popold6 = pop6;
            archive = updateArchive(archive, popold6(I == 2, :), valParents6(I == 2));
            popold6(I == 2, :) = ui(I == 2, :);
            valParents6(I == 2,:) = valOffspring6(I == 2, :);
            goodCR6 = CR6(I == 2);
            goodF6 = F6(I == 2);
            if max(valParents6) > overallBestVal
                overallBestVal = max(valParents6);
            end
%             arrayGbestChange(6) = arrayGbestChange(6) + sum(prevalParents6- valParents6);
            if prevalParents6(I2(1)) == valParents6(I2(1))
                for r = 1 : 3
                    popold6(I2(r),:) =  mixPop(arraySixth(I2(r)),:);
                    valParents6(I2(r)) = mixVal(arraySixth(I2(r)));
                end
            end 
            mixPop(arraySixth,:) = popold6;
            mixVal(arraySixth) = valParents6; 
           end 
        elseif CSGS == 7
           arraySeventh = permutation(1:mixPopSize);
          %% ============================mutation 7 =====================================%%%%
           if ~isempty(arraySeventh)
            pop7 = mixPop(arraySeventh,:); % the old population becomes the current population 
            valParents7 = mixVal(arraySeventh);
            popsize7 = length(arraySeventh);
            [~,I1]=sort(mixVal, 'descend');
            [~,I2]=sort(valParents7, 'ascend');
            for r = 1 : 3
                pop7(I2(r),:) = mixPop(I1(r),:);
                valParents7(I2(r)) = mixVal(I1(r));   
            end
            prevalParents7 = valParents7;
            if gen > 1 && ~isempty(goodCR7) && sum(goodF7) > 0 % If goodF and goodCR are empty, pause the update
                CRm7 = (1 - c) * CRm7 + c * mean(goodCR7);
                Fm7 = (1 - c) * Fm7 + c * sum(goodF7 .^ 2) / sum(goodF7); % Lehmer mean
            else
                CRm7 =  CR(paraIndex(7)); 
                Fm7 = F(paraIndex(7));
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F7, CR7] = randFCR(popsize7, CRm7, 0.1, Fm7, 0.1);
            r0 = [1 : popsize7];
            popAll = [pop7; archive.pop];
            [r1, r2] = gnR1R2(popsize7, size(popAll, 1), r0);
            % Find the p-best solutions
            [~, indBest] = sort(valParents7, 'descend');
            pNP = max(round(pj * popsize7), 5); % choose at least two best solutions
            randindex = ceil(rand(1, popsize7) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop7(indBest(randindex), :); % randomly choose one of the top 100p% solutions
            rot = (0:1:popsize7-1);
            ind = randperm(2);
            a1  = randperm(popsize7);             % shuffle locations of vectors
            rt = rem(rot+ind(1),popsize7);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),popsize7);
            a3  = a2(rt+1);
            pm1 = pop7(a1,:);             % shuffled population 1
            pm2 = pop7(a2,:);             % shuffled population 2
            pm3 = pop7(a3,:);             % shuffled population 3
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            %DE/best/2
            vi = pbest + F7(:, ones(1, D)) .* ((pop7(r1, :) - popAll(r2, :)) + (pm3 - pm2));
            
            mask = rand(popsize7, D) > CR7(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize7)'; cols = floor(rand(popsize7, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize7 D], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop7(mask);

            valOffspring7 = zeros(popsize7,1);
            for i = 1 : popsize7
              for j = 1 : D
                if ui(i,j) >= threshold
                   ui(i,j) = 1;
                else
                   ui(i,j) = 0;
                end   
                if size(find(ui(i,:)==0),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 1;
                end
                if size(find(ui(i,:)==1),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 0;
                end
              end  
            end
            parfor i = 1 : popsize7
               valOffspring7(i) = fun(dataset,ui(i,:),1);
            end
            FES = FES + popsize7;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents7, I] = max([valParents7, valOffspring7], [], 2);
            popold7 = pop7;
            archive = updateArchive(archive, popold7(I == 2, :), valParents7(I == 2));
            popold7(I == 2, :) = ui(I == 2, :);
            valParents7(I == 2,:) = valOffspring7(I == 2, :);
            goodCR7 = CR7(I == 2);
            goodF7 = F7(I == 2);
            if max(valParents7) > overallBestVal
                overallBestVal = max(valParents7);
            end
%             arrayGbestChange(7) = arrayGbestChange(7) + sum(prevalParents7- valParents7);
            if prevalParents7(I2(1)) == valParents7(I2(1))
                for r = 1 : 3
                    popold7(I2(r),:) =  mixPop(arraySeventh(I2(r)),:);
                    valParents7(I2(r)) = mixVal(arraySeventh(I2(r)));
                end
            end 
            mixPop(arraySeventh,:) = popold7;
            mixVal(arraySeventh) = valParents7; 
           end 
        elseif CSGS == 8
           arrayEighth = permutation(1:mixPopSize);
          %% ============================mutation 8 =====================================%%%%
           if ~isempty(arrayEighth)
            pop8 = mixPop(arrayEighth,:); % the old population becomes the current population 
            valParents8 = mixVal(arrayEighth);
            popsize8 = length(arrayEighth);
            [~,I1]=sort(mixVal, 'descend');
            [~,I2]=sort(valParents8, 'ascend');
            for r = 1 : 3
                pop8(I2(r),:) = mixPop(I1(r),:);
                valParents8(I2(r)) = mixVal(I1(r));   
            end
            prevalParents8 = valParents8;
            if gen > 1 && ~isempty(goodCR8) && sum(goodF8) > 0 % If goodF and goodCR are empty, pause the update
                CRm8 = (1 - c) * CRm8 + c * mean(goodCR8);
                Fm8 = (1 - c) * Fm8 + c * sum(goodF8 .^ 2) / sum(goodF8); % Lehmer mean
            else
                CRm8 =  CR(paraIndex(8)); 
                Fm8 = F(paraIndex(8));
            end
            % Generate CR according to a normal distribution with mean CRm, and std 0.1
            % Generate F according to a cauchy distribution with location parameter Fm, and scale parameter 0.1
            [F8, CR8] = randFCR(popsize8, CRm8, 0.1, Fm8, 0.1);
            r0 = [1 : popsize8];
            popAll = [pop8; archive.pop];
            [r1, r2] = gnR1R2(popsize8, size(popAll, 1), r0);
            % Find the p-best solutions
            [~, indBest] = sort(valParents8, 'descend');
            pNP = max(round(pj * popsize8), 5); % choose at least two best solutions
            randindex = ceil(rand(1, popsize8) * pNP); % select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); % to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop8(indBest(randindex), :); % randomly choose one of the top 100p% solutions
            rot = (0:1:popsize8-1);
            ind = randperm(2);
            a1  = randperm(popsize8);             % shuffle locations of vectors
            rt = rem(rot+ind(1),popsize8);        % rotate indices by ind(1) positions
            a2  = a1(rt+1);                 % rotate vector locations
            rt = rem(rot+ind(2),popsize8);
            a3  = a2(rt+1);
            pm1 = pop8(a1,:);             % shuffled population 1
            pm2 = pop8(a2,:);             % shuffled population 2
            pm3 = pop8(a3,:);             % shuffled population 3
            % == == == == == == == == == == == == == == == Mutation == == == == == == == == == == == == ==
            %DE/best/3
            vi = pbest + F8(:, ones(1, D)) .* ((pop8(r1, :) - popAll(r2, :)) + (pm3 - pm2)+(pop8 - pm1));
            
            mask = rand(popsize8, D) > CR8(:, ones(1, D)); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : popsize8)'; cols = floor(rand(popsize8, 1) * D)+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([popsize8 D], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop8(mask);

            valOffspring8 = zeros(popsize8,1);
             for i = 1 : popsize8
              for j = 1 : D
                if ui(i,j) >= threshold
                   ui(i,j) = 1;
                else
                   ui(i,j) = 0;
                end   
                if size(find(ui(i,:)==0),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 1;
                end
                if size(find(ui(i,:)==1),2)==D
                   r3 = randperm(D);
                   r4 = r3(1);
                   ui(i,r4) = 0;
                end
              end  
            end
            parfor i = 1 : popsize8
               valOffspring8(i) = fun(dataset,ui(i,:),1);
            end
            FES = FES + popsize8;
            % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
            % I == 1: the parent is better; I == 2: the offspring is better
            [valParents8, I] = max([valParents8, valOffspring8], [], 2);
            popold8 = pop8;
            archive = updateArchive(archive, popold8(I == 2, :), valParents8(I == 2));
            popold8(I == 2, :) = ui(I == 2, :);
            valParents8(I == 2,:) = valOffspring8(I == 2, :);
            goodCR8 = CR8(I == 2);
            goodF8 = F8(I == 2);
            if max(valParents8) > overallBestVal
                overallBestVal = max(valParents8);
            end
%             arrayGbestChange(8) = arrayGbestChange(8) + sum(prevalParents8- valParents8);
            if prevalParents8(I2(1)) == valParents8(I2(1))
                for r = 1 : 3
                    popold8(I2(r),:) =  mixPop(arrayEighth(I2(r)),:);
                    valParents8(I2(r)) = mixVal(arrayEighth(I2(r)));
                end
            end 
            mixPop(arrayEighth,:) = popold8;
            mixVal(arrayEighth) = valParents8; 
           end 
        end
        [outcome,best] = max(mixVal);
        if mod(gen,100) == 0 %disp the results per 1w fes
%             dlmwrite('genChange.txt',max(mixVal),'\t')
            countp = countp + 1;
            AccProcessionCSGS(countp, n) =  max(mixVal);
            SizeProcessionCSGS(countp, n) = size(find(mixPop(best, : )==1), 2);
            disp(max(mixVal));
            
        end   
    end
    time = time + 1;     
  end
%     figure(1); plot(1:MaxGen,curve); xlabel('Number of Iterations');
%     ylabel('Fitness Value'); title('MP-SaCoDE'); grid on;
    [DE_gbestval,ibest] = max(mixVal);
    DE_gbest = mixPop(ibest, : );
    CountSize = 0;
    parfor i = 1 : D
        if mixPop(ibest,i) == 1
            CountSize = CountSize + 1;
            selectedFeature(n, i) = 1;
        end    
    end
%     DE_Size(l) = CountSize;
%     Size_rate(l) = CountSize/D;
    overallBestVal = max(mixVal);
%     results(l) = overallBestVal;
    disp(n);
    disp(CountSize);
%     disp(DE_Size(l));
%     disp(Size_rate(l));
    disp(overallBestVal);
%    selectedFeature(n, :) = DE_gbest;
    AccResults(n,1) = overallBestVal;
    AccResults(n,2) = CountSize;
%     Feature = selectedFeature;
%    if n == 1 
%         %Data 1
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\grammatical_facial_expression01.txt');
%         dataset_train=grammatical_facial_expression01;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\grammatical_facial_expression01.txt');
%         dataset_test=grammatical_facial_expression01;
%         
%    elseif n == 2  
%         %Data 2
%          load('C:\Users\qywxb\Desktop\MP_MSaCoDE\SemeionHandwrittenDigit.txt');
%         dataset_train=SemeionHandwrittenDigit;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\SemeionHandwrittenDigit.txt');
%         dataset_test=SemeionHandwrittenDigit;
%         %load('SemeionHandwrittenDigit.txt');
%         %dataset=SemeionHandwrittenDigit;
%    elseif n == 3  
%         %Data 3
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\isolet5.txt');
%         dataset_train=isolet5;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\isolet5.txt');
%         dataset_test=isolet5;       
%    elseif n == 4  
%         %Data 4
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\MultipleFeaturesDigit.txt');
%         dataset_train=MultipleFeaturesDigit;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\MultipleFeaturesDigit.txt');
%         dataset_test=MultipleFeaturesDigit;
%    elseif n == 5   
%         %Data 5
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\HAPTDataSet.txt');
%         dataset_train=HAPTDataSet;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\HAPTDataSet.txt');
%         dataset_test=HAPTDataSet;
%    elseif n == 6
%         %Data 6
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\har.txt');
%         dataset_train=har;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\har.txt');
%         dataset_test=har;
%    elseif n == 7  
%         %Data 7
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\UJIIndoorLoc.txt');
%         dataset_train=UJIIndoorLoc;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\UJIIndoorLoc.txt');
%         dataset_test=UJIIndoorLoc;
%    elseif n == 8  
%         %Data 8
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\MadelonValid.txt');
%         dataset_train=MadelonValid;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\MadelonValid.txt');
%         dataset_test=MadelonValid;
%    elseif n == 9 
%         %Data 9
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\OpticalRecognitionofHandwritten.txt');
%         dataset_train=OpticalRecognitionofHandwritten;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\OpticalRecognitionofHandwritten.txt');
%         dataset_test=OpticalRecognitionofHandwritten;
%    elseif n == 10  
%         %Data 10
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\ConnectionistBenchData.txt');
%         dataset_train=ConnectionistBenchData;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\ConnectionistBenchData.txt');
%         dataset_test=ConnectionistBenchData;
%    elseif n == 11  
%         %Data 11
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\wdbc.txt');
%         dataset_train=wdbc;
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\wdbc.txt');
%         dataset_test=wdbc;
%    elseif n == 12  
%         %Data 12
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\LungCancer.txt');
%         dataset_train=LungCancer; 
%         
%         load('C:\Users\qywxb\Desktop\MP_MSaCoDE\test\LungCancer.txt');
%         dataset_test=LungCancer;
%    end
%    fit = zeros(1,100);        
%    for i = 1 :100
%        fit(i) = fun(dataset_train,dataset_test,Feature(n,:),1);
%    end
%    disp('Test Acc');
%    disp(mean(fit));
%    AccResults(n, 3) = mean(fit);
save('AccuracyPre.mat','AccResults');
save('SelectFeaturesPre.mat','selectedFeature');
%     l = l + 1;
  save('AccProcessionPreCSGS.mat','AccProcessionCSGS');
  save('SizeProcessionPreCSGS.mat','SizeProcessionCSGS');
  end
end
save('AccProcessionCSGS.mat','AccProcessionCSGS');
save('SizeProcessionCSGS.mat','SizeProcessionCSGS');
save('Accuracy.mat','AccResults');
save('SelectFeatures.mat','selectedFeature');
% disp(results);
% disp(DE_Size);
% disp(Size_rate);
toc;
