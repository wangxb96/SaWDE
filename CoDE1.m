function [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR2,goodF2,CRm2,Fm2,goodCR4,goodF4,CRm4,Fm4,weight] = ...
       CoDE1(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR2,goodF2,CRm2,Fm2,goodCR4,goodF4,CRm4,Fm4,paraIndex,D,D1,dataset,gen,array,weight)
       % Objective function
       fun=@jFitnessFunction1;
       arrayFirst = array;
       arraySecond = array;
       arrayFourth = array;
       F = [0.5 1 0.6 0.9 0.5 0.9 0.6 1];
       CR = [0.1 0.2 0.9 0.8 0.9 0.1 0.8 0.2];
       c = 1/10;
       pj = 0.1;
       threshold = 0.6;
   %%  %% ===========================mutation 1=====================================%%%%
      if ~isempty(arrayFirst) 
        pop1 = mixPop(arrayFirst,:); % the old population becomes the current population
        valParents1 = mixVal(arrayFirst);
        nfeat01 = nfeat(arrayFirst);
        fitness01 = fitness(arrayFirst);
        popsize = length(arrayFirst);
        [~,I1]=sort(fitness, 'ascend');%I1:the indices of xxpopsize,big->small
        [~,I2]=sort(fitness01, 'descend');%I2:the indices of arrayFirst,small->big
        for r = 1 : 3
            pop1(I2(r),:) = mixPop(I1(r),:);%put the overall best into arrayFirst's smallest part
            valParents1(I2(r)) = mixVal(I1(r));%put the overall best value into arrayFirst's smallest part
            nfeat01(I2(r)) = nfeat(I1(r));
            fitness01(I2(r)) = fitness(I1(r));
        end
        prefitness1 = fitness01;

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
        [~, indBest] = sort(fitness01, 'ascend');
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
        nfeatOffspring1 = zeros(popsize,1);
        fitnessOffspring1 = zeros(popsize,1);
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
        for i = 1 : popsize
           valOffspring1(i) = fun(dataset,ui(i,:),1);
           nfeatOffspring1(i) = size(find(ui(i,:) == 1 ),2);
           if FES >= 500000
                  fitnessOffspring1(i) = (0.1*(nfeatOffspring1(i)/D1))-0.9*valOffspring1(i);
          else
                  fitnessOffspring1(i) = 0.1 - 0.9*valOffspring1(i);
          end
        end
        for i = 1 : popsize 
           if nfeatOffspring1(i) <= selectD
              if valOffspring1(i) > mixValMax
%                  for k = 1 : 1
%                      valOffspring1(i) = valOffspring1(i) + fun(dataset,ui(i,:),1);
%                      FES = FES + 1;
%                  end
%                  valOffspring1(i) = valOffspring1(i)/2;
%                  if valOffspring1(i) > mixValMax
                    mixPopMax = ui(i,:);
                    mixValMax = valOffspring1(i);
                    nfeatMax = nfeatOffspring1(i);
                   if FES >= 500000
                      fitnessOffspring1(i) = (0.1*(nfeatOffspring1(i)/D1))-0.9*valOffspring1(i);
                  else
                      fitnessOffspring1(i) = 0.1 - 0.9*valOffspring1(i);
                  end
%                  end 
              end
           end 
           if fitnessOffspring1(i) < fitnessMin 
              fitnessMin = fitnessOffspring1(i);
              mixPopfitMin = ui(i,:);
              mixValfitMin = valOffspring1(i);
              nfeatfitMin = nfeatOffspring1(i);
          end               
        end
        FESj = FESj + popsize;
        FES = FES + popsize;
        % I == 1: the parent is better; I == 2: the offspring is better
        [fitness01, I] = min([fitness01, fitnessOffspring1], [], 2);
        popold1 = pop1;
        archive = updateArchive(archive, popold1(I == 2, :), fitness01(I == 2));
        t0 = find(I == 2);
        t = size(t0,1);
        if t > 0
            for i = 1 : t
                 A = find(ui(t0(i),:) == 1);
                 B = find(popold1(t0(i),:) == 1);
                 C = setdiff(A,B);
                 if size(C,2) == 0
                     break;
                 else
                     for j = 1 : size(C,2)
                         weight(1,C(j)) = weight(1,C(j)) + 1;
                     end
                 end
            end
        end
        popold1(I == 2, :) = ui(I == 2, :);
        fitness01(I == 2,:) = fitnessOffspring1(I == 2, :);
        nfeat01(I == 2,:) = nfeatOffspring1(I == 2,:);
        valParents1(I == 2,:) = valOffspring1(I == 2,:);
        goodCR = CRj(I == 2);
        goodF = Fj(I == 2);
        if min(fitness01) < overallBestVal
            overallBestVal = min(fitness01);
        end
       arrayGbestChange(1) = arrayGbestChange(1) + sum(prefitness1- fitness01);
       for r = 1 : 3
           if prefitness1(I2(r)) == fitness01(I2(r)) %if nothing changed at last,restore it.                
               popold1(I2(r),:) =  mixPop(arrayFirst(I2(r)),:);
               valParents1(I2(r)) = mixVal(arrayFirst(I2(r)));        
               nfeat01(I2(r)) = nfeat(arrayFirst(I2(r)));
               fitness01(I2(r)) = fitness(arrayFirst(I2(r)));
            end   
        end
        mixPop(arrayFirst,:) = popold1;
        mixVal(arrayFirst) = valParents1;
        nfeat(arrayFirst) = nfeat01;
        fitness(arrayFirst) = fitness01;
       end
        %% ===========================mutation 2=====================================%%%%
       if ~isempty(arraySecond)
        pop2 = mixPop(arraySecond,:); % the old population becomes the current population 
        valParents2 = mixVal(arraySecond);
        nfeat02 = nfeat(arraySecond);
        fitness02 = fitness(arraySecond);
        popsize2 = length(arraySecond);
        [~,I1]=sort(fitness, 'ascend');
        [~,I2]=sort(fitness02, 'descend');
        for r = 1 : 3
            pop2(I2(r),:) = mixPop(I1(r),:);
            valParents2(I2(r)) = mixVal(I1(r));
            nfeat02(I2(r)) = nfeat(I1(r));
            fitness02(I2(r)) = fitness(I1(r));
        end 
        prefitness2 = fitness02;

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
        nfeatOffspring2 = zeros(popsize2,1);
        fitnessOffspring2 = zeros(popsize2,1);
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
        for i = 1 : popsize2
           valOffspring2(i) = fun(dataset,ui(i,:),1);
           nfeatOffspring2(i) = size(find(ui(i,:) == 1 ),2);
           if FES >= 500000
                  fitnessOffspring2(i) = (0.1*(nfeatOffspring2(i)/D1))-0.9*valOffspring2(i);
          else
                  fitnessOffspring2(i) = 0.1 - 0.9*valOffspring2(i);
          end
        end
        for i = 1 : popsize2 
           if nfeatOffspring2(i) <= selectD
%               if valOffspring2(i) > mixValMax
%                  for k = 1 : 1
%                      valOffspring2(i) = valOffspring2(i) + fun(dataset,ui(i,:),1);
%                      FES = FES + 1;
%                  end
%                  valOffspring2(i) = valOffspring2(i)/2;
                 if valOffspring2(i) > mixValMax
                    mixPopMax = ui(i,:);
                    mixValMax = valOffspring2(i);
                    nfeatMax = nfeatOffspring2(i);
                    if FES >= 500000
                       fitnessOffspring2(i) = (0.1*(nfeatOffspring2(i)/D1))-0.9*valOffspring2(i);
                    else
                       fitnessOffspring2(i) = 0.1 - 0.9*valOffspring2(i);
                    end                    
                 end 
%               end
           end
           if fitnessOffspring2(i) < fitnessMin 
              fitnessMin = fitnessOffspring2(i);
              mixPopfitMin = ui(i,:);
              mixValfitMin = valOffspring2(i);
              nfeatfitMin = nfeatOffspring2(i);
          end            
        end   
        FES = FES + popsize2;
        % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
        % I == 1: the parent is better; I == 2: the offspring is better
        [fitness02, I] = min([fitness02, fitnessOffspring2], [], 2);
        popold2 = pop2;
        t0 = find(I == 2);
        t = size(t0,1);
        if t > 0
            for i = 1 : t
                 A = find(ui(t0(i),:) == 1);
                 B = find(popold2(t0(i),:) == 1);
                 C = setdiff(A,B);
                 if size(C,2) == 0
                     break;
                 else
                     for j = 1 : size(C,2)
                         weight(1,C(j)) = weight(1,C(j)) + 1;
                     end
                 end
            end
        end
        popold2(I == 2, :) = ui(I == 2, :);
        fitness02(I == 2,:) = fitnessOffspring2(I == 2, :);
        nfeat02(I == 2,:) = nfeatOffspring2(I == 2,:);
        valParents2(I == 2,:) = valOffspring2(I == 2,:);
        goodCR2 = CR2(I == 2);
        goodF2 = F2(I == 2);            
        if min(fitness02) < overallBestVal
            overallBestVal = min(fitness02);
        end 
        arrayGbestChange(2) = arrayGbestChange(2) + sum(prefitness2- fitness02);
        for r = 1 : 3
            if prefitness2(I2(r)) == fitness02(I2(r))
                popold2(I2(r),:) =  mixPop(arraySecond(I2(r)),:);
                valParents2(I2(r)) = mixVal(arraySecond(I2(r)));
                nfeat02(I2(r)) = nfeat(arraySecond(I2(r)));
                fitness02(I2(r)) = fitness(arraySecond(I2(r)));                    
            end
        end
        mixPop(arraySecond,:) = popold2;
        mixVa(arraySecond) = valParents2;
        nfeat(arraySecond) = nfeat02;
        fitness(arraySecond) = fitness02;            
       end
      %% ===========================mutation 4 =====================================%%%%
       if ~isempty(arrayFourth)
        pop4 = mixPop(arrayFourth,:); % the old population becomes the current population 
        valParents4 = mixVal(arrayFourth);
        nfeat04 = nfeat(arrayFourth);
        fitness04 = fitness(arrayFourth);
        popsize4 = length(arrayFourth);
        [~,I1]=sort(fitness, 'ascend');
        [~,I2]=sort(fitness04, 'descend');
        for r = 1 : 3
            pop4(I2(r),:) = mixPop(I1(r),:);
            valParents4(I2(r)) = mixVal(I1(r)); 
            nfeat04(I2(r)) = nfeat(I1(r));
            fitness04(I2(r)) = fitness(I1(r));
        end
        prefitness4 = fitness04;
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
        [~, indBest] = sort(fitness04, 'ascend');
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
        nfeatOffspring4 = zeros(popsize4,1);
        fitnessOffspring4 = zeros(popsize4,1);
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
        for i = 1 : popsize4
           valOffspring4(i) = fun(dataset,ui(i,:),1);
           nfeatOffspring4(i) = size(find(ui(i,:) == 1 ),2); 
           if FES >= 500000
                  fitnessOffspring4(i) = (0.1*(nfeatOffspring4(i)/D1))-0.9*valOffspring4(i);
          else
                  fitnessOffspring4(i) = 0.1 - 0.9*valOffspring4(i);
          end
        end
        for i = 1 : popsize4 
           if nfeatOffspring4(i) <= selectD
%               if valOffspring4(i) > mixValMax
%                  for k = 1 : 1
%                      valOffspring4(i) = valOffspring4(i) + fun(dataset,ui(i,:),1);
%                      FES = FES + 1;
%                  end
%                  valOffspring4(i) = valOffspring4(i)/2;
                 if valOffspring4(i) > mixValMax
                    mixPopMax = ui(i,:);
                    mixValMax = valOffspring4(i);
                    nfeatMax = nfeatOffspring4(i);
                    if FES >= 500000
                       fitnessOffspring4(i) = (0.1*(nfeatOffspring4(i)/D1))-0.9*valOffspring4(i);
                    else
                       fitnessOffspring4(i) = 0.1 - 0.9*valOffspring4(i);
                    end                                
                 end 
%               end
           end 
           if fitnessOffspring4(i) < fitnessMin 
              fitnessMin = fitnessOffspring4(i);
              mixPopfitMin = ui(i,:);
              mixValfitMin = valOffspring4(i);
              nfeatfitMin = nfeatOffspring4(i);
          end               
        end
        FES = FES + popsize4;
        % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
        % I == 1: the parent is better; I == 2: the offspring is better
        [fitness04, I] = min([fitness04, fitnessOffspring4], [], 2);
        popold4 = pop4;
        archive = updateArchive(archive, popold4(I == 2, :), fitness04(I == 2));
        t0 = find(I == 2);
        t = size(t0,1);
        if t > 0
            for i = 1 : t
                 A = find(ui(t0(i),:) == 1);
                 B = find(popold4(t0(i),:) == 1);
                 C = setdiff(A,B);
                 if size(C,2) == 0
                     break;
                 else
                     for j = 1 : size(C,2)
                         weight(1,C(j)) = weight(1,C(j)) + 1;
                     end
                 end
            end
        end
        popold4(I == 2, :) = ui(I == 2, :);
        fitness04(I == 2,:) = fitnessOffspring4(I == 2, :);
        nfeat04(I == 2, :) = nfeatOffspring4(I == 2, :);
        valParents4(I == 2, :) = valOffspring4(I == 2, :);
        goodCR4 = CR4(I == 2);
        goodF4 = F4(I == 2);
        if min(fitness04) < overallBestVal
            overallBestVal = min(fitness04);
        end
        arrayGbestChange(4) = arrayGbestChange(4) + sum(prefitness4- fitness04);
        for r = 1 : 3
            if prefitness4(I2(r)) == fitness04(I2(r))                
                popold4(I2(r),:) =  mixPop(arrayFourth(I2(r)),:);
                valParents4(I2(r)) = mixVal(arrayFourth(I2(r)));                   
                nfeat04(I2(r)) = nfeat(arrayFourth(I2(r)));
                fitness04(I2(r)) = fitness(arrayFourth(I2(r)));            
            end
        end 
        mixPop(arrayFourth,:) = popold4;
        mixVal(arrayFourth) = valParents4; 
        nfeat(arrayFourth) = nfeat04;
        fitness(arrayFourth) = fitness04;
       end   