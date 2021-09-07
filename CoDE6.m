function [mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR5,goodF5,CRm5,Fm5,goodCR6,goodF6,CRm6,Fm6,weight] = ...
       CoDE6(mixPop,mixVal,nfeat,mixPopMax,mixValMax,nfeatMax,fitness,fitnessMin,mixPopfitMin,mixValfitMin,nfeatfitMin,selectD,overallBestVal,arrayGbestChange,archive,FES,FESj,goodCR,goodF,CRm1,Fm1,goodCR5,goodF5,CRm5,Fm5,goodCR6,goodF6,CRm6,Fm6,paraIndex,D,D1,dataset,gen,array,weight)
       % Objective function
       fun=@jFitnessFunction1;
       arrayFirst = array;
       arrayFifth = array;
       arraySixth = array;
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
%               if valOffspring1(i) > mixValMax
%                  for k = 1 : 1
%                      valOffspring1(i) = valOffspring1(i) + fun(dataset,ui(i,:),1);
%                      FES = FES + 1;
%                  end
%                  valOffspring1(i) = valOffspring1(i)/2;
                 if valOffspring1(i) > mixValMax
                    mixPopMax = ui(i,:);
                    mixValMax = valOffspring1(i);
                    nfeatMax = nfeatOffspring1(i);
                 if FES >= 500000
                    fitnessOffspring1(i) = (0.1*(nfeatOffspring1(i)/D1))-0.9*valOffspring1(i);
                else
                    fitnessOffspring1(i) = 0.1 - 0.9*valOffspring1(i);
                end
                 end 
%               end
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
      %% ============================mutation 5 =====================================%%%%
       if ~isempty(arrayFifth)
        pop5 = mixPop(arrayFifth,:); % the old population becomes the current population 
        valParents5 = mixVal(arrayFifth);
        nfeat05 = nfeat(arrayFifth);
        fitness05 = fitness(arrayFifth);
        popsize5 = length(arrayFifth);
        [~,I1]=sort(mixVal, 'ascend');
        [~,I2]=sort(valParents5, 'descend');
        for r = 1 : 3
            pop5(I2(r),:) = mixPop(I1(r),:);
            valParents5(I2(r)) = mixVal(I1(r));
            nfeat05(I2(r)) = nfeat(I1(r));
            fitness05(I2(r)) = fitness(I1(r));
        end
        prefitness5 = fitness05;
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
        [~, indBest] = sort(fitness05, 'ascend');
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
        nfeatOffspring5 = zeros(popsize5,1);
        fitnessOffspring5 = zeros(popsize5,1);
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
        for i = 1 : popsize5
           valOffspring5(i) = fun(dataset,ui(i,:),1);
           nfeatOffspring5(i) = size(find(ui(i,:) == 1 ),2);
           if FES >= 500000
                  fitnessOffspring5(i) = (0.1*(nfeatOffspring5(i)/D1))-0.9*valOffspring5(i);
          else
                  fitnessOffspring5(i) = 0.1 - 0.9*valOffspring5(i);
          end
        end
        for i = 1 : popsize5 
           if nfeatOffspring5(i) <= selectD
%               if valOffspring5(i) > mixValMax
%                  for k = 1 : 1
%                      valOffspring5(i) = valOffspring5(i) + fun(dataset,ui(i,:),1);
%                      FES = FES + 1;
%                  end
%                  valOffspring5(i) = valOffspring5(i)/2;
                 if valOffspring5(i) > mixValMax
                    mixPopMax = ui(i,:);
                    mixValMax = valOffspring5(i);
                    nfeatMax = nfeatOffspring5(i);
                    if FES >= 500000
                       fitnessOffspring5(i) = (0.1*(nfeatOffspring5(i)/D1))-0.9*valOffspring5(i);
                    else
                       fitnessOffspring5(i) = 0.1 - 0.9*valOffspring5(i);
                    end 
                 end 
%               end
            end 
          if fitnessOffspring5(i) < fitnessMin 
              fitnessMin = fitnessOffspring5(i);
              mixPopfitMin = ui(i,:);
              mixValfitMin = valOffspring5(i);
              nfeatfitMin = nfeatOffspring5(i);
          end  
        end   
        FES = FES + popsize5;
        % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
        % I == 1: the parent is better; I == 2: the offspring is better
        [fitness05, I] = min([fitness05,  fitnessOffspring5], [], 2);
        popold5 = pop5;
        archive = updateArchive(archive, popold5(I == 2, :), fitness05(I == 2));
        t0 = find(I == 2);
        t = size(t0,1);
        if t > 0
            for i = 1 : t
                 A = find(ui(t0(i),:) == 1);
                 B = find(popold5(t0(i),:) == 1);
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
        popold5(I == 2, :) = ui(I == 2, :);
        fitness05(I == 2,:) = fitnessOffspring5(I == 2, :);
        nfeat05(I == 2, :) = nfeatOffspring5(I == 2, :);
        valParents5(I == 2, :) = valOffspring5(I == 2, :);
        goodCR5 = CR5(I == 2);
        goodF5 = F5(I == 2);
        if min(fitness05) < overallBestVal
            overallBestVal = min(fitness05);
        end
        arrayGbestChange(5) = arrayGbestChange(5) + sum(prefitness5- fitness05);
        for r = 1 : 3
            if prefitness5(I2(r)) == fitness05(I2(r))                
                popold5(I2(r),:) =  mixPop(arrayFifth(I2(r)),:);
                valParents5(I2(r)) = mixVal(arrayFifth(I2(r)));
                nfeat05(I2(r)) = nfeat(arrayFifth(I2(r)));
                fitness05(I2(r)) = fitness(arrayFifth(I2(r)));
            end
        end 
        mixPop(arrayFifth,:) = popold5;
        mixVal(arrayFifth) = valParents5; 
        nfeat(arrayFifth) = nfeat05;
        fitness(arrayFifth) = fitness05;
       end 
      %% ============================mutation 6 =====================================%%%%
       if ~isempty(arraySixth)
        pop6 = mixPop(arraySixth,:); % the old population becomes the current population 
        valParents6 = mixVal(arraySixth);
        nfeat06 = nfeat(arraySixth);
        fitness06 = fitness(arraySixth);
        popsize6 = length(arraySixth);
        [~,I1]=sort(fitness, 'ascend');
        [~,I2]=sort(fitness06, 'descend');
        for r = 1 : 3
            pop6(I2(r),:) = mixPop(I1(r),:);
            valParents6(I2(r)) = mixVal(I1(r));
            nfeat06(I2(r)) = nfeat(I1(r));
            fitness06(I2(r)) = fitness(I1(r));
        end
        prefitness6 = fitness06;
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
        [~, indBest] = sort(fitness06, 'ascend');
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
        nfeatOffspring6 = zeros(popsize6,1);
        fitnessOffspring6 = zeros(popsize6,1);
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
        for i = 1 : popsize6
           valOffspring6(i) = fun(dataset,ui(i,:),1);
           nfeatOffspring6(i) = size(find(ui(i,:) == 1 ),2);
           if FES >= 500000
                  fitnessOffspring6(i) = (0.1*(nfeatOffspring6(i)/D1))-0.9*valOffspring6(i);
          else
                  fitnessOffspring6(i) = 0.1 - 0.9*valOffspring6(i);
          end
        end
        for i = 1 : popsize6 
           if nfeatOffspring6(i) <= selectD
%               if valOffspring6(i) > mixValMax
%                  for k = 1 : 1
%                      valOffspring6(i) = valOffspring6(i) + fun(dataset,ui(i,:),1);
%                      FES = FES + 1;
%                  end
%                  valOffspring6(i) = valOffspring6(i)/2;
                 if valOffspring6(i) > mixValMax
                    mixPopMax = ui(i,:);
                    mixValMax = valOffspring6(i);
                    nfeatMax = nfeatOffspring6(i);
                    if FES >= 500000
                       fitnessOffspring6(i) = (0.1*(nfeatOffspring6(i)/D1))-0.9*valOffspring6(i);
                    else
                       fitnessOffspring6(i) = 0.1 - 0.9*valOffspring6(i);
                    end    
                 end 
%               end
           end
           if fitnessOffspring6(i) < fitnessMin 
              fitnessMin = fitnessOffspring6(i);
              mixPopfitMin = ui(i,:);
              mixValfitMin = valOffspring6(i);
              nfeatfitMin = nfeatOffspring6(i);
          end  
        end          
        FES = FES + popsize6;
        % == == == == == == == == == == == == == == == Selection == == == == == == == == == == == == ==
        % I == 1: the parent is better; I == 2: the offspring is better
        [fitness06, I] = min([fitness06,  fitnessOffspring6], [], 2);
        popold6 = pop6;
        archive = updateArchive(archive, popold6(I == 2, :), fitness06(I == 2));
        t0 = find(I == 2);
        t = size(t0,1);
        if t > 0
            for i = 1 : t
                 A = find(ui(t0(i),:) == 1);
                 B = find(popold6(t0(i),:) == 1);
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
        popold6(I == 2, :) = ui(I == 2, :);
        fitness06(I == 2,:) = fitnessOffspring6(I == 2, :);
        nfeat06(I == 2, :) = nfeatOffspring6(I == 2, :);
        valParents6(I == 2, :) = valOffspring6(I == 2, :);
        goodCR6 = CR6(I == 2);
        goodF6 = F6(I == 2);
        if min(fitness06) < overallBestVal
            overallBestVal = min(fitness06);
        end
        arrayGbestChange(3) = arrayGbestChange(3) + sum(prefitness6- fitness06);
        for r = 1 : 3
            if prefitness6(I2(r)) == fitness06(I2(r))
                popold6(I2(r),:) =  mixPop(arraySixth(I2(r)),:);
                valParents6(I2(r)) = mixVal(arraySixth(I2(r)));
                nfeat06(I2(r)) = nfeat(arraySixth(I2(r)));
                fitness06(I2(r)) = fitness(arraySixth(I2(r)));                    
            end
        end 
        mixPop(arraySixth,:) = popold6;
        mixVal(arraySixth) = valParents6;
        nfeat(arraySixth) = nfeat06;
        fitness(arraySixth) = fitness06;
       end 
   