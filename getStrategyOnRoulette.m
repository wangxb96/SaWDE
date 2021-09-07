function strategySeq = getStrategyOnRoulette()
    num = 15;
    numOfCoDE = 10;
    P = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]./numOfCoDE;
    m = length(P);
    Select = zeros(1,num); 
    r = rand(1,num); 
    for i=1:num
        sumP = 0;
        r2 = randperm(m);
        j = r2(1); 
        while sumP < r(i)
            sumP = sumP + P(mod(j-1,m)+1);
            j = j+1;
        end
        Select(i) = mod(j-2,m)+1;
    end

    tempSelect = zeros(1,m);
    for i=1:m
        tempSelect(i) = sum(Select==i);
    end
    [~,ind_CSA] = sort(tempSelect,'descend');
    strategySeq = ind_CSA(1);