'Original PSO'

%% Problem
Problem = {'Alizadeh-2000-v1','Alizadeh-2000-v2','Bittner-2000','Garber-2001','West-2001','Nutt-2003-v2',...
    'Pomeroy-2002-v1','Pomeroy-2002-v2','Shipp-2002-v1','Armstrong-2002-v1','Dyrskjot-2003','Liang-2005'};

%% Objective function
fun2 = @jFitnessFunction2;

for l = 1 : length(Problem)   
    tic;
    p_name = Problem{l};
    results.p_name = p_name;   
    dataset = load(['C:\Users\c\Desktop\SaWDE\train\',p_name]);
    dataset = dataset.Train;
    feat=dataset(:,1:end-1); label=dataset(:,end);
    D = size(feat,2);
    opts.N = 100;
    opts.T = 1000000;
    PSO = jParticleSwarmOptimization(feat,label,opts);
    results.trainacc = 1 - PSO.fitG; 
    results.selectedfeatures = length(PSO.ff);
    testdataset = load(['C:\Users\c\Desktop\SaWDE\test\',p_name]);
    testdataset = testdataset.Test;
    fit = fun2(dataset, testdataset, PSO.sf ,1);
    results.testacc = fit;
    toc;
    time = num2str(toc);
    disp(time);
    results.time = time;
    saveResults(results);
end


