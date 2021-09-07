%% Data Partition
% addpath(genpath('data'));
Problem = {'Alizadeh-2000-v1','Alizadeh-2000-v2','Bittner-2000','Garber-2001','West-2001'};

 %% MAIN LOOP
for j = 1:length(Problem)
    p_name = Problem{j};
    data = load(p_name);
    Data = [data.fea, data.gnd];
    filepathtest = 'C:\Users\qywxb\Desktop\SaWDE\test\';
    filepathtrain = 'C:\Users\qywxb\Desktop\SaWDE\train\';
    cv = cvpartition(Data(:,end), 'holdout', 0.3);
    idxs = cv.test;
    Test = Data(idxs,:);
    Train = Data(~idxs,:);
    save([filepathtest,p_name],'Test');
    save([filepathtrain,p_name],'Train');
end