function rs = saveResults(results)
  fid = [];
    if (exist([pwd filesep 'results.csv'], 'file') == 0)
        fid = fopen([pwd filesep 'results.csv'], 'w');
        fprintf(fid, '%s, %s, %s, %s, %s\n', ...
            'Data Set','Training Accuracy', 'Test Accuracy', 'Selected Features', 'Cost Time');
    elseif (exist([pwd filesep 'results.csv'], 'file') == 2)
        fid = fopen([pwd filesep 'results.csv'], 'a');
    end
    fprintf(fid, '%s, ', results.p_name);
    fprintf(fid, '%f, %f, %f, %s\n', ...
          results.trainacc, results.testacc, results.selectedfeatures, results.time);
    fclose(fid);
end



