function exportGLMtoTex(stats,features,analysisStr)

    if nargin < 3
        analysisStr = 'unknownAnalysis';
    end

    % prepare the data for export
    clear tex
    tex.data = [stats.beta(2:end) stats.p(2:end)];
    tex.tableRowLabels = features;
    tex.tableColLabels = {'beta','p'};
    tex.dataFormat = {'%.3f',1,'%.4f',1};
    tex.tableCaption = strcat('GLM',analysisStr);
    latex = latexTable(tex); % and create the obj

    fid=fopen(strcat(analysisStr,'GLM','.tex'),'w');
    [nrows,ncols] = size(latex);
    for row = 1:nrows
        fprintf(fid,'%s\n',latex{row,:});
    end
    fclose(fid);