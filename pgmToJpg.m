files = dir;
idx = ~cellfun(@isempty,strfind({files.name},'.pgm'));

files = {files(idx).name};

for i = 1:length(files)
    tmp = imread(files{i});
    fname = strtok(files{i},'.');
    fname = strcat(fname,'.jpg');
    imwrite(tmp,fname,'jpg');
end