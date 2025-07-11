function [X_all] = collectSimData(Path, filename)

file = dir(fullfile(Path,[filename,'_*']));
numFiles = length({file.name});

if numFiles == 0
    error(['Unable to find file or directory ', fullfile(Path,filename),'*'])
end


X = {};
for i = 1:numFiles
    load(fullfile(Path, [filename,'_',num2str(i),'.mat']))
    X = cat(3,X,X_all);
end

X_all = X;
end
