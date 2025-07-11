function [X_all, tracker, refVarInfo] = collectMPSAdata(Path, filename)

file = dir(fullfile(Path,[filename,'_*']));
numFiles = length({file.name});

if numFiles == 0
    error(['Unable to find file or directory ', fullfile(Path,filename),'*'])
end

X = {};
track = [];
for i = 1:numFiles
    load(fullfile(Path, [filename,'_',num2str(i),'.mat']))
    X = cat(3,X,X_all);
    track = [track;tracker];
end

X_all = X;
tracker = track;
refVarInfo = refVarInfo;
end
