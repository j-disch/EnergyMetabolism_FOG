function [] = addsppar(file, Inputdata)
% function creates new dXdT file for MPSA with variable specified
% parameters based on the input dXdT file
% The new file name is the name of the input dXdT file with the suffix '_MPSA'

% file: dXdT file ('dXdTMuscleMetabolism_OxPhos_FT.m')
% Inputdata: name and unspecified parameter (spparams_all)

%% Read dXdT file
filetext = fileread(file);

newstr = filetext;

%% new line style
if any(filetext == 13)
    line = [13, 10];
else
    line = 10;
end

%% get specified parameter names
vars = table2array(Inputdata(:,1));

%%
check_sppar = strfind(newstr,'sppar');

%%
if isempty(check_sppar)

    idx_start = strfind(newstr,file(1:end-2));
    idx_start = idx_start(2)+length(file(1:end-2));

    idx_end = strfind(newstr,')');
    idx_end = idx_end(idx_start < idx_end);
    idx_end = idx_end(1);

    % add new mandatory input parameter
    newstr = strrep(newstr,newstr(idx_start:idx_end),['_MPSA' newstr(idx_start:idx_end-1) ',sppar)']);

    % add description for new mandatory input parameter
    chr = 'parameter vector for the free parameters';
    chr_new = [chr line '%   sppar   parameter vector for the specified parameters'];
    newstr = strrep(newstr,chr,chr_new);

  
    for i = 1:length(vars)

        v_idx = [];


        v_idx = strfind(newstr, append(vars(i),"=")); % find variable

        if isempty(v_idx)
            v_idx = strfind(newstr, append(vars(i)," ="));
        end

        s_idx = strfind(newstr(v_idx(1):v_idx(1)+100), ";"); % find semicolon after variable
        s_idx = v_idx(1)+s_idx(1)-1;

        newstr = strrep(newstr,newstr(v_idx(1):s_idx),char(append(vars(i),'=sppar(',num2str(i),');'))); % make specified parameter variable

    end

    %% Create, edit and save new dXdT file
    % create (or open) Outputfile
    Outputfile = [file(1:end-2),'_MPSA.m'];
    fid = fopen(Outputfile, 'w');
    % write content into Outputfile
    fprintf(fid, '%s', newstr);
    % close file
    fclose(fid);
    %%
    fprintf('The dXdT file for the MPSA (%s) was (re)generated based on %s.\n',Outputfile,file)
else
    fprintf('Exit: the input dXdT file (%s) already contains the mandatory input parameters sppar.\n',file)
end
