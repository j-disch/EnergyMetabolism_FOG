%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to compute the loss and the KS statistics for the MPSA and to
% genereate the associated Fig. 4 in 
% "Dynamic balance of myoplasmic energetics and redox state in a 
% fast-twitch oxidative glycolytic skeletal muscle fiber"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

%% Setting paths
path(path,'../ExperimentalData');
path(path,'../Model');
path(path,'../mySimData');
path(path,'../refSimData');
path(path,'../Simulations');
path(path,'../Util');

%% Date of MPSA
DateMPSA = '2025_03_28';

%% Select options

COMPUTEoption = 0; % compute loss function and KS statistic value?

SAVEoption = 0; % save loss function values and KS statistic?

LOADref = 1; % if COMPUTEoption = 0 load data from \refSimDate (1) or \mySimData (0)?

%% Load data if computeERROR == 0 and/or computeKS == 0

% define folder
if LOADref == 0
    nameDataFolder = 'mySimData';
else
    nameDataFolder = 'refSimData';
end

% load results of Kolmogorov-Smirnov test
if COMPUTEoption == 0
    load(fullfile(nameDataFolder,[DateMPSA,'_MPSA_KS_all']))
end

%% save results in the following folder:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

%%
%%%%%%%%%%%%%%%%%%%%%%%% Compute loss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
load(fullfile(nameDataFolder,[DateMPSA,'_MPSA_LHC_and_var_matrix.mat'])) % load LHC matrix
num_strat = size(LHC,1);
tracker_all = 1;

num_phase = 3; % simulation phases
ATPase = [10 20 40]; % ATPase rates during muscle activation

%% Compute loss
if COMPUTEoption == 1
    Err_all = cell(num_phase,1,num_strat+1,length(ATPase));
    for k = 1:length(ATPase)
        [X_all, tracker, refVarInfo] = collectMPSAdata(fullfile(parent,nameDataFolder),[DateMPSA,'_MPSA_ATPase_',num2str(ATPase(k))]); % load simulation results of MPSA
        tracker_all = tracker_all.*tracker;


        %% inizialize reference simulation
        num_ref_metab = length(fieldnames(refVarInfo));
        SIM_ref = cell(num_phase,1); % simulation with reference parameterisation
        torig = cell(num_phase,1);
        SIM_strat = cell(num_phase,1); % simulations with

        %%
        for idx=1:num_strat+1

            if tracker(idx) > 1000
                % case: simulation crashed before reaching target t
                % or case: simulation exceeded stopTime
                Err_all(:,1,idx,:) = {1000*ones(1,num_ref_metab)};
            else

                for i = 1:num_phase

                    x = X_all{i,2,idx};
                    t = X_all{i,1,idx};


                    if idx == 1 % reference simulation
                        Err_all(i,1,idx,k) = {0*ones(1,num_ref_metab)};
                        x(end,:) = [];
                        t(end,:) = [];
                        SIM_ref(i,1) = {x(t>=301,:)};
                        torig(i,1) = {t(t>=301,:)};

                    else

                        SIM_strat(i,1) = {interp1(t,x,torig{i,1})};

                        if i ==1 % resting phase (baseline)
                            Err_all(i,1,idx,k) = {sum((SIM_strat{i,1}-SIM_ref{i,1}).^2)};
                            Delta = SIM_strat{1,1}(end,:)-SIM_ref{1,1}(end,:);
                        else

                            Err_all(i,1,idx,k) = {sum((SIM_strat{i,1}-SIM_ref{i,1}-Delta).^2)};

                        end
                    end
                end
            end
        end
    end
    Err_all = cell2mat(Err_all);
    if SAVEoption == 1
        %% Save errors
        name = [DateMPSA,'_MPSA_err_all'];
        save(fullfile(FolderPathData,name),'Err_all','tracker_all','refVarInfo','num_ref_metab');
    end


%%
%%%%%%%%%%%%%%%%%%%%%%%% K-S statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute KS statistics
% if computeKS == 1
    %% matrix with parameter sets par_sets
    var_name = par_sets.Properties.VariableNames;
    num_v = size(var_name,2);
    %% Combine active and recovery phase to dynamic phase
    err = Err_all;
    err = [err(1,:,:,:);sum(err([2:3],:,:,:),1)];
    %% K-S statistics
    % replace error values (2000 and 3000) with max calculated error

    % err:
    % 1st dimention: 2 states (baseline, dynamic response)
    % 2nd dimention: 10 ref state variables
    % 3rd dimention: 3001 stratifications
    % 4th dimention: 3 ATPase rates

    idx = find(tracker_all<1000);

    M1 = max(err(:,:,idx,:),[],[1 3]); % max err for every reference state variable
    idx_err = find(tracker_all>=1000);
    fprintf( '# err: %d.\n', length(idx_err) )
    new_err = err;
    new_err(:,:,idx_err,:) = ones(num_phase-1,1,size(idx_err,1)).*M1; % replace error value >= 1000 with max err of each reference variable

    %% calculate threshold for K-S-statistic

    % threshold(1x10x1x3):
    % 1st dimention: -
    % 2nd dimention: (10) reference state variables
    % 3rd dimention: -
    % 4th dimention: (3) ATPase rates
    threshold = mean(new_err(2,:,:,:),[1 3]); % mean of each ref variable over dynamic phase and across all data sets
    threshold_rest = mean(new_err(1,:,:,1),[1 3]); % mean of each ref state variable over resting phase

    % store KS-statistic in
    results_all = ones(num_v,size(new_err,4)+1,num_ref_metab);

    for k = 1:num_ref_metab % Pi, PCr, pH, ATP, ADP, AMP, FBP, G3P, Pyr, Lac
        for i = 1:size(new_err,4)+1 % rest, 10, 20, 40 mM/min

            % split parameter sets into two groups (acceptable Sa with err< threshold and
            % unacceptable Su with err >= threshold)
            if i == 1
                % Rest
                split = new_err(1,k,:,1)<threshold_rest(1,k,:,1);
            else
                % dynamic response
                split = new_err(2,k,:,i-1)<threshold(:,k,:,i-1); % error of dynamic phase for the kth metabolic state variable smaller than mean err of dynamic response across all param sets
            end

            Sa = find(split == 1); % acceptable
            Su = find(split == 0); % unacceptable

            % Perform the Kolmogorov-Smirnov test
            for j = 1: num_v
                if par_sets{1,j} == 0
                    results_all(j,i,k) = 0;
                else

                    [h,p,ks2stat] = kstest2(par_sets{Sa,j},par_sets{Su,j});

                    results_all(j,i,k) = ks2stat;

                end
            end
        end
    end

    if SAVEoption == 1
        %% Save KS statistics
        name = [DateMPSA,'_MPSA_KS_all'];
        save(fullfile(FolderPathData,name),'var_name','results_all');
        matname = [DateMPSA,'_MPSA_LHC_and_var_matrix'];
        save(fullfile(FolderPathData,matname),'LHC','par_sets','thestate')
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fig. 4 (heatmap MPSA)

indices = [];
for k = 1:size(results_all,2) % rest, 10, 20, 40 mM/min
    fig = figure();
    fig.Position = [384.3333 216.3333 143.3333 222.0000];

    results = results_all(:,k,:);
    results = round(results,2);

    idx = find(sum(maxk(results,1,3)>=0.1,3)>=1); % only consider parameters with K-S values > 0.1
    indices = unique([indices;idx]);

    names = var_name(1,sort(indices));


    h = heatmap(names,{'Pi', 'PCr', 'pH', 'ATP', 'ADP', 'AMP', 'FBP', 'G3P', 'PYR', 'LAC'}',results(sort(indices),:)');
    h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
    h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';


    cm = [1 1 1 ; 1 1 1;  1 0.6 0.1 ; 0.9 0.2 0; 0.9 0.2 0];         % Colormap (RGB)
    cmi = interp1([0; 0.08; 0.5; 0.8; 1], cm, (0:0.001:1));          % interpolated Colormap
    colormap(cmi); % use costume colormap

    n = 3; % should be odd
    caxis((n-1)/2*[0,1]) % align colour axis properly

    colorbar
    activity = {'\bf{rest}','\bf{low intensity}','\bf{medium intensity}','\bf{high intensity}'};
    h.NodeChildren(3).Title.Interpreter = 'latex';
    title(activity{k})
end

