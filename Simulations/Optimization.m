%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run the optimization for the parameter identification 
% in "Dynamic balance of myoplasmic energetics and redox state in a 
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

%% save optimization results in following folder:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

%% Load Experimental Data
load('Kushmerick1985_Fig4.mat')
Data(1,1) = {pH_2Hz};
Data(1,2)= {pH_4Hz};
Data(1,3) = {pH_10Hz};
Data(2,1) = {PCr_2Hz};
Data(2,2) = {PCr_4Hz};
Data(2,3) = {PCr_10Hz};
Data(3,1)= {Pi_2Hz};
Data(3,2) = {Pi_4Hz};
Data(3,3) = {Pi_10Hz};

%% load model info
load dXdTMuscleMetabolism_OxPhos_FT.mat

%% load inital conditions
load('Optimization_x0_init.mat')

%% Set proton buffer sizes
n = 3; % number of compartments
BX(1:n) = 0; % optimized value

%% Set proton buffer binding constants
K_BX = ones(1,n)*1e-7;

%% Set parameters for simulations

params = [];

% Vmax in M/min
params(getfield(modelInfo.ParID, 'x_PGLM_fiber')) =  3*0.9;
params(getfield(modelInfo.ParID, 'x_PGI_fiber')) =   1.2;
params(getfield(modelInfo.ParID, 'x_PFKa_fiber')) =  0.13;
params(getfield(modelInfo.ParID, 'x_FBA_fiber')) =  0.2;
params(getfield(modelInfo.ParID, 'x_TPI_fiber')) =  24;
params(getfield(modelInfo.ParID, 'x_GAPDH_fiber')) =  3.3;
params(getfield(modelInfo.ParID, 'x_G3PDH_fiber')) =  1.3;
params(getfield(modelInfo.ParID, 'x_PGK_fiber')) =  2;
params(getfield(modelInfo.ParID, 'x_PGYM_fiber')) =  1.6;
params(getfield(modelInfo.ParID, 'x_ENO_fiber')) =  0.3;
params(getfield(modelInfo.ParID, 'x_PYK_fiber')) =  0.2;
params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) =  0; % model input (optimized)
params(getfield(modelInfo.ParID, 'x_CK_fiber')) =  1;
params(getfield(modelInfo.ParID, 'x_AK_fiber')) =  1.2;
params(getfield(modelInfo.ParID, 'x_LDH_fiber')) =  4;

params(getfield(modelInfo.ParID, 'x_OxPhosO2_fiber')) =  0; % optimized value
params(getfield(modelInfo.ParID, 'Kop_ADP_fiber')) =  0; % optimized value
params(getfield(modelInfo.ParID, 'Kop_Pi_fiber')) =  1e-3;
params(getfield(modelInfo.ParID, 'Kop_PYR_fiber')) =  0.15e-3;
params(getfield(modelInfo.ParID, 'nH_op_fiber')) =  0; % optimized value

params(getfield(modelInfo.ParID, 'x_MCT_fiber_to_extracellular')) =  0.0877;
params(getfield(modelInfo.ParID, 'Kmct_lac_fiber_to_extracellular')) =  0.016;
params(getfield(modelInfo.ParID, 'x_CO2Diff_fiber_to_extracellular')) = 25;

params(getfield(modelInfo.ParID, 'C_fiber_to_extracellular')) = 1;


params(getfield(modelInfo.ParID, 'C_extracellular_to_capillary')) = 1;
params(getfield(modelInfo.ParID, 'x_CO2weg_extracellular_to_capillary')) = 15;
params(getfield(modelInfo.ParID, 'x_HCO3weg_extracellular_to_capillary')) = 1.5;
params(getfield(modelInfo.ParID, 'x_LacWeg_extracellular_to_capillary')) = 0;

%% define Phases
tend = [300.375 1.8 2.825]; % vector containing duration of each Phase in min

%% Settings for ODE loop
settings.x0 = x0;
settings.params = params;
settings.K_BX = K_BX;
settings.BX = BX;
settings.modelInfo = modelInfo;
settings.clamp_idx = [modelInfo.SVarID.lactate_extracellular]; % default: clamped ex. lactate
settings.ATPase_on = [0 1 0];
settings.tend = tend;

%% Define optimization function
my_opt_fun = @(p)costfunction_optimization(p,settings,Data,[],[]);

%% Set up random initial conditions
% nH_OXPHOS, Vmax_OXPHOS, Kadp_OXPHOS, ATPase(2Hz), ATPase(4Hz), ATPase(10Hz), XB and Pi_init (TPP)
p_range = [2,10,1e-5,10,10,10,0.031,0; ...
    3,20,10e-5,60,70,80,0.084,4.0e-3];

numOfVariables = size(p_range,2);

%% Do the optimization
n_run    = 20; % number of optimization runs
pvals    = cell(1,n_run);
err      = cell(1,n_run);
exitFlag = cell(1,n_run);
out      = cell(1,n_run);
pop      = cell(1,n_run);
scores   = cell(1,n_run);

rng(0)
rand_seed = randi([0 100000],1,n_run);


parfor idx=1:n_run
    disp(['Running optimization iteration ', num2str(idx)]);
    rng(rand_seed(idx))
    [pvals{idx}, err{idx},exitFlag{idx}, out{idx}, pop{idx}, scores{idx}] = ...
        ga(my_opt_fun,numOfVariables,[],[],[],[],p_range(1,:),p_range(2,:));
end

%% Save the results
formatOut = 'yyyy_mm_dd';
matname = [datestr(now,formatOut),'_opt_res'];
save(fullfile(FolderPathData,matname),'pvals','err','exitFlag','out','pop', 'scores','rand_seed')

%% End of optimization
display( 'Done!' );