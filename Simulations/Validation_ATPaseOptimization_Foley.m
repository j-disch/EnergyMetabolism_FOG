%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run ATPase rate optimization for model validation
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

%% save results in the following folder:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

%% Load Experimental Data
load('Foley1991.mat')

Data(1,1) = {Foley1991_pH};
Data(2,1) = {Foley1991_PCr};

%% load model info and input (x0, params, BX, K_BX and clamp_idx)
getModelInput

%% Make sure all the variables except the membrane potentials cannot become negative
varlist = [ 1 : length( modelInfo.SVarList ) ];

%% define Phases
tend = [1 8 7]; % duration of phases in min

%% initialize for ODE loop
settings.x0 = x0;
settings.params = params;
settings.K_BX = K_BX;
settings.BX = BX;
settings.modelInfo = modelInfo;
settings.ATPase_on = [0 1 0];
settings.clamp_idx = clamp_idx; % default: clamped ex. lactate
settings.tend = tend;

%% Define optimization function 2Hz
my_opt_fun = @(p)costfunction_ATPaseOptimization_Foley(p,settings,Data,[],[]);

%% Set up random initial conditions

% ATPase in M/min
p_range = [0; ...
            60*1e-3];

numOfVariables = size(p_range,2);

%% Do the optimization
tic
thestate = rng;
rng(thestate)
rng default
opts = optimoptions(@ga, 'Display','iter','PlotFcn', {@gaplotbestf,@gaplotstopping});
[p, err,exitFlag,Output, population, scores] = ga(my_opt_fun,numOfVariables,[],[],[],[],p_range(1,:),p_range(2,:),[],opts);
toc
fprintf('The optimized ATPase rate during exercise is %.2f mM/min \n',p*1e3+0.48)

%% Save the results
formatOut = 'yyyy_mm_dd';
matname = [datestr(now,formatOut),'_opt_res_Foley'];
save(fullfile(FolderPathData,matname),'p','err','exitFlag','Output','population', 'scores','thestate')

%% End of optimization
display( 'Done!' );