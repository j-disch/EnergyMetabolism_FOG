%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate latin hypercube sample matrix and to run MPSA simulations
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

%% load model info and input (x0, params, BX, K_BX and clamp_idx)
getModelInput

%% save results in:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

splitData = 1500; % number of simulation runs saved in one file

%% Load tables with parameter values
load('MPSA_parameter.mat')

% parameter values
sp = table2array(spparams(:,2))'; % specified parameters

up = table2array(uspparams(:,2))'; % unspecified parameters

% (if necessary) update 'up' values to the current optimized unspecified parameter values
for i = length(up)
    name = uspparams{i,1};
    up(i) = params(getfield(modelInfo.ParID,name));
end

%% Generate dXdT file for MPSA to include specified parameter vector as mandatory input argument
file = 'dXdTMuscleMetabolism_OxPhos_FT.m'; % model
Varnames = spparams(:,1); % list of specified parameters to be included in the MPSA

addsppar(file, Varnames)

%% Determine pool sizes to be included in MPSA
pool_names = ["TPP", "TCr", "TAN", "Redox", "pH", "buffer"];
% total phosphate pool
TPP = 2*x0(:,modelInfo.SVarID.ATP_fiber) + ...
    1*x0(:,modelInfo.SVarID.ADP_fiber) + ...
    x0(:,modelInfo.SVarID.Pi_fiber) + ...
    x0(:,modelInfo.SVarID.phosphocreatine_fiber) + ...
    x0(:,modelInfo.SVarID.glycerol3phos_fiber) + ...
    x0(:,modelInfo.SVarID.glucose1phos_fiber) + ...
    x0(:,modelInfo.SVarID.glucose6phos_fiber) + ...
    x0(:,modelInfo.SVarID.fructose6phos_fiber) + ...
    2*x0(:,modelInfo.SVarID.fructose16phos_fiber) + ...
    x0(:,modelInfo.SVarID.glyceraldehydephos_fiber) + ...
    x0(:,modelInfo.SVarID.dihydroxyacetonephos_fiber)+ ...
    2*x0(:,modelInfo.SVarID.bpg_fiber) + ...
    x0(:,modelInfo.SVarID.pg3_fiber)+ ...
    x0(:,modelInfo.SVarID.pg2_fiber)+ ...
    x0(:,modelInfo.SVarID.pep_fiber);

% total creatine pool
TCr = x0(:,modelInfo.SVarID.phosphocreatine_fiber) + ...
    x0(:,modelInfo.SVarID.creatine_fiber);

% total adenine nucleotide pool
TAN = x0(:,modelInfo.SVarID.ATP_fiber) + ...
    x0(:,modelInfo.SVarID.ADP_fiber);
x0(:,modelInfo.SVarID.AMP_fiber);

% Nicotinamide adenine dinucleotide
Redox =  x0(:,modelInfo.SVarID.NAD_fiber) + ...
    x0(:,modelInfo.SVarID.NADH_fiber);

% pH
pH = -log10(x0(:,modelInfo.SVarID.H_fiber));

% intrinisc buffer size
buffer = BX(1);


% initial pool sizes
pools = [TPP, TCr, TAN, Redox, pH, buffer];


% initial concentrations of variables used to adjust the pool sizes
pools_var = [x0(:,modelInfo.SVarID.Pi_fiber), x0(:,modelInfo.SVarID.creatine_fiber),...
    x0(:,modelInfo.SVarID.ATP_fiber), x0(:,modelInfo.SVarID.NAD_fiber),...
    pH, buffer];

%
ATP_orig = x0(:,modelInfo.SVarID.ATP_fiber);
PCr_orig = x0(:,modelInfo.SVarID.phosphocreatine_fiber);

%% Number of variables and stratifications
h_sp = height(spparams);
h_up = height(uspparams);
h_pools = length(pools);
num_v = h_sp+h_up+h_pools;
disp( ['# variables: ', num2str(num_v)]);
num_strat = 3000; % number of stratification/parameter sets
disp( ['# stratifications: ', num2str(num_strat)]);

%% Latin hypercube sampling
rng default
thestate = rng;
LHC = lhsdesign(num_strat,num_v); % Latin hypercube sample matrix

%% Perturbation
pert = 0.01;

%% get parameter sets
% set up matrix
X = ones(num_strat,num_v).*[sp, up, pools]; % matrix containing all reference parameter values
par_names = [table2array(spparams(:,1))', table2array(uspparams(:,1))', pool_names]; % all parameter names

% Generate matrix with all parameter sets assuming a uniform distribution
X = X-pert.*X+LHC.*(X+pert.*X-(X-pert.*X));

% change pools by changing the initial concentration of Pi, Cr, ATP, NAD, pH and BX repectively
X_con = ones(num_strat,h_pools).*[pools_var]; % ref concentration of state variable used to adjust pool size
X_pools = ones(num_strat,h_pools).*[pools]; % ref pool size

% edit matrix
X(:,h_up+h_sp+1:end) = X_con-pert.*X_pools+LHC(:,h_up+h_sp+1:end).*(2*pert.*X_pools);


par_sets = [sp, up, pools_var; X]; % adding reference parameter values
sp = par_sets(:,1:h_sp); % all specified parameter sets

%% create table with all parameter sets
par_sets= array2table(par_sets,"VariableNames",par_names);

%% Save LHC matrix parameter sets
formatOut = 'yyyy_mm_dd';
date = datestr(now,formatOut);
matname = [date,'_MPSA_LHC_and_var_matrix'];
save(fullfile(FolderPathData,matname),'LHC','par_sets','thestate')

%% Make sure all the variables except the membrane potentials cannot become negative
varlist = [ 1 : length( modelInfo.SVarList ) ];

%% define Phases
tend = [305 5 5]; % vector containing duration of each Phase in min
num_phase = length(tend);

%% define ATPase rates
ATPase_on = [0 1 0];
ATPase_Hz = [10 20 40]-0.48; % in mM/min

%% define clamped state variables for each of the three simulation phases
clamp_idx = {};
clamp_idx(1) = {[modelInfo.SVarID.lactate_extracellular, modelInfo.SVarID.H_fiber]}; % default: clamped ex. lactate
clamp_idx(2) = {[modelInfo.SVarID.lactate_extracellular]}; % default: clamped ex. lactate
clamp_idx(3) = {[modelInfo.SVarID.lactate_extracellular]}; % default: clamped ex. lactate

%% select reference state variables
% save only simulated time courses of the selected reference state variables
refVar = ["Pi_fiber","phosphocreatine_fiber","H_fiber","ATP_fiber","ADP_fiber","AMP_fiber","fructose16phos_fiber","glycerol3phos_fiber","pyruvate_fiber","lactate_fiber"];
idx_refVar = [];
refVarInfo = struct;

for m= 1:length(refVar)
    idx_refVar = [idx_refVar,getfield(modelInfo.SVarID,refVar(m))];
    refVarInfo = setfield(refVarInfo,refVar(m),m);
end
%% Run simulations

for j = 1:length(ATPase_Hz)
    fprintf('ATPase rate = %d mM/min \n', ATPase_Hz(j)+0.48);
    X_all = cell(num_phase,2);
    tracker = [];
    num_split = 1;
    idx = 1;

    for k=1:num_strat+1
        fprintf('step #   : %d \t', k);
        
        tendold = 0;

        tracker(idx,1) = 1;

        t = [];

        try

            % buffer size
            BX(:)= par_sets.buffer(k);
            % Vmax in M/min
            params(getfield(modelInfo.ParID, 'x_PGLM_fiber')) =  par_sets.x_PGLM_fiber(k);
            params(getfield(modelInfo.ParID, 'x_PGI_fiber')) =   par_sets.x_PGI_fiber(k);
            params(getfield(modelInfo.ParID, 'x_PFKa_fiber')) =  par_sets.x_PFKa_fiber(k);
            params(getfield(modelInfo.ParID, 'x_FBA_fiber')) =  par_sets.x_FBA_fiber(k);
            params(getfield(modelInfo.ParID, 'x_TPI_fiber')) =  par_sets.x_TPI_fiber(k);
            params(getfield(modelInfo.ParID, 'x_GAPDH_fiber')) =  par_sets.x_GAPDH_fiber(k);
            params(getfield(modelInfo.ParID, 'x_G3PDH_fiber')) =  par_sets.x_G3PDH_fiber(k);
            params(getfield(modelInfo.ParID, 'x_PGK_fiber')) =  par_sets.x_PGK_fiber(k);
            params(getfield(modelInfo.ParID, 'x_PGYM_fiber')) =  par_sets.x_PGYM_fiber(k);
            params(getfield(modelInfo.ParID, 'x_ENO_fiber')) =  par_sets.x_ENO_fiber(k);
            params(getfield(modelInfo.ParID, 'x_PYK_fiber')) =  par_sets.x_PYK_fiber(k);
            params(getfield(modelInfo.ParID, 'x_CK_fiber')) =  par_sets.x_CK_fiber(k);
            params(getfield(modelInfo.ParID, 'x_AK_fiber')) =  par_sets.x_AK_fiber(k);
            params(getfield(modelInfo.ParID, 'x_LDH_fiber')) =  par_sets.x_LDH_fiber(k);

            params(getfield(modelInfo.ParID, 'x_OxPhosO2_fiber')) =  par_sets.x_OxPhosO2_fiber(k);
            params(getfield(modelInfo.ParID, 'Kop_ADP_fiber')) =  par_sets.Kop_ADP_fiber(k);
            params(getfield(modelInfo.ParID, 'Kop_PYR_fiber')) =  par_sets.Kop_PYR_fiber(k);
            params(getfield(modelInfo.ParID, 'nH_op_fiber')) = par_sets.nH_op_fiber(k);


            params(getfield(modelInfo.ParID, 'x_MCT_fiber_to_extracellular')) = par_sets.x_MCT_fiber_to_extracellular(k);
            params(getfield(modelInfo.ParID, 'Kmct_lac_fiber_to_extracellular')) = par_sets.Kmct_lac_fiber_to_extracellular(k);
            params(getfield(modelInfo.ParID, 'x_CO2Diff_fiber_to_extracellular')) = par_sets.x_CO2Diff_fiber_to_extracellular(k);


            params(getfield(modelInfo.ParID, 'x_CO2weg_extracellular_to_capillary')) = par_sets.x_CO2weg_extracellular_to_capillary(k);
            params(getfield(modelInfo.ParID, 'x_HCO3weg_extracellular_to_capillary')) = par_sets.x_HCO3weg_extracellular_to_capillary(k);


            % edit initial conditions to change total pool sizes
            x0(getfield(modelInfo.SVarID, 'creatine_fiber')) =  par_sets.TCr(k)-2*(ATP_orig-par_sets.TAN(k));
            x0(getfield(modelInfo.SVarID, 'NAD_fiber')) =  par_sets.Redox(k);
            x0(getfield(modelInfo.SVarID, 'ATP_fiber')) =  par_sets.TAN(k);
            x0(getfield(modelInfo.SVarID, 'Pi_fiber')) =  par_sets.TPP(k);
            x0(getfield(modelInfo.SVarID, 'phosphocreatine_fiber')) = PCr_orig+2*(ATP_orig-par_sets.TAN(k));

            x0(getfield(modelInfo.SVarID, 'H_fiber')) =  10^(-par_sets.pH(k));



            %% unspecified parameters
            sppar = sp(k,:)';

            %% initial conditions
            x0_phase = x0;

            %% solve ode
            fprintf('...solving ode...   ');
            for i =1:num_phase
                t = [];
                y = [];
                stopTime = now + 5*60 / 86400;  % 5*60 seconds - Error: stops sim after max. 5min and goes to next parameter set
                options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-9,'AbsTol',1e-9,'OutputFcn', @(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx,sppar, flag) myOutputFcn(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx,sppar, flag,stopTime));

                params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) =  ATPase_on(i)*ATPase_Hz(j).*1e-3;


                [t,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT_MPSA, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params, clamp_idx{i},sppar);


                if t(end) < tend(i)+tendold
                    tracker(idx,1) = 2000;
                    fprintf( 'error 2000: model crashed.\n' )
                    break
                end
                % store simulation results for reference state variables
                X_all(i,1,idx)={t(t>=300,:)};
                X_all(i,2,idx)={y(t>=300,idx_refVar)};
                x0_phase = y(end,:);

                tendold = tendold+tend(i);
            end
            
        catch

            tracker(idx,1) = 3000;
            fprintf( 'error 3000.\n' )
        end

        idx = idx+1;
        fprintf( 'completed.\n' )

        if idx-1 == splitData
            %% save simulation results
            fprintf( '...saving simulation results...\n' )
            nameData = [date,'_MPSA_ATPase_',num2str(ATPase_Hz(j)+0.48),'_',num2str(num_split)];
            save(fullfile(FolderPathData,nameData),'X_all','tracker', 'refVarInfo');
            num_split = num_split+1;
            X_all = cell(num_phase,2);
            tracker = [];
            idx = 1;
        end

      
    end
    %% save simulation results
    fprintf( '...saving simulation results... \n' )
    nameData = [date,'_MPSA_ATPase_',num2str(ATPase_Hz(j)+0.48),'_',num2str(num_split)];
    save(fullfile(FolderPathData,nameData),'X_all','tracker', 'refVarInfo');

end

%% End of simulation
display( 'Done!' );

function status = myOutputFcn(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx,sppar, flag, stopTime)
status = double(now > stopTime);
end
