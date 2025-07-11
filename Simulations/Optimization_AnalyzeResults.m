%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to analyze the optimization results for the parameter identification 
% and to genereate the associated figures
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

%% Date of optimization run
DateOptimization = '2025_03_28';

%% Select options

SIMoption = 0; % run simulation?

SAVEoption = 0; % save simulation results?

SAVEmodelInput = 0; % save optimized model input (parameterization and new baseline)?

LOADref = 1; % if SIMoption = 0 load data from \refSimDate (1) or \mySimData (0)?

%% Load data if SIMoption == 0

% define folder
if LOADref == 0
    nameDataFolder = 'mySimData';
else
    nameDataFolder = 'refSimData';
end


% if SIMoption == 0 load simulation results and optimized parameterization and baseline
if SIMoption == 0
    load(fullfile(nameDataFolder,'Optimization.mat')) % simulation results
    load(fullfile(nameDataFolder,[DateOptimization,'_opt_modelInput.mat'])) % optimized parameterization and baseline
end

%% save results in the following folder:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

%% Load optimization results
load(fullfile(nameDataFolder,[DateOptimization,'_opt_res.mat']))

%% get optimization results
p_all = [];
err_all = [];
for i = 1:size(pvals,2)
    p_all = [p_all;pvals{i}];
    err_all = [err_all;err{i}];
end

num_sets = size(p_all,1); % number of optimization runs

% optimized parameters
nH_OP = p_all(:,1);
Voxphos = p_all(:,2);
Kop_ADP = p_all(:,3);
ATPase_Hz = p_all(:,4:6);
BX_sol = p_all(:,7);
Pi_init = p_all(:,8);

% find optimiztaion run with min loss
best = find(min(err_all)==err_all);
disp( ['The optimized parameter set from run #',num2str(best), ' has the smalles loss function value'] );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load model info
load dXdTMuscleMetabolism_OxPhos_FT.mat

%% define Phases
tend = [300.375 1.8 4]; % vector containing duration of each Phase in min
num_phase = length(tend);

% ATPase on/off
ATPase_on = [0 1 0];

%% stimulation frequencies
Hz = [2 4 10];
num_Hz = length(Hz);

%% clamped state variables
clamp_idx = [modelInfo.SVarID.lactate_extracellular]; % default: clamped ex. lactate

%% Run simulation
if SIMoption == 1
    %% load inital conditions used for optimization
    load('Optimization_x0_init.mat')

    %% Set proton buffer sizes
    n = 3;
    BX(1:n) = 0; % optimized parameter

    %% Set proton buffer binding constants
    K_BX = ones(1,n)*1e-7;

    %% Set parameters to values used for optimization

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


    %% set optimized parameters
    params(getfield(modelInfo.ParID, 'nH_op_fiber')) = nH_OP(best); % set nH to opztimized value
    params(getfield(modelInfo.ParID, 'Kop_ADP_fiber')) =  Kop_ADP(best); % set Kop_ADP to optimized value
    BX(1:n) = BX_sol(best); % set BX to optimized value
    x0(getfield(modelInfo.SVarID, 'Pi_fiber')) =  Pi_init(best); %set TPP to optimized value

    
    %% Make sure all the variables except the membrane potentials cannot become negative
    varlist = [ 1 : length( modelInfo.SVarList ) ];

    %% initialize for ODE loop
    % store simulation results in X_all
    X_all = cell(num_phase,3,num_Hz);

    tendold = 0;

    value = [];
    val = 0;

    %% run simulation and plot results
    for k = 1:num_Hz
        x0_phase = x0;
        tendold = 0;
        t_all = [];
        x = [];
        for i =1:num_phase

            options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-12,'AbsTol',1e-12);

            % set ATPase rate to optimized value
            params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) =  ATPase_on(i)*ATPase_Hz(best,k)*1e-3;
            % set Vmax of OXPHOS model to optimized value
            params(getfield(modelInfo.ParID, 'x_OxPhosO2_fiber')) =  Voxphos(best)*1e-3/16;

            % solve ODE
            [t,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params,clamp_idx);


            tendold = tendold+tend(i);
            x0_phase = y(length(t),:);
            X_all(i,1,k)={t};
            X_all(i,2,k)={y};

        end
    end

    if SAVEoption == 1
        save(fullfile(FolderPathData,'Optimization'),'X_all')
    end


    %% Save optimized parameterization and new baseline
    x0 = X_all{1,2,1}(end,:); % new baseline with best optimized parameter set
    if SAVEmodelInput == 1
        name_modelInput = [DateOptimization,'_opt_modelInput.mat']; % name of file
        save(fullfile(FolderPathData,name_modelInput),'x0','params','BX','K_BX','clamp_idx')
    end

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute TPP from optimized init Pi values

TPP_woPi = 2*x0(:,modelInfo.SVarID.ATP_fiber) + ...
    1*x0(:,modelInfo.SVarID.ADP_fiber) + ...
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

TPP_all = TPP_woPi+p_all(:,8);

p_all(:,8) = TPP_all;

%% new order: Vmax, nH, Kadp, BX, TPP, ATPase(2Hz), ATPase(4Hz), ATPase(10Hz)
p_new = p_all;

p_new(:,1) = p_all(:,2)./16;
p_new(:,2) = p_all(:,1);
p_new(:,3) = p_all(:,3).*1e6;
p_new(:,4:5) = p_all(:,7:8).*1e3;
p_new(:,6:8) = p_all(:,4:6);



%% optimized parameter set
fprintf('\n')
fprintf('Optimized parameter set:\n')
fprintf('Vmax: %.2f mM/min \n',p_new(best,1))
fprintf('nH: %.1f (-) \n',p_new(best,2))
fprintf(['Kadp: %.3f',char(181),'M \n'],p_new(best,3))
fprintf('BX: %.1f mM \n',p_new(best,4))
fprintf('TPP: %.1f mM \n',p_new(best,5))
fprintf('ATPase 2Hz: %.1f mM/min \n',p_new(best,6))
fprintf('ATPase 4Hz: %.1f mM/min \n',p_new(best,7))
fprintf('ATPase 10Hz: %.1f mM/min \n',p_new(best,8))
fprintf('\n')

%% Coefficient of variation

S = std(p_new);
M = mean(p_new);

COV = S./M.*100; % Coefficient of variation in percent

fprintf('COV of Vmax: %.1f %% \n',COV(1))
fprintf('COV of nH: %.1f %% \n',COV(2))
fprintf('COV of Kadp: %.1f %% \n',COV(3))
fprintf('COV of BX: %.1f %% \n',COV(4))
fprintf('COV of TPP: %.1f %% \n',COV(5))
fprintf('COV of ATPase 2Hz: %.1f %% \n',COV(6))
fprintf('COV of ATPase 4Hz: %.1f %% \n',COV(7))
fprintf('COV of ATPase 10Hz: %.1f %% \n',COV(8))


%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load experimental Data (Kushmerick, 1985)
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

%% Plots
color_vector = [0.31 0.58 0.9; 0.8500 0.3250 0.0980; 0.6 0.6 0.6];

% x-Axis limits
t1 = -0.375;
t2 = 5;

%% Fig. 2 (parameter identification)
fig = figure(200);
set(fig,'Units', 'normalized', 'Position', [0.2, 0, 0.6, 0.9]);
N1 = 3;
N2 = 3;

for k = 1:size(Data,2)
    color = color_vector(k,:);


    Pi = [];
    PCr = [];
    pH = [];
    t = [];

    for j = 1:num_phase
        x = X_all{j,2,k};
        t_phase = X_all{j,1,k}-tend(1);
        Pi = [Pi; x(1:end-1,getfield(modelInfo.SVarID,'Pi_fiber'))];
        PCr = [PCr; x(1:end-1,getfield(modelInfo.SVarID,'phosphocreatine_fiber'))];
        pH = [pH; -log10(x(1:end-1,getfield(modelInfo.SVarID,'H_fiber')))];

        t = [t; t_phase(1:end-1)];
    end

    % Kushmerick1985
    pH_data = Data{1,k};
    PCr_data = Data{2,k};
    Pi_data = Data{3,k};

    % Model baseline
    Pi_rest = X_all{1,2,k}(end,getfield(modelInfo.SVarID,'Pi_fiber'))*1e3;
    PCr_rest = X_all{1,2,k}(end,getfield(modelInfo.SVarID,'phosphocreatine_fiber'))*1e3;
    pH_rest = -log10(X_all{1,2,k}(end,getfield(modelInfo.SVarID,'H_fiber')));
    Delta_pH = pH_rest-pH_data(1,2);


    % Fig. 2A-C (PCr)
    ax = subplot(N1,N2,k);
    box on; hold on;
    plot(PCr_data(:,1)-0.375,(PCr_data(:,2))./(PCr_data(1,2)).*(PCr_rest),'-k*','MarkerSize',4,'LineWidth',0.4) % experimental data
    plot(t,PCr*1e3,'color',color,'LineWidth',0.75); % simulation 
    xlim([t1 t2])
    ylim([0 40])
    xlabel('$t$ (min)','Interpreter', 'latex')
    ylabel('PCr (mM)','Interpreter', 'latex')
    title(sprintf('\\textbf{%d Hz}', Hz(k)),'Interpreter', 'latex')
    % start and end of muscle activation
    plot([0 0],[0 150],'--k')
    plot([1.8 1.8],[0 150],'--k')

    % Fig. 2D-F (Pi)
    ax = subplot(N1,N2,k+3);
    box on; hold on;
    plot(Pi_data(:,1)-0.375,(Pi_data(:,2))./(Pi_data(1,2)).*(Pi_rest+1),'-k*','MarkerS',4,'LineWidth',0.4) % experimental data
    plot(t,Pi*1e3+1,'color',color,'LineWidth',0.75); % simulation 
    xlabel('$t$ (min)','Interpreter', 'latex')
    ylabel('Pi (mM)','Interpreter', 'latex')
    xlim([t1 t2])
    ylim([0 45])
    % start and end of muscle activation
    plot([0 0],[0 150],'--k')
    plot([1.8 1.8],[0 150],'--k')

    % Fig. 2G-I (pH)
    ax = subplot(N1,N2,k+6);
    box on; hold on;
    plot(pH_data(:,1)-0.375,pH_data(:,2)+Delta_pH,'-k*','MarkerS',4,'LineWidth',0.4) % experimental data
    plot(t,pH,'color',color,'LineWidth',0.75); % simulation 
    xlabel('$t$ (min)','Interpreter', 'latex')
    ylabel('pH (-)','Interpreter', 'latex')
    xlim([t1 t2])
    ylim([6.4 7.2])
    % arrow
    Xlim=[-0.375 5];
    Ylim=[6.4 7.2];
    P1 = [tend(2)/2+tend(2)*0.31 Ylim(1)+(Ylim(2)-Ylim(1))*0.05]; % from point
    P2 = [tend(2) Ylim(1)+(Ylim(2)-Ylim(1))*0.05]; % to point
    Pos = ax.Position;
    X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
    X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
    Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
    Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
    an = annotation('textarrow', X_conv, Y_conv,'String',' act. ', 'HeadStyle','vback3','Interpreter', 'latex');
    an.HeadLength = 5;
    hold on
    P1 = [tend(2)/2-tend(2)*0.31 Ylim(1)+(Ylim(2)-Ylim(1))*0.05];
    P2 = [0 Ylim(1)+(Ylim(2)-Ylim(1))*0.05];
    X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
    X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
    Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
    Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
    an = annotation('textarrow', X_conv, Y_conv,'String',' ', 'HeadStyle','vback3');
    an.HeadLength = 5;
    % start and end of muscle activation
    plot([0 0],[6.4 7.2],'--k')
    plot([tend(2) tend(2)],[6.4 7.2],'--k')



    %% calculate MAE of Pi time course
    Pi_sim = interp1(t,Pi*1e3+1,Pi_data(:,1)-0.375);
    MAE_Pi = sum(abs((Pi_data(:,2))./(Pi_data(1,2)).*(Pi_rest+1) - Pi_sim))./length(Pi_data(:,1));
    fprintf('MAE of Pi time course at %d Hz: %.2f mM \n', Hz(k), MAE_Pi);

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Correlation coefficients (r_values) and p_values
[r_value p_value] = corrcoef(p_new)
Corr = r_value.*tril(ones(8));

%% Suppl. Fig 2 (Correlation matrix)
figure()

Corr_round = round(Corr,2);

names = {'$V_\mathrm{max}^\mathrm{OxPhos}$','$n_\mathrm{H}^\mathrm{OxPhos}$','$K_\mathrm{ADP}^\mathrm{OxPhos}$','[X]$_{\mathrm{T}}$','TPP','ATPase(2Hz)','ATPase(4Hz)','ATPase(10Hz)'};
h = heatmap(names,names,Corr_round);
h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';


cm = [0.1 0.2 0.5; 0.5 0.8 1 ; 1 1 1; 1 0.6 0.1 ; 0.9 0.2 0];   % Colormap (RGB)
cmi = interp1([-1;-0.6; 0; 0.7; 1], cm, (-1:0.001:1));          % interpolated Colormap
colormap(cmi); % use costume colormap

n = 3; % should be odd
caxis((n-1)/2*[-1,1]) % align colour axis properly
colorbar
