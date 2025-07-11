%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run validation simulations and to genereate the associated Fig. 3
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

%% Date of ATPase optimization run
DateATPaseOpt = '2025_03_28';

%% Select options
SIMoption = 0; % run simulation?
SAVEoption = 0; % save simulation results?

LOADref = 1; % if SIMoption = 0 load data from \refSimDate (1) or \mySimData (0)?

%% Load data if SIMoption == 0

% define folder
if LOADref == 0
    nameDataFolder = 'mySimData';
else
    nameDataFolder = 'refSimData';
end

% if SIMoption == 0 load simulation results
if SIMoption == 0
    load(fullfile(nameDataFolder,'Validation.mat')) % simulation results
end

%% save results in the following folder:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load model info and input (x0, params, BX, K_BX and clamp_idx)
getModelInput

%% define Phases
tend = [1 8 7]; % duration of phase in min
num_phase = length(tend);

% ATPase rates
% reported dPCr/dt(t=0) based on monoexponential fit through PCr in
% µmol/g/min (Foley at al. 1991)
ATPase_lit = 11.51;
ATPase_lit_low = 11.51-1.47;
ATPase_lit_up = 11.51+1.47;

% optimized ATPase rate
load(fullfile(nameDataFolder,[DateATPaseOpt,'_opt_res_Foley']),'p')
ATPase_opt = p;

ATPase_act = [ATPase_lit_low/0.66, ATPase_lit_up/0.66, ATPase_lit/0.66, ATPase_opt*1e3+0.48]-0.48;
ATPase_on = [0 1 0];

%%
if SIMoption == 1
    %% initialize for ODE loop

    X_all = cell(num_phase,2,length(ATPase_act));


    % Make sure all the variables except the membrane potentials cannot become negative
    varlist = [ 1 : length( modelInfo.SVarList ) ];

    options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-9,'AbsTol',1e-9);



    for j = 1:length(ATPase_act)
        x0_phase = x0;
        tendold = 0;
        for i =1:num_phase
            t = [];
            y = [];
            params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) =  ATPase_on(i)*ATPase_act(j)*1e-3;

            [t,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params,clamp_idx);

            X_all(i,1,j)={t};
            X_all(i,2,j)={y};

            x0_phase = y(end,:);
            tendold = tendold+tend(i);

        end
    end
    if SAVEoption == 1
        save(fullfile(FolderPathData,'Validation'),'X_all')
    end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load experimental data
load('Foley1991.mat')

%% Plots
LineW1 = 0.95;
color_shade = [0.9 0.9 0.9];
color_mean = [0.6 0.6 0.6];
color_opt = [0 0.1 0.8];

% x-axis limits
t1 = -1;
t2 = 14;

%% collect simulation results
pH=cell(length(ATPase_act),1);
PCr=cell(length(ATPase_act),1);
t=cell(length(ATPase_act),1);

for i = 1:length(ATPase_act)
    x = [];
    t_phase = [];
    for j = 1:num_phase
        x = [x;X_all{j,2,i}];
        t_phase = [t_phase;X_all{j,1,i}];
        x(end,:)=[];
        t_phase(end,:)=[];
    end
    pH(i) = {x(:,modelInfo.SVarID.H_fiber)};
    PCr(i) = {x(:,modelInfo.SVarID.phosphocreatine_fiber)};
    t(i)={t_phase};
end

t_lower = t{1,1}(1:end-1)-tend(1);
t_upper = t{2,1}(1:end-1)-tend(1);
t_mean = t{3,1}-tend(1);
t_opt = t{4,1}-tend(1);

PCr_lower = PCr{1,1}(1:end-1)*1e3;
PCr_upper = PCr{2,1}(1:end-1)*1e3;
PCr_mean = PCr{3,1}*1e3;
PCr_opt = PCr{4,1}*1e3;

pH_lower = -log10(pH{1,1}(1:end-1));
pH_upper = -log10(pH{2,1}(1:end-1));
pH_mean = -log10(pH{3,1});
pH_opt = -log10(pH{4,1});

%% Fig. 3 (Validation)
fig = figure(300);
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;

% Fig. 3A (PCr)
ax = subplot(N1,N2,1);
ax.Position = [0.1100 0.5838 0.3628 0.3412];
hold on
box on
% reported ATPase rate
patch([t_upper; flip(t_lower)],[PCr_upper; flip(PCr_lower)],color_shade,'EdgeColor','none') % shaded area
plot(t_mean,PCr_mean,'color',color_mean,'LineStyle','-','LineWidth',0.5); % mean
% optimized ATPase rate
plot(t_opt,PCr_opt,'color',color_opt,'LineStyle','-','LineWidth',LineW1);
% Foley data
PCr_steady = PCr_opt(1);
plot(Foley1991_PCr(:,1),Foley1991_PCr(:,2)./100*PCr_steady(1),'k*','MarkerSize',5)
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('PCr (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([5 40])
arrows(t1, t2, tend, ax)
lgd = legend('...','calculated ATPase rate $\pm$ SEM','optimised ATPase rate' ,'experimental data','','','Interpreter', 'latex');
lgd.Position = [0.7139 0.11+0.3412-0.1045 0.2218 0.1045];

% Fig. 3B (pH)
ax = subplot(N1,N2,2);
ax.Position = [0.5703 0.5838 0.3628 0.3412];
hold on
box on
% reported ATPase rate
patch([t_upper; flip(t_lower)],[pH_upper; flip(pH_lower)],color_shade,'EdgeColor','none') % shaded area
plot(t_mean,pH_mean,'color',color_mean,'LineStyle','-','LineWidth',0.5); % mean
% optimized ATPase rate
plot(t_opt,pH_opt,'color',color_opt,'LineStyle','-','LineWidth',LineW1);
% Foley data
plot(Foley1991_pH(:,1),Foley1991_pH(:,2),'k*','MarkerSize',5)
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('pH (-)','Interpreter', 'latex')
xlim([t1 t2])
ylim([6.8 7.15])
arrows(t1, t2, tend, ax)

%% calculate MAE (error between prediction and in vivio data)

Foley1991_PCr(1,1) =0;
Foley1991_PCr(1,2) =100;
Foley1991_pH(1,1) = 0;

t_pH  = Foley1991_pH(:,1);
t_PCr = Foley1991_PCr(:,1);

PCr_opt_int = interp1(t_opt,PCr_opt,t_PCr);
pH_opt_int = interp1(t_opt,pH_opt,t_pH);
PCr_mean_int = interp1(t_mean,PCr_mean,t_PCr);
pH_mean_int = interp1(t_mean,pH_mean,t_pH);

PCr_data = Foley1991_PCr(:,2)./100*PCr_steady(1);
pH_data = Foley1991_pH(:,2);

% mean absolute error (MAE)
err_PCr_mean = sum(abs([PCr_mean_int- PCr_data]))./length(PCr_data);
err_pH_mean = sum(abs([pH_mean_int- pH_data]))./length(pH_data);
err_PCr_opt = sum(abs([PCr_opt_int- PCr_data]))./length(PCr_data);
err_pH_opt = sum(abs([pH_opt_int- pH_data]))./length(pH_data);

fprintf('The MAE for PCr using the calculated ATPase rate is %.2f mM \n',err_PCr_mean)
fprintf('The MAE for pH using the calculated ATPase rate is %.2f \n',err_pH_mean)
fprintf('The MAE for PCr using the optimized ATPase rate is %.2f mM \n',err_PCr_opt)
fprintf('The MAE for pH using the optimized ATPase rate is %.2f \n',err_pH_opt)
