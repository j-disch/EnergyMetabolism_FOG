%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run the simulations for Suppl. Fig. 7
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
    load(fullfile(nameDataFolder,'Appendix_UnderdamedSystem.mat')) 
end

%% save results in the following folder:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load model info
load dXdTMuscleMetabolism_OxPhos_FT.mat

%% load x0, params, BX, K_BX and clamp_idx
load 2025_03_28_opt_modelInput.mat

%% define Phases
tend = [1 10 2]; % duration of phase in min
num_phase = length(tend);

ATPase_on = [0 1 0];

ATPase_act = 10-0.48; % in mM/min

% three model configurations: control, 1.5*Vmax_OxPhos, 2*Vmax_G3PDH
n_config = 3;
OxPhos_Vmax = [1 1.5 1];
OxPhos_control = params(getfield(modelInfo.ParID,'x_OxPhosO2_fiber'));
G3PDH_Vmax = [1 1 2];
G3PDH_control = params(getfield(modelInfo.ParID,'x_G3PDH_fiber'));

%%
if SIMoption == 1
    %% initialize for ODE loop

    varlist = [ 1 : length( modelInfo.SVarList ) ];
    options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-12,'AbsTol',1e-12);

    X_all = cell(num_phase,5,length(ATPase_act));
    tic
    for j = 1:n_config
        x0_phase = x0;
        tendold = 0;
        for i =1:num_phase
            t = [];
            y = [];
            params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) = ATPase_on(i)*ATPase_act*1e-3;
            params(getfield(modelInfo.ParID,'x_OxPhosO2_fiber')) = OxPhos_Vmax(j)*OxPhos_control;
            params(getfield(modelInfo.ParID,'x_G3PDH_fiber')) = G3PDH_Vmax(j)*G3PDH_control;


            [t,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params,clamp_idx);

            X_all(i,1,j)={t};
            X_all(i,2,j)={y};

            x0_phase = y(end,:);
            tendold = tendold+tend(i);

            %% compute fluxes
            J = [];
            Kapp = [];
            for m = 1:length(t)
                [junk,J(m,:),Kapp(m,:)] = dXdTMuscleMetabolism_OxPhos_FT( t(m), y(m,:).', BX, K_BX, [],[],params,clamp_idx );
            end
            X_all(i,3,j)={J};
            X_all(i,4,j)={Kapp};

            % THERMODYNAMIC DATA
            T =310.15;
            RT = 8.314*T/1e3; % kJ  mol^{-1}
            X_all(i,5,j)={-RT.*log(Kapp)};

        end
    end
    %%
    if SAVEoption == 1
        save(fullfile(FolderPathData,'Appendix_UnderdamedSystem'),'X_all')
    end
    toc
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots
LineW = 0.75;
color = [0.31 0.58 0.9; 0.8500 0.3250 0.0980; 0.6 0.6 0.6;0 0 0];

% x-Axis limits
t1 = -1;
t2 = tend(2)+2;

%% Suppl. Fig. 7 (PCr, NADH and fluxes associated with NADH)

fig = figure();
fig.Position = [142.3333 56.3333 1040 571.3333];
N1 = 3;
N2 = 3;

% Suppl. Fig. 7A-C (PCr)
for k = 1:n_config
    ax = subplot(N1,N2,k); box on; hold on;
    ax.Position = [0.1300+(k-1)*0.2808 0.7160 0.2134 0.2090];
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'phosphocreatine_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW)
    end
    xlabel('$t$ (min)','Interpreter', 'latex')
    ylabel('PCr (mM)','Interpreter', 'latex')
    xlim([t1 t2])
    ylim([-5 40])
    arrows(t1, t2, tend, ax)
end

% Suppl. Fig. 7D-F (NADH)
for k = 1:n_config
    ax = subplot(N1,N2,3+k); box on; hold on;
    ax.Position = [0.1300+(k-1)*0.2808 0.4163 0.2134 0.2090];
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'NADH_fiber'))*1e6,'color',color(k,:),'LineWidth',LineW)
    end
    xlabel('$t$ (min)','Interpreter', 'latex')
    ylabel('NADH ($\mathrm{\mu}$M)','Interpreter', 'latex')
    xlim([t1 t2])
    ylim([0 3])
    arrows(t1, t2, tend, ax)
end


% Suppl. Fig. 7G-I (fluxes associated with NADH)
for k = 1:n_config
    ax = subplot(N1,N2,6+k); box on; hold on;
    ax.Position = [0.1300+(k-1)*0.2808 0.1167 0.2134 0.2090];
    for j = 1:num_phase
        x = X_all{j,3,k};
        t = X_all{j,1,k}-tend(1);
        plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
        plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
        plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
        plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
        plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW/3)
    end
    xlabel('$t$ (min)','Interpreter', 'latex')
    ylabel('$J^k_{\mathrm{NADH}}$ (mM/min)','Interpreter', 'latex')
    xlim([t1 t2])
    ylim([-2.5 2.5])
    arrows(t1, t2, tend, ax)
end