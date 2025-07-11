%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run in silico Experiment 2 and to genereate the associated figures
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
parent = fileparts(pwd);
if SIMoption == 0
    [X_all] = collectSimData(fullfile(parent,nameDataFolder),'Experiment2_steady_control'); % simulation results control
    [X_all_LDH_KO] = collectSimData(fullfile(parent,nameDataFolder),'Experiment2_steady_LDH_KO'); % simulation results LDH KO
end

%% save results in the following folder:
FolderPathData = fullfile(parent,'mySimData');

splitData = 50; % number of simulation runs saved in one file

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define Phases
tend = [300 20 1]; % vector containing duration of each Phase in min
num_phase = length(tend);

% ATPase rates
ATPase_on = [0 1 0];
ATPase_act = linspace(0,20,81)-0.48; % in mM/min
ATPase_act(1:2) = [];
[tf,loc]=ismember([0.5,10,20],ATPase_act+0.48);

%% Run simulation
% load model info and input (x0, params, BX, K_BX and clamp_idx)
getModelInput


if SIMoption == 1

    % initialize for ODE loop
    varlist = [ 1 : length( modelInfo.SVarList ) ];
    options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-9,'AbsTol',1e-9);

    %% Control simulations
    idx = 1;
    num_split = 1;
    X = cell(num_phase,3);
    fprintf( 'Control simulations:\n' );
    for j = 1:length(ATPase_act)
        fprintf('step #   : %d \t', j);
        x0_phase = x0;
        tendold = 0;
        fprintf('...solving ode...   ');
        for i =1:num_phase
            t = [];
            y = [];
            params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) = ATPase_on(i)*ATPase_act(j)*1e-3;
            [t,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params,clamp_idx);
            % toc
            X(i,1,j)={t};
            X(i,2,j)={y};

            x0_phase = y(end,:);
            tendold = tendold+tend(i);

            % compute fluxes and Kapp
            J = [];
            Kapp = [];
            for m = 1:length(t)
                [junk,J(m,:),Kapp(m,:)] = dXdTMuscleMetabolism_OxPhos_FT( t(m), y(m,:).', BX, K_BX, [],[],params,clamp_idx );
            end
            X(i,3,j)={J};

        end
        fprintf( 'completed.\n' )

        idx = idx+1;

        if SAVEoption == 1 && idx-1 == splitData
            fprintf( '...saving simulation results...\n' )
            X_all = X(:,:,1+(num_split-1)*splitData:end);
            save(fullfile(FolderPathData,['Experiment2_steady_control_',num2str(num_split)]),'X_all')
            num_split = num_split+1;
            idx = 1;
        end

    end

    if SAVEoption == 1
        %% save simulation results
        fprintf( '...saving simulation results...\n' )
        X_all = X(:,:,1+(num_split-1)*splitData:end);
        save(fullfile(FolderPathData,['Experiment2_steady_control_',num2str(num_split)]),'X_all')
    end


    %% LDH KO
    params(getfield(modelInfo.ParID, 'x_LDH_fiber')) = 0;
    % define clamped state variables for each of the three simulation phases
    clamp_idx = {};
    clamp_idx(1) = {[modelInfo.SVarID.lactate_extracellular, modelInfo.SVarID.H_fiber]}; % default: clamped ex. lactate
    clamp_idx(2) = {[modelInfo.SVarID.lactate_extracellular]}; % default: clamped ex. lactate
    clamp_idx(3) = {[modelInfo.SVarID.lactate_extracellular]}; % default: clamped ex. lactate

    idx = 1;
    num_split = 1;
    X_LDH_KO = cell(num_phase,3);
    fprintf( 'LDH KO simulations:\n' );
    for j = 1:length(ATPase_act)
        fprintf('step #   : %d \t', j);
        x0_phase = x0;
        tendold = 0;
        fprintf('...solving ode...   ');
        for i =1:num_phase
            t = [];
            y = [];
            params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) = ATPase_on(i)*ATPase_act(j)*1e-3;
            [t,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params,clamp_idx{i});
            % toc
            X_LDH_KO(i,1,j)={t};
            X_LDH_KO(i,2,j)={y};

            x0_phase = y(end,:);
            tendold = tendold+tend(i);

            % compute fluxes
            J = [];
            Kapp = [];
            for m = 1:length(t)
                [junk,J(m,:),Kapp(m,:)] = dXdTMuscleMetabolism_OxPhos_FT( t(m), y(m,:).', BX, K_BX, [],[],params,clamp_idx{i} );
            end
            X_LDH_KO(i,3,j)={J};

        end

        fprintf( 'completed.\n' )

        idx = idx+1;

        if SAVEoption == 1 && idx-1 == splitData
            %% save simulation results
            fprintf( '...saving simulation results...\n' )
            X_all = X_LDH_KO(:,:,1+(num_split-1)*splitData:end);
            save(fullfile(FolderPathData,['Experiment2_steady_LDH_KO_',num2str(num_split)]),'X_all')
            num_split = num_split+1;
            idx = 1;
        end
    end

    if SAVEoption == 1
        %% save simulation results
        fprintf( '...saving simulation results...\n' )
        X_all = X_LDH_KO(:,:,1+(num_split-1)*splitData:end);
        save(fullfile(FolderPathData,['Experiment2_steady_LDH_KO_',num2str(num_split)]),'X_all')
    end

    X_all = X;
    X_all_LDH_KO = X_LDH_KO;

end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots
color = [0.31 0.58 0.9; 0.8500 0.3250 0.0980; 0 0 0.6; 0 0 0; 0.6 0.6 0.6];
LineW = [0.75,0.75];
style={'-','--','-.'};

% x-Axis limits
t1 = -1;
t2 = tend(2)+1;

%% collect steady states
J_ATPase_steady = [];

J_CK_steady = [];
J_OxP_steady = [];
J_LDH_steady = [];
J_GAPDH_steady = [];
J_G3PDH_steady = [];
J_Gly_steady = [];

J_CK_steady_LDH_KO = [];
J_Gly_steady_LDH_KO = [];
J_OxP_steady_LDH_KO = [];
J_LDH_steady_LDH_KO = [];
J_GAPDH_steady_LDH_KO = [];
J_G3PDH_steady_LDH_KO = [];


for k = 1: length(ATPase_act)
    J = X_all{2,3,k};
    J_LDH_KO = X_all_LDH_KO{2,3,k};
    t = X_all{2,1,k}-tend(1);
    idx = find(t <=10);
    idx = idx(end);
    x = X_all{2,2,k};

    % CONTROL
    J_CK_steady= [J_CK_steady;J(end,getfield(modelInfo.FluxID,'CK_fiber'))*1e3];
    J_ATPase_steady= [J_ATPase_steady;J(end,getfield(modelInfo.FluxID,'ATPASE_fiber'))*1e3];
    J_OxP_steady = [J_OxP_steady;J(end,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3];
    J_LDH_steady = [J_LDH_steady;J(end,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3];
    J_GAPDH_steady = [J_GAPDH_steady;J(end,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3];
    J_G3PDH_steady = [J_G3PDH_steady;J(end,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3];
    J_Gly_steady = [J_Gly_steady;J(end,getfield(modelInfo.FluxID,'PYK_fiber'))*1e3+J(end,getfield(modelInfo.FluxID,'PGK_fiber'))*1e3-J(end,getfield(modelInfo.FluxID,'PFKa_fiber'))*1e3];

    % LDH KO
    J_Gly_steady_LDH_KO = [J_Gly_steady_LDH_KO;J_LDH_KO(end,getfield(modelInfo.FluxID,'PYK_fiber'))*1e3+J_LDH_KO(end,getfield(modelInfo.FluxID,'PGK_fiber'))*1e3-J_LDH_KO(end,getfield(modelInfo.FluxID,'PFKa_fiber'))*1e3];
    J_CK_steady_LDH_KO = [J_CK_steady_LDH_KO;J_LDH_KO(end,getfield(modelInfo.FluxID,'CK_fiber'))*1e3];
    J_OxP_steady_LDH_KO = [J_OxP_steady_LDH_KO;J_LDH_KO(end,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3];
    J_LDH_steady_LDH_KO = [J_LDH_steady_LDH_KO;J_LDH_KO(end,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3];
    J_GAPDH_steady_LDH_KO = [J_GAPDH_steady_LDH_KO;J_LDH_KO(end,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3];
    J_G3PDH_steady_LDH_KO = [J_G3PDH_steady_LDH_KO;J_LDH_KO(end,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3];

end



%% Fig. 7 (ATP, NADH, pH and Lac dynamics at ATPase = 20 mM/min)

fig = figure(700);
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;

% Fig. 7A (ATP)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1100 0.5838 0.3628 0.3412];
for j = 1:num_phase
    x = X_all{j,2,end};
    t = X_all{j,1,end}-tend(1);
    x_KO = X_all_LDH_KO{j,2,end};
    t_KO = X_all_LDH_KO{j,1,end}-tend(1);
    plot(t,x(:,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle',style{1})
    plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle',style{1})
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('ATP (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-3 15])
arrows(t1, t2, tend, ax)

% rectangle
Pos = ax.Position;
Xlim=[t1 t2];
Ylim=[-3 15];
dim = [18.0 10 3.0 0.5];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom:
ax1 = axes('Position',[0.32 0.7 0.3628/3.5 0.3412/3.2],'Box','on');
plot([0 0],[-0.07 0.07],'--k')
for j = 2:3
    x = X_all{j,2,end};
    t = X_all{j,1,end}-tend(1);
    x_KO = X_all_LDH_KO{j,2,end};
    t_KO = X_all_LDH_KO{j,1,end}-tend(1);
    hold on
    plot(t,x(:,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle',style{1})
    plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle',style{1})
end
plot([0 0],[-500 5000],'--k')
plot([tend(2) tend(2)],[-500 5000],'--k')
xlim([18 21])
ylim([10 10.5])


% Fig. 7B (NADH)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.5703 0.5838 0.3628 0.3412];
for j = 1:num_phase
    x = X_all{j,2,end};
    t = X_all{j,1,end}-tend(1);
    x_KO = X_all_LDH_KO{j,2,end};
    t_KO = X_all_LDH_KO{j,1,end}-tend(1);
    plot(t,x(:,getfield(modelInfo.SVarID,'NADH_fiber'))*1e6,'color',color(1,:),'LineWidth',LineW(1),'LineStyle',style{1})
    plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'NADH_fiber'))*1e6,'color',color(2,:),'LineWidth',LineW(1),'LineStyle',style{1})
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('NADH ($\mathrm{\mu}$M)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-3 10])
arrows(t1, t2, tend, ax)
lgd = legend('Control','LDH KO','Interpreter','latex');


% Fig. 7D (lactate)
ax = subplot(N1,N2,4); box on; hold on;
ax.Position = [0.5703 0.11 0.3628 0.3412];
for j = 1:num_phase
    x = X_all{j,2,end};
    t = X_all{j,1,end}-tend(1);
    x_KO = X_all_LDH_KO{j,2,end};
    t_KO = X_all_LDH_KO{j,1,end}-tend(1);
    plot(t,x(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle',style{1})
    plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle',style{1})

end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('LAC (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-3 10])
arrows(t1, t2, tend, ax)

% Fig. 7C (pH)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.11 0.1100 0.3628 0.3412];
for j = 1:num_phase
    x = X_all{j,2,end};
    t = X_all{j,1,end}-tend(1);
    x_KO = X_all_LDH_KO{j,2,end};
    t_KO = X_all_LDH_KO{j,1,end}-tend(1);
    plot(t,-log10(x(:,getfield(modelInfo.SVarID,'H_fiber'))),'color',color(1,:),'LineWidth',LineW(1),'LineStyle',style{1})
    plot(t_KO,-log10(x_KO(:,getfield(modelInfo.SVarID,'H_fiber'))),'color',color(2,:),'LineWidth',LineW(1),'LineStyle',style{1})

end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('pH (-)','Interpreter', 'latex')
xlim([t1 t2])
ylim([6.4 7.2])
arrows(t1, t2, tend, ax)


%% Fig. 8

fig = figure(800);
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;

% Fig. 8B (Control: ATP synthesis)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1100 0.5838 0.3628 0.3412];
area(J_ATPase_steady',[(J_OxP_steady*16)./J_ATPase_steady.*100,(J_Gly_steady)./J_ATPase_steady.*100,(J_CK_steady)./J_CK_steady.*100])
colororder([0.31 0.58 0.9; 0.6 0.6 0.6; 0.8500 0.3250 0.0980;0 0 0.6])
area(J_ATPase_steady',-[(J_ATPase_steady)./J_ATPase_steady.*100])
xlabel('ATP demand (mM/min)','Interpreter', 'latex')
ylabel({'ATP synthesis (\%)'},'Interpreter', 'latex')
ylim([-100 100])
lgd = legend('$J^{OxPhos}_{\mathrm{ATP}}$','$J^{Gly.}_{\mathrm{ATP}}$','$J^{CK}_{\mathrm{ATP}}$','Interpreter','latex');
yticks(linspace(-100,100,5))
title('\bf{Control}','Interpreter','latex')
lgd = legend('$J^{OxPhos}_{\mathrm{ATP}}$','$J^{Gly.}_{\mathrm{ATP}}$','$J^{CK}_{\mathrm{ATP}}$','$J^{ATPase}_{\mathrm{ATP}}$','Interpreter','latex');

for k = 1:length(loc)
    fprintf('Control: The OxPhos contribution to ATP synthesis at ATPase = %.1f mM/min is: %.0f %% \n', J_ATPase_steady(loc(k)),16.*J_OxP_steady(loc(k))./J_ATPase_steady(loc(k)).*100)
end

% Fig. 8E (LDH KO: ATP synthesis)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.11 0.1100 0.3628 0.3412];
colororder([0.31 0.58 0.9; 0.6 0.6 0.6; 0.8500 0.3250 0.0980;0 0 0.6])
area(J_ATPase_steady',[(J_OxP_steady_LDH_KO*16)./J_ATPase_steady.*100,(J_Gly_steady_LDH_KO)./J_ATPase_steady.*100,(J_CK_steady_LDH_KO)./J_CK_steady.*100])
area(J_ATPase_steady',-[(J_ATPase_steady)./J_ATPase_steady.*100])
xlabel('ATP demand (mM/min)','Interpreter', 'latex')
ylabel({'ATP synthesis (\%)'},'Interpreter', 'latex')
ylim([-100 100])
lgd = legend('$J^{OxPhos}_{\mathrm{ATP}}$','$J^{Gly.}_{\mathrm{ATP}}$','$J^{CK}_{\mathrm{ATP}}$','$J^{ATPase}_{\mathrm{ATP}}$','Interpreter','latex');
yticks(linspace(-100,100,5))
title('\bf{LDH KO}','Interpreter','latex')


for k = 1:length(loc)
    fprintf('LDH KO: The OxPhos contribution to ATP synthesis at ATPase = %.1f mM/min is: %.0f %% \n', J_ATPase_steady(loc(k)), 16.*J_OxP_steady_LDH_KO(loc(k))./J_ATPase_steady(loc(k)).*100)
end


% CONTROL
all = [J_GAPDH_steady,-J_OxP_steady,-J_LDH_steady,J_G3PDH_steady];
all_sum = sum(max(all,0),2);
% contribution relative to total NAD reduction
contribution = [-J_OxP_steady./all_sum.*100, -J_LDH_steady./all_sum.*100,(J_GAPDH_steady)./all_sum.*100 ];

for k = 1:length(loc)
    fprintf('Control: The LDH contribution (flux: %.2f mM/min) to pyruvate uptake at ATPase = %.1f mM/min is: %.0f %% \n',J_LDH_steady(loc(k)), J_ATPase_steady(loc(k)), J_LDH_steady(loc(k))./all_sum(loc(k)).*100)
end

% LDH KO
all_LDH_KO = [J_GAPDH_steady_LDH_KO,-J_OxP_steady_LDH_KO,-J_LDH_steady_LDH_KO,J_G3PDH_steady_LDH_KO];
all_sum_LDH_KO = sum(max(all_LDH_KO,0),2);
% contribution relative to total NAD reduction
contribution_LDH_KO = [-J_OxP_steady_LDH_KO./all_sum_LDH_KO.*100, -J_LDH_steady_LDH_KO./all_sum_LDH_KO.*100,(J_GAPDH_steady_LDH_KO)./all_sum_LDH_KO.*100 ];

for k = 1:length(loc)
    fprintf('LDH KO: The LDH contribution to pyruvate uptake at ATPase = %.1f mM/min is: %.0f %% \n', J_ATPase_steady(loc(k)), J_LDH_steady_LDH_KO(loc(k))./all_sum_LDH_KO(loc(k)).*100)
end


% Fig. 8C (Control: NAD reduction)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.5703 0.5838 0.3628 0.3412];
patch([J_ATPase_steady;J_ATPase_steady(end);J_ATPase_steady(1)],[contribution(:,2);0;0],[0.8500 0.3250 0.0980]) % LDH
hold on
patch([J_ATPase_steady;J_ATPase_steady(end);J_ATPase_steady(1)],[-1.*(100+contribution(:,1));-100;-100],[0.31 0.58 0.9]) % OxPhos
hold on
patch([J_ATPase_steady;J_ATPase_steady(end);J_ATPase_steady(1)],[100-contribution(:,3);100;100],[0.6 0.6 0.6]) % GAPDH
xlabel('ATP demand (mM/min)','Interpreter', 'latex')
ylabel({'NAD reduction (\%)'},'Interpreter', 'latex')
title('Control','Interpreter', 'latex')
lgd = legend('$J^{LDH}_{\mathrm{ATP}}$','$J^{OxPhos}_{\mathrm{ATP}}$','$J^{GAPDH}_{\mathrm{ATP}}$','Interpreter','latex');


% Fig. 8F (LDH KO: NAD reduction)
ax = subplot(N1,N2,4); box on; hold on;
ax.Position = [0.5703 0.11 0.3628 0.3412];
patch([J_ATPase_steady;J_ATPase_steady(end);J_ATPase_steady(1)],[contribution_LDH_KO(:,2);0;0],[0.8500 0.3250 0.0980]) % LDH
hold on
patch([J_ATPase_steady;J_ATPase_steady(end);J_ATPase_steady(1)],[-1.*(100+contribution_LDH_KO(:,1));-100;-100], [0.31 0.58 0.9]) % OxPhos
hold on
patch([J_ATPase_steady;J_ATPase_steady(end);J_ATPase_steady(1)],[100-contribution_LDH_KO(:,3);100;100],[0.6 0.6 0.6]) % GAPDH
xlabel('ATP demand (mM/min)','Interpreter', 'latex')
ylabel({'NAD reduction (\%)'},'Interpreter', 'latex')
title('LDH KO','Interpreter', 'latex')
lgd = legend('$J^{LDH}_{\mathrm{ATP}}$','$J^{OxPhos}_{\mathrm{ATP}}$','$J^{GAPDH}_{\mathrm{ATP}}$','Interpreter','latex');



%% MAE pH
% control
j = 2; % phase 2
x = X_all{j,2,end};
torig = X_all{j,1,end}-tend(1);
pH_control = -log10(x(:,getfield(modelInfo.SVarID,'H_fiber')));

% LDH KO
j = 2; % phase 2
x = X_all_LDH_KO{j,2,end};
t = X_all_LDH_KO{j,1,end}-tend(1);
pH_LDH_KO = interp1(t,-log10(x(:,getfield(modelInfo.SVarID,'H_fiber'))),torig);

MAE_pH = sum(abs(pH_control - pH_LDH_KO))./length(torig);

fprintf('The MAE of the pH time course is %0.3f \n',MAE_pH)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Suppl. Fig. 5 (Reaction rates associated with NADH and Lac and Pyr dynamics)
[tf,loc]=ismember([10,20],ATPase_act+0.48);

fig = figure();
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;

% Suppl. Fig. 5A (Low intensity: fluxes associated with NADH)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1100 0.5838 0.3628 0.3412];
k =1;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    x_KO = X_all_LDH_KO{j,3,loc(k)};
    t_KO = X_all_LDH_KO{j,1,loc(k)}-tend(1);

    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1)*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle','-')

    plot(t_KO,x_KO(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1)*1.5,'LineStyle',':')
    plot(t_KO,x_KO(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle','-.')
    plot(t_KO,-x_KO(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle','--')
    plot(t_KO,-x_KO(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle','-')

end

xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$J^k_{\mathrm{NADH}}$ (mM/min)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-2 2])
title('\bf{Low intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% Suppl. Fig. 5C (Low intensity : LAC and PYR)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.11 0.1100 0.3628 0.3412];
k = 1;
for j = 1:num_phase
    x = X_all{j,2,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    x_KO = X_all_LDH_KO{j,2,loc(k)};
    t_KO = X_all_LDH_KO{j,1,loc(k)}-tend(1);

    plot(t,x(:,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1))
    plot(t,x(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(1,:),'LineWidth',1.5*LineW(1),'LineStyle',':')
    plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1))
    plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(2,:),'LineWidth',1.5*LineW(1),'LineStyle',':')

end

xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('LAC and PYR (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-0.3 2])
title('\bf{Low intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% Suppl. Fig. 5B (Medium intensity: fluxes associated with NADH)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.5703 0.5838 0.3628 0.3412];
k =2;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    x_KO = X_all_LDH_KO{j,3,loc(k)};
    t_KO = X_all_LDH_KO{j,1,loc(k)}-tend(1);

    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1)*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1),'LineStyle','-')


    plot(t_KO,x_KO(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1)*1.5,'LineStyle',':')
    plot(t_KO,x_KO(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle','-.')
    plot(t_KO,-x_KO(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle','--')
    plot(t_KO,-x_KO(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1),'LineStyle','-')

end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$J^k_{\mathrm{NADH}}$ (mM/min)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-7.5 7.5])
title('\bf{Medium intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)


% Suppl. Fig. 5D (medium intensity: LAC and PYR)
ax = subplot(N1,N2,4); box on; hold on;
ax.Position = [0.5703 0.11 0.3628 0.3412];
k = 2;
for j = 1:num_phase
    x = X_all{j,2,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    x_KO = X_all_LDH_KO{j,2,loc(k)};
    t_KO = X_all_LDH_KO{j,1,loc(k)}-tend(1);

    plot(t,x(:,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3,'color',color(1,:),'LineWidth',LineW(1))
    plot(t,x(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(1,:),'LineWidth',1.5*LineW(1),'LineStyle',':')
    plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3,'color',color(2,:),'LineWidth',LineW(1))
    plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(2,:),'LineWidth',1.5*LineW(1),'LineStyle',':')

end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('LAC and PYR (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-1.5 10])
title('\bf{Medium intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)


%% Suppl. Fig. 6
fig = figure();
fig.Position = [142.3333 56.3333 1040 571.3333];
N1 = 3;
N2 = 3;


% generate colormaps
selection = [1,7:8:length(ATPase_act)];
N = length(selection);
[tf,loc2]=ismembertol([ATPase_act(selection)],[0:0.001:20],1e-6);

cm = [0 0 0; 0.1 0.2 0.5; 0.5 0.8 1];                 % CONTROL: Colormap (RGB)
cmi = interp1([0; 3; 20], cm, (0:0.001:20));          % interpolated Colormap
colormatrix1 = cmi(loc2,:);

cm2 = [0 0 0; 0.6500 0.1250 0.0980;0.8500 0.3250 0.0980; 1 0.6250 0.1980]; % LDH KO: Colormap (RGB)
cmi2 = interp1([0;4; 8; 20], cm2, (0:0.001:20));          % interpolated Colormap
colormatrix2 = cmi2(loc2,:);


% Suppl. Fig. 6A (CONTROL: NADH)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1300 0.7160 0.2134 0.2090];
idx = 1;
for k = selection
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'NADH_fiber'))*1e6,'color',colormatrix1(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlim([t1 t2])
ylim([-0.5 4])
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('NADH ($\mathrm{\mu}$M)','Interpreter', 'latex')
title('\bf{Control}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% Suppl. Fig. 6B (CONTROL: ATP stoichiometry factor of OxPhos)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.1300+1*0.2808 0.7160 0.2134 0.2090];
idx = 1;
for k = selection
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);

        Ka_NADH = 0.1e-6;

        NADH_fiber = x(:,getfield(modelInfo.SVarID,'NADH_fiber'));
        s_NADH = (1./(1+(Ka_NADH./NADH_fiber).^4));
        s_ATP = (46+9.*s_NADH)./(11/3)+1;

        plot(t,s_ATP,'color',colormatrix1(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$\mathrm{s_{ATP}}$(-)','Interpreter', 'latex')
xlim([t1 t2])
ylim([15.77 16.01])
title('Control')
arrows(t1, t2, tend, ax)


% rectangle
Pos = ax.Position;
Xlim=[t1 t2];
Ylim=[15.77 16.01];
dim = [0 15.993 5.0 0.01];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.1300+1*0.2808+0.12 0.7160+0.07 0.2134/4 0.2090/4],'Box','on');
hold on
idx = 1;
for k = selection
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);

        Ka_NADH = 0.1e-6;

        NADH_fiber = x(:,getfield(modelInfo.SVarID,'NADH_fiber'));
        s_NADH = (1./(1+(Ka_NADH./NADH_fiber).^4));
        s_ATP = (46+9.*s_NADH)./(11/3)+1;

        plot(t,s_ATP,'color',colormatrix1(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlim([0 5])
ylim([15.99940 16])
yticks([15.99940 16])


% Suppl. Fig. 6C (CONTROL: NADH stoichiometry factor of OxPhos)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.1300+2*0.2808 0.7160 0.2134 0.2090];
idx = 1;
for k = selection
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);

        Ka_NADH = 0.1e-6;

        NADH_fiber = x(:,getfield(modelInfo.SVarID,'NADH_fiber'));
        s_NADH = (1./(1+(Ka_NADH./NADH_fiber).^4));
        s_ATP = (46+9.*s_NADH)./(11/3)+1;

        plot(t,s_NADH,'color',colormatrix1(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$\mathrm{s_{NADH}}$(-)','Interpreter', 'latex')
xlim([t1 t2])
ylim([0.9 1.005])
title('Control')
arrows(t1, t2, tend, ax)
colormap(ax, cmi); % use costume colormap
n = 3;
caxis((n-1)/2*[0,20]) % align colour axis properly
cb1 =colorbar(ax,'Position', [0.93 0.7160 0.015 0.2090]);


% rectangle
Pos = ax.Position;
Xlim=[t1 t2];
Ylim=[0.9 1.005];
dim = [0 0.997 5.0 0.004];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.1300+2*0.2808+0.12 0.7160+0.07 0.2134/4 0.2090/4],'Box','on');
hold on
idx = 1;
for k = selection
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);

        Ka_NADH = 0.1e-6;

        NADH_fiber = x(:,getfield(modelInfo.SVarID,'NADH_fiber'));
        s_NADH = (1./(1+(Ka_NADH./NADH_fiber).^4));
        s_ATP = (46+9.*s_NADH)./(11/3)+1;

        plot(t,s_NADH,'color',colormatrix1(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlim([0 5])
ylim([0.9997 1])
yticks([0.9997 1])


% Suppl. Fig. 6D (LDH KO: NADH)
ax = subplot(N1,N2,4); box on; hold on;
ax.Position = [0.1300 0.4163 0.2134 0.2090];
idx = 1;
for k = selection
    for j = 1:num_phase
        x_KO = X_all_LDH_KO{j,2,k};
        t_KO = X_all_LDH_KO{j,1,k}-tend(1);
        plot(t_KO,x_KO(:,getfield(modelInfo.SVarID,'NADH_fiber'))*1e6,'color',colormatrix2(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlim([t1 t2])
ylim([-0.5 4])
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('NADH ($\mathrm{\mu}$M)','Interpreter', 'latex')
title('\bf{LDH KO}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% Suppl. Fig. 6E (LDH KO: ATP stoichiometry factor of OxPhos)
ax = subplot(N1,N2,5); box on; hold on;
hold on
ax.Position = [0.1300+1*0.2808 0.4163 0.2134 0.2090];
idx = 1;
for k = selection
    for j = 1:num_phase
        x_KO = X_all_LDH_KO{j,2,k};
        t_KO = X_all_LDH_KO{j,1,k}-tend(1);

        Ka_NADH = 0.1e-6;

        NADH_fiber_KO = x_KO(:,getfield(modelInfo.SVarID,'NADH_fiber'));
        s_NADH_KO = (1./(1+(Ka_NADH./NADH_fiber_KO).^4));
        s_ATP_KO = (46+9.*s_NADH_KO)./(11/3)+1;

        plot(t_KO,s_ATP_KO,'color',colormatrix2(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$\mathrm{s_{ATP}}$(-)','Interpreter', 'latex')
xlim([t1 t2])
ylim([15.77 16.01])
title('LDH KO')
arrows(t1, t2, tend, ax)


% rectangle
Pos = ax.Position;
Xlim=[t1 t2];
Ylim=[15.77 16.01];
dim = [0 15.993 5.0 0.01];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.1300+1*0.2808+0.12 0.4163+0.07 0.2134/4 0.2090/4],'Box','on');
hold on
idx = 1;
for k = selection
    for j = 1:num_phase
        x_KO = X_all_LDH_KO{j,2,k};
        t_KO = X_all_LDH_KO{j,1,k}-tend(1);

        Ka_NADH = 0.1e-6;

        NADH_fiber_KO = x_KO(:,getfield(modelInfo.SVarID,'NADH_fiber'));
        s_NADH_KO = (1./(1+(Ka_NADH./NADH_fiber_KO).^4));
        s_ATP_KO = (46+9.*s_NADH_KO)./(11/3)+1;

        plot(t_KO,s_ATP_KO,'color',colormatrix2(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlim([0 5])
ylim([15.99940 16])
yticks([15.99940 16])

% Suppl. Fig. 6F (LDH KO: NADH stoichiometry factor of OxPhos)
ax2 = subplot(N1,N2,6); box on; hold on;
ax2.Position = [0.1300+2*0.2808 0.4163 0.2134 0.2090];
idx = 1;
for k = selection
    for j = 1:num_phase
        x_KO = X_all_LDH_KO{j,2,k};
        t_KO = X_all_LDH_KO{j,1,k}-tend(1);

        Ka_NADH = 0.1e-6;

        NADH_fiber_KO = x_KO(:,getfield(modelInfo.SVarID,'NADH_fiber'));
        s_NADH_KO = (1./(1+(Ka_NADH./NADH_fiber_KO).^4));

        plot(t_KO,s_NADH_KO,'color',colormatrix2(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$\mathrm{s_{NADH}}$(-)','Interpreter', 'latex')
xlim([t1 t2])
ylim([0.9 1.005])
title('LDH KO')
arrows(t1, t2, tend, ax2)
colormap(ax2, cmi2); % use costume colormap
n = 3;
caxis((n-1)/2*[0,20]) % align colour axis properly
cb2 =colorbar(ax2,'Position',[0.93 0.4163 0.015 0.2090]);
% cb.Position = [0.93 0.4163 0.015 0.2090];

% rectangle
Pos = ax2.Position;
Xlim=[t1 t2];
Ylim=[0.9 1.005];
dim = [0 0.997 5.0 0.004];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.1300+2*0.2808+0.12 0.4163+0.07 0.2134/4 0.2090/4],'Box','on');
hold on
idx = 1;
for k = selection
    for j = 1:num_phase
        x_KO = X_all_LDH_KO{j,2,k};
        t_KO = X_all_LDH_KO{j,1,k}-tend(1);

        Ka_NADH = 0.1e-6;

        NADH_fiber_KO = x_KO(:,getfield(modelInfo.SVarID,'NADH_fiber'));
        s_NADH_KO = (1./(1+(Ka_NADH./NADH_fiber_KO).^4));

        plot(t_KO,s_NADH_KO,'color',colormatrix2(idx,:),'LineWidth',LineW(2),'LineStyle',style{1})
    end
    idx = idx+1;
end
xlim([0 5])
ylim([0.9997 1])
yticks([0.9997 1])
