%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run in silico Experiment 3 and to genereate the associated figures
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
    load(fullfile(nameDataFolder,'Experiment3_LDH_KO_ischaemia.mat')) 
end

%% save results in the following folder:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load model info and input (x0, params, BX, K_BX and clamp_idx)
getModelInput

%% define Phases
tend = [300 1 0.05]; % vector containing duration of each Phase in min
num_phase = length(tend);

ATPase_act = 20-0.48;
ATPase_on = [0 1 0];

Voxphos = params(getfield(modelInfo.ParID,'x_OxPhosO2_fiber'));
ox_on = [1 0 0];

VCO2_weg = params(getfield(modelInfo.ParID,'x_CO2weg_extracellular_to_capillary'));
VHCO3_weg = params(getfield(modelInfo.ParID,'x_HCO3weg_extracellular_to_capillary'));
washout_on = [1 0 0];

%%
if SIMoption == 1
    %% Model configurations
    num_conf = 3;
    conf = ["control","LDH KO","LDH KO + clamped redox"];

    % LDH KO for model config 2 and 3
    Vldh = params(getfield(modelInfo.ParID,'x_LDH_fiber'));
    LDH_on = [1 0 0];

    % define clamped state variables for each of the three simulation phases
    clamp_idx = {};

    % control
    j = 1;
    clamp_idx(1,j) = {[modelInfo.SVarID.lactate_extracellular]}; % default: clamped ex. lactate
    clamp_idx(2,j) = {[]};
    clamp_idx(3,j) = {[]};

    % LDH KO
    j = 2;
    clamp_idx(1,j) = {[modelInfo.SVarID.lactate_extracellular, modelInfo.SVarID.H_fiber]}; % default: clamped ex. lactate
    clamp_idx(2,j) = {[]};
    clamp_idx(3,j) = {[]};

    % LDH KO and clamped redox
    j = 3;
    clamp_idx(1,j) = {[modelInfo.SVarID.lactate_extracellular, modelInfo.SVarID.H_fiber]}; % default: clamped ex. lactate
    clamp_idx(2,j) = {[modelInfo.SVarID.NAD_fiber, modelInfo.SVarID.NADH_fiber]};
    clamp_idx(3,j) = {[modelInfo.SVarID.NAD_fiber, modelInfo.SVarID.NADH_fiber]};

    %% initialize for ODE loop
    params_init = params;
    %% Control simulations
    X_all = cell(num_phase,5,num_conf);


    varlist = [ 1 : length( modelInfo.SVarList ) ];
    options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-9,'AbsTol',1e-9);

    for j = 1:num_conf
        params = params_init;
        fprintf( '%s simulation \n', conf(j))
        x0_phase = x0;
        tendold = 0;
        params(getfield(modelInfo.ParID, 'x_LDH_fiber')) = LDH_on(j)*Vldh;
        for i =1:num_phase
            t = [];
            y = [];
            params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) = ATPase_on(i)*ATPase_act*1e-3;
            params(getfield(modelInfo.ParID, 'x_OxPhosO2_fiber')) = ox_on(i)*Voxphos;
            params(getfield(modelInfo.ParID, 'x_CO2weg_extracellular_to_capillary')) = washout_on(i)*VCO2_weg;
            params(getfield(modelInfo.ParID, 'x_HCO3weg_extracellular_to_capillary')) = washout_on(i)*VHCO3_weg;


            [t,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params,clamp_idx{i,j});
            % toc
            X_all(i,1,j)={t};
            X_all(i,2,j)={y};

            x0_phase = y(end,:);
            tendold = tendold+tend(i);

            %% compute fluxes and Kapp
            J = [];
            Kapp = [];
            for m = 1:length(t)
                [junk,J(m,:),Kapp(m,:)] = dXdTMuscleMetabolism_OxPhos_FT( t(m), y(m,:).', BX, K_BX, [],[],params,clamp_idx{i,j});
            end
            X_all(i,3,j)={J};


        end
    end
    %
    if SAVEoption == 1
        save(fullfile(FolderPathData,'Experiment3_LDH_KO_ischaemia'),'X_all')
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots
% x-Axis limits
t1 = -0.05;
t2 = tend(2)+0.05;

%% Fig. 9F (accumulation of FBP, G3P, PYR and LAC)

fig = figure(904);
fig. Position = [509 182.3333 332*2 314];
ax = axes('Position',[0.1300 0.2127 0.6750 0.4023]);

deltas=[];
conf = ["Control","LDH KO","LDH KO with clamped redox"];
for k = 1:size(X_all,3)
    x = X_all{2,2,k};
    deltas(1,k) = x(end,getfield(modelInfo.SVarID,'fructose16phos_fiber'))*1e3-x(1,getfield(modelInfo.SVarID,'fructose16phos_fiber'))*1e3;
    deltas(2,k) = x(end,getfield(modelInfo.SVarID,'glycerol3phos_fiber'))*1e3-x(1,getfield(modelInfo.SVarID,'glycerol3phos_fiber'))*1e3;
    deltas(3,k) = x(end,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3-x(1,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3;
    deltas(4,k) = x(end,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3-x(1,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3;

    fprintf('%s:\t \x0394 FBP = %.1f mM;\t \x0394 G3P = %.1f mM;\t \x0394 PYR = %.1f mM;\t \x0394 LAC = %.1f mM;\t\n', conf(k),deltas(1,k),deltas(2,k),deltas(3,k),deltas(4,k))
end
fprintf('\n')

X = categorical({'FBP','G3P','PYR','LAC'});
X = reordercats(X,{'FBP','G3P','PYR','LAC'});
b = bar(X,deltas','FaceColor','flat');
b(1).CData = [0.31 0.58 0.9]; % Control
b(2).CData = [0.8500 0.3250 0.0980]; % LDH KO
b(3).CData = [0.6 0.6 0.6]; % LDH KO with clamped redox

ylabel('$\Delta$[$\mathrm{L_{total}}$] (mM)','Interpreter','latex')
title('\bf{Metabolite accumlation during exercise}','Interpreter','latex')

%% Fig. 9D-E (ATP and NADH)
color = [0.31 0.58 0.9; 0.8500 0.3250 0.0980;0.6 0.6 0.6];
LineW = [0.75];
style = {'-','-','-','-.','-'};

% Fig. 9D ATP
fig = figure(902);
fig. Position = [509 182.3333 332 314];
ax = axes('Position',[0.1300 0.2127 0.6750 0.4023]);
box on; hold on;
for k = 1:size(X_all,3)
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);

        plot(t,x(:,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle',style{k})

        if j == 2
            fprintf('%s: ATP decreases by %.1f mM \n', conf(k),x(1,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3-x(end,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3)
        end
    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('ATP (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([5 12])
arrows(t1, t2, tend, ax)

% rectangle
Xlim=[t1 t2];
Ylim=[-3 15];
Pos = ax.Position;
dim = [0.91 9.9 0.1 1];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.50 0.35 0.3628/2 0.3412/3],'Box','on');
plot([0 0],[-0.07 0.07],'--k')
hold on
for k = 1:size(X_all,3)

    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);

        plot(t,x(:,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle',style{k})
    end
end
plot([tend(2) tend(2)],[-5000 5000],'--k')
xlim([0.91 1.01])
ylim([10 10.5])


fprintf('\n')

% Fig. 9E NADH
fig = figure(903);
fig. Position = [509 182.3333 332 314];
ax = axes('Position',[0.1300 0.2127 0.6750 0.4023]);
box on; hold on;
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('NADH ($\mathrm{\mu}$M)','Interpreter', 'latex')

for k = 1:size(X_all,3)
    for j = 1:num_phase
        x = X_all{j,2,k};
        t = X_all{j,1,k}-tend(1);

        plot(t,x(:,getfield(modelInfo.SVarID,'NADH_fiber'))*1e6,'color',color(k,:),'LineWidth',LineW,'LineStyle',style{k})

        if j == 2
            fprintf('%s: NADH increases %.1f -fold \n', conf(k),max(x(:,getfield(modelInfo.SVarID,'NADH_fiber')))./x(1,getfield(modelInfo.SVarID,'NADH_fiber')))
        end
    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-2 25])
arrows(t1, t2, tend, ax)

fprintf('\n')

%% Fig. 9A (flux heatmap)
fluxes = ones(12,3);
for k = 1:size(X_all,3)
    for j = 2
        J = X_all{j,3,k};
        t = X_all{j,1,k}-tend(1);
        fluxes(1,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'PGLM_fiber'))*1e3)./t(end);
        fluxes(2,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'PGI_fiber'))*1e3)./t(end);
        fluxes(3,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'PFKa_fiber'))*1e3)./t(end);
        fluxes(4,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'FBA_fiber'))*1e3)./t(end);
        fluxes(5,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3)./t(end);
        fluxes(6,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'TPI_fiber'))*1e3)./t(end);
        fluxes(7,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3)./t(end);
        fluxes(8,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'PGK_fiber'))*1e3)./t(end);
        fluxes(9,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'PGYM_fiber'))*1e3)./t(end);
        fluxes(10,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'ENO_fiber'))*1e3)./t(end);
        fluxes(11,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'PYK_fiber'))*1e3)./t(end);
        fluxes(12,k) = trapz(t,J(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3)./t(end);
    end
    fprintf('%s: the glyocolytic ATP synthesis rate at t = 1 min is %.1f mM/min \n', conf(k),-J(end,getfield(modelInfo.FluxID,'PFKa_fiber'))*1e3+J(end,getfield(modelInfo.FluxID,'PGK_fiber'))*1e3+J(end,getfield(modelInfo.FluxID,'PYK_fiber'))*1e3)

end
fprintf('\n')

fprintf('LDH KO: %.0f %% increase in PFK flux compared to the control system\n',fluxes(3,2)./fluxes(3,1)*100-100)
fprintf('LDH KO: %.0f %% increase in PYK flux compared to the control system\n',fluxes(11,2)./fluxes(11,1)*100-100)
fprintf('LDH KO: %.0f %% increase in G3PDH flux compared to the control system\n',fluxes(5,2)./fluxes(5,1)*100-100)

fprintf('LDH KO with clamped redox: %.0f %% increase in PFK flux compared to the control system\n',fluxes(3,3)./fluxes(3,1)*100-100)
fprintf('LDH KO with clamped redox: %.0f %% increase in PYK flux compared to the control system\n',fluxes(11,3)./fluxes(11,1)*100-100)
fprintf('LDH KO with clamped redox: %.0f %% increase in G3PDH flux compared to the control system\n \n',fluxes(5,3)./fluxes(5,1)*100-100)


fluxes_abs = abs(fluxes);

figure(901)
names = {'$PGLM$','$PGI$','$PFKa$','$FBA$','$G3PDH$','$TPI$','$GAPDH$','$PGK$','$PGYM$','$ENO$','$PYR$','$LDH$'};
names2 = {'control','LDH def','LDH def, clamped redox'};
h = heatmap(names2,names,fluxes_abs);
h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';

cm = [0.1 0.2 0.5; 0.5 0.8 1 ; 1 1 1; 1 0.6 0.1 ; 0.9 0.2 0; 0.7 0 0.1]; % complete LDH KO
cmi = interp1([-1;-0.6; -0.4; 0.1; 0.5; 1], cm, (-1:0.001:1));          % interpolated Colormap 20 mM/min

colormap(cmi); % use costume colormap

n = 3; % should be odd
caxis((n-1)/2*[0,16]) % align colour axis properly 20 mM/min complete LDH KO

colorbar;
