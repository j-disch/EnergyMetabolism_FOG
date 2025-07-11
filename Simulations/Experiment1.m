%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to run in silico Experiment 1 and to genereate the associated figures
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
    load(fullfile(nameDataFolder,'Experiment1.mat')) % simulation results
end

%% save results in the following folder:
parent = fileparts(pwd);
FolderPathData = fullfile(parent,'mySimData');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define Phases
tend = [1 10 2]; % duration of phase in min
num_phase = length(tend);

% ATPase rates
ATPase_on = [0 1 0];
ATPase_act = linspace(0,30,61)-0.48; % in mM/min
ATPase_act(1) = [];
[tf,loc]=ismember([10,20,30],ATPase_act+0.48);

%% load model info and input (x0, params, BX, K_BX and clamp_idx)
getModelInput

%% Run simulation
if SIMoption == 1
    %% initialize for ODE loop

    varlist = [ 1 : length( modelInfo.SVarList ) ];
    options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-12,'AbsTol',1e-12);

    X_all = cell(num_phase,5,length(ATPase_act));
    tic
    for j = 1:length(ATPase_act)
        x0_phase = x0;
        tendold = 0;
        for i =1:num_phase
            t = [];
            y = [];
            params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) = ATPase_on(i)*ATPase_act(j)*1e-3;

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
        save(fullfile(FolderPathData,'Experiment1'),'X_all')
    end
    toc
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots
LineW = 0.75;
color = [0.31 0.58 0.9; 0.8500 0.3250 0.0980; 0.6 0.6 0.6;0 0 0];

% x-Axis limits
t1 = -1;
t2 = tend(2)+2;

%% Fig. 5 (PCr, ATP, pH and NADH dynamics)

fig = figure(500);
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;

% Fig. 5A (PCr and Pi)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1100 0.5838 0.3628 0.3412];
for k = 1:length(loc)
    for j = 1:num_phase
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'phosphocreatine_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW)
        plot(t,x(:,getfield(modelInfo.SVarID,'Pi_fiber'))*1e3,'color',color(k,:),'LineWidth',1.5*LineW,'LineStyle',':')
    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('PCr and Pi (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-5 45])
plot([0 0],[-500 5000],'--k')
plot([tend(2) tend(2)],[-500 5000],'--k')
arrows(t1, t2, tend, ax)

% Fig. 5B (ATP)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.5703 0.5838 0.3628 0.3412];
for k = 1:length(loc)
    for j = 1:num_phase
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW)
    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('ATP (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-1.875 15])
arrows(t1, t2, tend, ax)


% rectangle
Pos = ax.Position;
Xlim=[t1 t2];
Ylim=[-1.875 15];
dim = [8.0 9.5 3.0 1.0];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.5703+0.17 0.5838+0.09 0.3628/3 0.3412/3],'Box','on');
hold on
plot([0 0],[-0.07 0.07],'--k')
for k = 1:length(loc)
    for j = 2:3
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'ATP_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW)
    end

end
plot([tend(2) tend(2)],[-500 5000],'--k')
xlim([8 11])
ylim([9.5 10.5])


% Fig. 5C (pH)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.11 0.1100 0.3628 0.3412];
for k=1:length(loc)
    for j = 1:num_phase
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        plot(t,-log10(x(:,getfield(modelInfo.SVarID,'H_fiber'))),'color',color(k,:),'LineWidth',LineW)

    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('pH (-)','Interpreter', 'latex')
xlim([t1 t2])
ylim([6.0 7.2])
arrows(t1, t2, tend, ax)

% Fig. 5D (NADH)
ax = subplot(N1,N2,4); box on; hold on;
ax.Position = [0.5703 0.11 0.3628 0.3412];
for k = 1:length(loc)
    for j = 1:num_phase
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'NADH_fiber'))*1e6,'color',color(k,:),'LineWidth',LineW)
    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('NADH ($\mathrm{\mu}$M)','Interpreter', 'latex')
xlim([t1 t2])
ylim([0 3])
arrows(t1, t2, tend, ax)

%% Analysis of simulated dynamic response
for k = 1:length(loc)
    for j = 2
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        hold on
        PCr_drop = 100-(x(end,getfield(modelInfo.SVarID,'phosphocreatine_fiber'))./x(1,getfield(modelInfo.SVarID,'phosphocreatine_fiber'))*100);
        fprintf('PCr drops by %.2f %% at ATPase = %.0f mM/min \n',PCr_drop, ATPase_act(loc(k)))
        Pi_incr = x(end,getfield(modelInfo.SVarID,'Pi_fiber'))./x(1,getfield(modelInfo.SVarID,'Pi_fiber'));
        fprintf('Pi increases %.0f -fold at ATPase = %.0f mM/min \n',Pi_incr, ATPase_act(loc(k)))
        pH = -log10(x(end,getfield(modelInfo.SVarID,'H_fiber')));
        fprintf('pH decreases to %.1f at ATPase = %.0f mM/min \n',pH, ATPase_act(loc(k)))
    end
end

%% Compute transformed standard free energy of half reaction NAD_ox + 2e- = NAD_red
% ionic strength
I_new = 0.1;
I_tabulated_reactants = 0;
% temperature
T_old = 298.15;
T_new = 310.15;

RT = 8.314*T_new/1e3; % kJ  mol^{-1}


% ionic strength correction of dGf0
P = I_new^0.5;
Q = I_tabulated_reactants^0.5;
f = ( P / ( 1 + 1.6 * P ) - Q / ( 1 + 1.6 * Q ) );

RTalpha  = 9.20483e-3 * T_old - 1.28467e-5 * T_old^2 + 4.95199e-8 * T_old^3;
beta_G = RTalpha * f;

% standard free energy of formation of reference species (Alberty 1993)
DfG0_H = 0;
DfG0_NAD = 0;
DfG0_NADH  = 22.65;

% Stochiometry of species
v = [-1 -1 1];
% charge of species
z = [1 -1 -2];

% standard Gibbs reaction energy corrected to ionic strength
DeltaG0_redox = v(1)*DfG0_H + v(2)*DfG0_NAD + v(3)*DfG0_NADH - beta_G*(v(1)*z(1)^2+v(2)*z(2)^2+v(3)*z(3)^2);

% equilibrium constant of the reference reaction (only corrected to I)
K_ref = exp(-DeltaG0_redox/RT);

%% Binding polynomials
% dissociation constants of NAD and NADH
% Kh = Inf;
% Km = Inf;
% Kk = Inf;
% Binding polynomials for NAD- and NADH2-
P_NAD = 1;
P_NADH = 1;


%% collect steady state Gibbs energies of Redox half reaction and ATP hydrolysis
Q_ATPase_end = [];
DeltaG0_ATPase_end = [];
Q_redox_end = [];
DeltaG0_redox_end = [];

% Redox half reaction
for k = 1: length(ATPase_act)
    x = X_all{2,2,k};
    J = X_all{2,3,k};
    t = X_all{2,1,k};

    % equilibrium constant using the reference reaction corrected to I
    K_app = K_ref.*x(:,getfield(modelInfo.SVarID,'H_fiber')).*P_NADH./P_NAD; % apparent equilibrium constant
    DeltaG0_app = -RT.*log(K_app); % standard transformed free energy, DrG'0

    % RT*ln(Q)
    Q_redox = 2.5786*log(x(:,getfield(modelInfo.SVarID,'NADH_fiber'))./x(:,getfield(modelInfo.SVarID,'NAD_fiber')));

    % collect steady state values
    Q_redox_end = [Q_redox_end; Q_redox(end)];
    DeltaG0_redox_end = [DeltaG0_redox_end; DeltaG0_app(end)];
end
fprintf('The range of the transformed gibbs free energy of the half reaction is: %.2f kJ/mol \n',range(DeltaG0_redox_end+Q_redox_end))

% ATP hydrolysis
for k = 1: length(ATPase_act)
    x = X_all{2,2,k};
    J = X_all{2,3,k};
    t = X_all{2,1,k};
    DeltaG0 = X_all{2,5,k};

    DeltaGo_ATPase = DeltaG0(:,getfield(modelInfo.FluxID,'ATPASE_fiber')); % standard transformed free energy, DrG'0

    % RT*ln(Q)
    Q_ATPase = 2.5786*log((x(:,getfield(modelInfo.SVarID,'ADP_fiber')).*x(:,getfield(modelInfo.SVarID,'Pi_fiber')))./x(:,getfield(modelInfo.SVarID,'ATP_fiber')));

    %collect steady state values
    Q_ATPase_end = [Q_ATPase_end; Q_ATPase(end)];
    DeltaG0_ATPase_end = [DeltaG0_ATPase_end; DeltaGo_ATPase(end) ];

end
fprintf('The range of the transformed gibbs free energy of ATP hydrolysis is: %.2f kJ/mol \n',range(DeltaG0_ATPase_end+Q_ATPase_end))

%% collect steady states
J_ATPase_steady = [];
J_OxP_steady = [];
J_LDH_steady = [];
J_PYK_steady = [];
J_Gly_steady = [];
J_CK_steady = [];

for k = 1: length(ATPase_act)
    J = X_all{2,3,k};

    J_ATPase_steady= [J_ATPase_steady;J(end,getfield(modelInfo.FluxID,'ATPASE_fiber'))*1e3];
    J_OxP_steady = [J_OxP_steady;J(end,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3];
    J_LDH_steady = [J_LDH_steady;J(end,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3];
    J_PYK_steady = [J_PYK_steady;J(end,getfield(modelInfo.FluxID,'PYK_fiber'))*1e3];
    J_CK_steady = [J_CK_steady;J(end,getfield(modelInfo.FluxID,'CK_fiber'))*1e3];
    J_Gly_steady = [J_Gly_steady;J(end,getfield(modelInfo.FluxID,'PYK_fiber'))*1e3+J(end,getfield(modelInfo.FluxID,'PGK_fiber'))*1e3-J(end,getfield(modelInfo.FluxID,'PFKa_fiber'))*1e3];
end

%% Fig. 6
fig = figure(600);
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;

% Fig. 6A (ATP synthesis)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1100 0.5838 0.3628 0.3412];
area(J_ATPase_steady',[16.*J_OxP_steady,J_Gly_steady, J_CK_steady])
yline(20,'k--', 'LineWidth', 0.75)
legend({'OxPhos' 'Glyc.' 'CK'},'Interpreter','latex')
ylabel('$J^k_{\mathrm{ATP}}$ (mM/min)','Interpreter','latex')
title('\bf{ATP synthesis}','Interpreter', 'latex')

for k = 1:length(loc)
    fprintf('The OxPhos contribution to ATP synthesis at ATPase = %.0f mM/min is: %.0f %% \n', ATPase_act(loc(k)), 16.*J_OxP_steady(loc(k))./J_ATPase_steady(loc(k)).*100)
end

% Fig. 6B (Transformed Gibbs energy of ATP hydrolysis)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.5703 0.5838 0.3628 0.3412];
hold on
plot(ATPase_act+0.48,DeltaG0_ATPase_end+Q_ATPase_end,'k','LineWidth',LineW )
for k = 1:length(loc)
    plot(ATPase_act(loc(k))+0.48,DeltaG0_ATPase_end(loc(k))+Q_ATPase_end(loc(k)),'o','MarkerFaceColor',color(k,:),'MarkerEdgeColor','none')
end
xlabel({'ATP demand (mM/min)'},'Interpreter', 'latex')
ylabel({'$\Delta _\mathrm{r} G^{\prime}$ (kJ/mol)'},'Interpreter', 'latex')
ylim([-70 -40])
title({'\bf{Transformed Gibbs energy of ATP hydrolysis}'; '(t = 10 min)'},'Interpreter', 'latex')
legend({'' 'low' 'medium' 'high'},'Interpreter','latex')


% Fig. 6C (NADH oxidation/PYR uptake)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.11 0.1100 0.3628 0.3412];
area(J_ATPase_steady',[J_OxP_steady,J_LDH_steady])
xlabel('ATP demand (mM/min)','Interpreter', 'latex')
ylabel('$J^{k}_\mathrm{NADH,PYR}$ (mM/min)','Interpreter','latex')
title('\bf{NADH oxidation/PYR uptake}','Interpreter', 'latex')
yline(1.25,'k--', 'LineWidth', 0.75)
legend({'OxPhos' 'LDH'},'Interpreter','latex')

% Fig. 6D (Transformed Gibbs energy of the redox half reaction)
ax = subplot(N1,N2,4); box on; hold on;
ax.Position = [0.5703 0.11 0.3628 0.3412];
plot(ATPase_act+0.48,DeltaG0_redox_end+Q_redox_end,'k','LineWidth',LineW )
for k = 1:length(loc)
    plot(ATPase_act(loc(k))+0.48,DeltaG0_redox_end(loc(k))+Q_redox_end(loc(k)),'o','MarkerFaceColor',color(k,:),'MarkerEdgeColor','none')
end
xlabel({'ATP demand (mM/min)'},'Interpreter', 'latex')
ylabel({'$\Delta _\mathrm{r} G^{\prime}$ (kJ/mol)'},'Interpreter', 'latex')
ylim([20 50])
title({'\bf{Transformed Gibbs energy of the redox half reaction}'; '(t = 10 min)'},'Interpreter', 'latex')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%% Appendix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Suppl. Fig. 3A-C (Reaction rates associated with NADH)
fig = figure();
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;

% Suppl. Fig. 3A (Reaction rates associated with NADH at low intensity)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1100 0.5838 0.3628 0.3412];
k =1; % low intensity
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW/3)

end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$J^k_{\mathrm{NADH}}$ (mM/min)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-5 5])
title('\bf{Low intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% rectangle
Pos = ax.Position;
Xlim=[t1 t2];
Ylim=[-3 3];
dim = [0 -0.2 1.0 0.4];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.11+0.245 0.5838+0.235 0.3628/4 0.36/4],'Box','on');
hold on
% low intensity
k =1;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW/3)

end
plot([0 0],[-30 30],'--k')
xlim([-0.1 1])
ylim([-0.2 0.2])


% Suppl. Fig. 3B (Reaction rates associated with NADH at medium intensity)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.5703 0.5838 0.3628 0.3412];
k=2;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW/3)
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$J^k_{\mathrm{NADH}}$ (mM/min)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-10 10])
title('\bf{Medium intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% rectangle
Pos = ax.Position;
Xlim=[t1 t2];
Ylim=[-3 3];
dim = [0 -0.2 1.0 0.4];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.5703+0.245 0.5838+0.235 0.3628/4 0.36/4],'Box','on');
hold on
% medium intensity
k =2;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW/3)


end
plot([0 0],[-30 30],'--k')
xlim([-0.1 1])
ylim([-0.2 0.2])



% Suppl. Fig. 3C (Reaction rates associated with NADH at high intensity)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.11 0.1100 0.3628 0.3412];
k=3;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW/3)

end
xlim([t1 t2]);
ylim([-15 15]);
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel({'$J^k_{\mathrm{NADH}}$ (mM/min)'},'Interpreter', 'latex')
title('\bf{High intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% rectangle
Pos = ax.Position;
Xlim=[t1 t2];
Ylim=[-3 3];
dim = [0 -0.2 1.0 0.4];
dim_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(dim(1)-Xlim(1));
dim_conv(3)=(Pos(3))/(Xlim(2)-Xlim(1))*(dim(3));
dim_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(dim(2)-Ylim(1));
dim_conv(4)=(Pos(4))/(Ylim(2)-Ylim(1))*(dim(4));
annotation('rectangle',dim_conv,'Color','k','LineWidth',0.6)

% Zoom
ax1 = axes('Position',[0.11+0.245 0.11+0.235 0.3628/4 0.36/4],'Box','on');
hold on
% high intensity
k =3;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW/3)
end
plot([0 0],[-30 30],'--k')
xlim([-0.1 1])
ylim([-0.2 0.2])


%% Suppl. Fig. 3 D-F (Zoom: Reaction rates associated with NADH)

fig = figure();
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;

% ZOOM: Suppl. Fig. 3D (Reaction rates associated with NADH at low intensity)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1100 0.5838 0.3628 0.3412];
% low intensity
k =1;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW)

end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$J^{\mathrm{k}}_{\mathrm{NADH}}$ (mM/min)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-10e-3 10e-3])
title('\bf{Zoom: Low intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)


% ZOOM: Suppl. Fig. 3E (Reaction rates associated with NADH at medium intensity)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.5703 0.5838 0.3628 0.3412];
k=2;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW)

end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('$J^{\mathrm{k}}_{\mathrm{NADH}}$ (mM/min)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-10e-3 10e-3])
title('\bf{Zoom: Medium intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)


% ZOOM: Suppl. Fig. 3F (Reaction rates associated with NADH at high intensity)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.11 0.1100 0.3628 0.3412];
k=3;
for j = 1:num_phase
    x = X_all{j,3,loc(k)};
    t = X_all{j,1,loc(k)}-tend(1);
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW*1.5,'LineStyle',':')
    plot(t,x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-.')
    plot(t,-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','--')
    plot(t,-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW,'LineStyle','-')
    plot(t,x(:,getfield(modelInfo.FluxID,'G3PDH_fiber'))*1e3+x(:,getfield(modelInfo.FluxID,'GAPDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'LDH_fiber'))*1e3-x(:,getfield(modelInfo.FluxID,'OxPhosO2_fiber'))*1e3,'color',color(4,:),'LineWidth',LineW)

end
xlim([t1 t2]);
ylim([-10e-3 10e-3]);
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel({'$J^{\mathrm{k}}_{\mathrm{NADH}}$ (mM/min)'},'Interpreter', 'latex')
title('\bf{Zoom: High intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)



%% Suppl. Fig. 4 (lactate and pyruvate dynamics)
fig = figure();
fig.Position = [287.6667 104.3333 720.6667 420.0000];
N1 = 2;
N2 = 2;


% Suppl. Fig. 4A (low intensity)
ax = subplot(N1,N2,1); box on; hold on;
ax.Position = [0.1100 0.5838 0.3628 0.3412];
for k = 1
    for j = 1:num_phase
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW)
        plot(t,x(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(k,:),'LineWidth',1.5*LineW,'LineStyle',':')
    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('LAC and PYR (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-0.3 2])
title('\bf{Low intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% Suppl. Fig. 4B (medium intensity)
ax = subplot(N1,N2,2); box on; hold on;
ax.Position = [0.5703 0.5838 0.3628 0.3412];
for k = 2
    for j = 1:num_phase
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW)
        plot(t,x(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(k,:),'LineWidth',1.5*LineW,'LineStyle',':')
    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('LAC and PYR (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-1.5 10])
title('\bf{Medium intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)

% Suppl. Fig. 4C (high intensity)
ax = subplot(N1,N2,3); box on; hold on;
ax.Position = [0.11 0.1100 0.3628 0.3412];
for k = 3
    for j = 1:num_phase
        x = X_all{j,2,loc(k)};
        t = X_all{j,1,loc(k)}-tend(1);
        plot(t,x(:,getfield(modelInfo.SVarID,'pyruvate_fiber'))*1e3,'color',color(k,:),'LineWidth',LineW)
        plot(t,x(:,getfield(modelInfo.SVarID,'lactate_fiber'))*1e3,'color',color(k,:),'LineWidth',1.5*LineW,'LineStyle',':')
    end
end
xlabel('$t$ (min)','Interpreter', 'latex')
ylabel('LAC and PYR (mM)','Interpreter', 'latex')
xlim([t1 t2])
ylim([-6 40])
title('\bf{High intensity}','Interpreter', 'latex')
arrows(t1, t2, tend, ax)