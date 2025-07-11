function [err] = costfunction_optimization_ga(p,settings,Data,w1,w2)
% loss function for parameter identification

nH_OP = p(1);
Voxphos = p(2);
Kop_ADP = p(3);
ATPase_Hz(1) = p(4); % 2 Hz
ATPase_Hz(2) = p(5); % 4 Hz
ATPase_Hz(3) = p(6); % 10 Hz
buffer = p(7);
Pi_init = p(8);

%% define Phases
tend = settings.tend;
num_phase = length(tend);

% ATPase on/off
ATPase_on = settings.ATPase_on;

%% model info
modelInfo = settings.modelInfo;

%% inital conditions and parameters
x0 = settings.x0;
params = settings.params;
clamp_idx = settings.clamp_idx;
BX = settings.BX;
K_BX = settings.K_BX;

%% fitted parameters
params(getfield(modelInfo.ParID, 'nH_op_fiber')) = nH_OP ;
params(getfield(modelInfo.ParID, 'Kop_ADP_fiber')) = Kop_ADP ;
x0(getfield(modelInfo.SVarID, 'Pi_fiber')) = Pi_init;
n = length(BX);
BX(1:n) = buffer;

%% Make sure all the variables except the membrane potentials cannot become negative
varlist = [ 1 : length( modelInfo.SVarList ) ];

%% initialize for ODE loop
% store simulation results in X_all
X_all = cell(num_phase,3,size(Data,2));

value = [];
val = 0;

for k = 1:size(Data,2)
    tendold = 0;
    x0_phase = x0;
    x=[];
    t = [];

    for i =1:num_phase
        t1 = [];
        y = [];
        stopTime = now + 1*60 / 86400;
        options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-9,'AbsTol',1e-9,'OutputFcn', @(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx,flag) myOutputFcn(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx,flag, stopTime));

        % fitted parameters
        params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) = ATPase_on(i)*ATPase_Hz(k)*1e-3;
        params(getfield(modelInfo.ParID, 'x_OxPhosO2_fiber')) = Voxphos*1e-3/16;

        % solving ODEs
        [t1,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params,clamp_idx);

        X_all(i,1,k)={t1};
        X_all(i,2,k)={y};
        x0_phase = y(length(t1),:);

        x = [x;X_all{i,2,k}];
        t = [t;X_all{i,1,k}];

        if t(end) < tend(i)+tendold
            value = 1000; % model crashed
            val =1;
            break
        end

        tendold = tendold+tend(i);
        x(end,:) = [];
        t(end,:) = [];

    end
    %%

    if val == 0
        t = t-300;
        x(t<-0.1,:)=[];
        t(t<-0.1)=[];



        Pi = 0;
        PCr = 0;
        pH = 0;
        err_Pi = [];
        err_PCr = [];
        err_pH = [];

        t_pH  = Data{1,k}(:,1);
        t_PCr = Data{2,k}(:,1);
        t_Pi = Data{3,k}(:,1);

        Pi = interp1(t,x(:,getfield(modelInfo.SVarID,'Pi_fiber')),t_Pi);
        PCr = interp1(t,x(:,getfield(modelInfo.SVarID,'phosphocreatine_fiber')),t_PCr);
        pH = interp1(t,-log10(x(:,getfield(modelInfo.SVarID,'H_fiber'))),t_pH);


        % Cost function

        Pi_data = Data{3,k}(:,2);
        PCr_data = Data{2,k}(:,2);
        pH_data = Data{1,k}(:,2);

        pH_steady = -log10(X_all{1,2,k}(end,getfield(modelInfo.SVarID,'H_fiber')));
        Delta_pH = pH_steady-pH_data(1);


        % assume same initial intracellular concentrations
        Pi_data_adj = Pi_data./Pi_data(1).*(Pi(1).*1e3+1);
        PCr_data_adj= PCr_data./PCr_data(1).*(PCr(1).*1e3);

        err_Pi = sum(abs([(Pi.*1e3+1)- Pi_data_adj]./(max(Pi_data_adj)-min(Pi_data_adj))))./length(Pi_data);
        err_PCr = sum(abs([PCr.*1e3 - PCr_data_adj]./(max(PCr_data_adj)-min(PCr_data_adj))))./length(PCr_data);
        err_pH = sum(abs([pH- pH_data-Delta_pH]./(max(pH_data)-min(pH_data))))./length(pH_data);

        value = [value;err_Pi;err_pH;err_PCr];
    else
        break
    end


end
err = sum(abs(value))/length(value);
end



%% End of simulation

function status = myOutputFcn(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx, flag, stopTime)
status = double(now > stopTime);
end