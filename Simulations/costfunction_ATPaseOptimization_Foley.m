function [err] = costfunction_optimization_ga_Foley(p,settings,Data,w1,w2)

ATPase_Hz = p; % in M/min

%% model info
modelInfo = settings.modelInfo;

%% define Phases
tend = settings.tend;
num_phase = length(tend);

%% mandatory input parameters
x0 = settings.x0;
params = settings.params;
BX = settings.BX;
K_BX = settings.K_BX;
clamp_idx = settings.clamp_idx;

%% initialize for ODE loop
x0_phase = x0;
ATPase_on = settings.ATPase_on;

varlist = [ 1 : length( modelInfo.SVarList ) ];

X_all = cell(num_phase,3);
value = [];
val = 0;
tendold = 0;


x=[];
t = [];
%% Solving ODE
    for i =1:num_phase
        t1 = [];
        y = [];
        stopTime = now + 1*60 / 86400;
        options = odeset('MaxStep',5e-2,'NonNegative', varlist,'RelTol',1e-9,'AbsTol',1e-9,'OutputFcn', @(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx,flag) myOutputFcn(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx,flag, stopTime));
        params(getfield(modelInfo.ParID, 'x_ATPASE_fiber')) = ATPase_on(i)*ATPase_Hz;

        [t1,y] = ode15s(@dXdTMuscleMetabolism_OxPhos_FT, [tendold tendold+tend(i)], x0_phase, options, BX, K_BX, [],[],params,clamp_idx);
 
        X_all(i,1)={t1};
        X_all(i,2)={y};
        x0_phase = y(length(t1),:);


        x = [x;X_all{i,2}];
        t = [t;X_all{i,1}];

        if t(end) < tend(i)+tendold
            value = 1000;
            val =1;
            break
        end

           
        tendold = tendold+tend(i);
        x(end,:) = [];
        t(end,:) = [];

    end
    %%

    if val == 0
        t = t-tend(1); 
        x(t<-0.1,:)=[];
        t(t<-0.1)=[];


        PCr = 0;
        err_PCr = [];

  
        t_PCr = Data{2,1}(:,1); % sampling points data

        % simulated PCr
        PCr = interp1(t,x(:,getfield(modelInfo.SVarID,'phosphocreatine_fiber')),t_PCr); 
        PCr_steady = X_all{1,2}(end,getfield(modelInfo.SVarID,'phosphocreatine_fiber'))*1e3;

        % assume same inital intracellular concentrations
        PCr_data = Data{2,1}(:,2)./100*PCr_steady(1);

        % Loss
        err_PCr = sum(abs([PCr.*1e3 - PCr_data]./(max(PCr_data)-min(PCr_data))))./length(PCr_data);

        value = [value;err_PCr];
    end

err = sum(abs(value))/length(value);
end

%% End of simulation

function status = myOutputFcn(t, x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx, flag, stopTime)
status = double(now > stopTime);
end