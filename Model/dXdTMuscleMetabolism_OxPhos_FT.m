% function [f,J,Kapp,] = dXdTMuscleMetabolism_OxPhos_FT(t,x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx)
%
%
%
% Output parameters:
%   f     time derivatives of the model
%   J     flux
%   Kapp    apparent equilibrium constants
%
%
% Mandatory input parameters:
%   t     time
%   x     state variables at t=0
%   BX    buffer sizes
%   K_BX  proton buffer dissociation constants ( fiber extracellular capillary )
%   Xcp0_NADH    protein-NADH binding capacity
%   Kn_NADH    protine-NADH binding dissociation constant
%   par   parameter vector for the free parameters
%   clamp_idx   index vector to clamp selected state variabes
%   
%
% State Variables:
% [glucose1phos_fiber , glucose6phos_fiber , fructose6phos_fiber , ATP_fiber , fructose16phos_fiber , ADP_fiber , glyceraldehydephos_fiber , dihydroxyacetonephos_fiber , Pi_fiber , NAD_fiber , bpg_fiber , NADH_fiber , glycerol3phos_fiber , pg3_fiber , pg2_fiber , pep_fiber , pyruvate_fiber , phosphocreatine_fiber , creatine_fiber , AMP_fiber , lactate_fiber , O2aq_fiber , co2g_fiber , coaq_fiber , co2g_extracellular , coaq_extracellular , lactate_extracellular , co2g_capillary , coaq_capillary , lactate_capillary ]
%
%
% Free Parameters:
% [x_PGLM_fiber , x_PGI_fiber , x_PFKa_fiber , x_FBA_fiber , x_TPI_fiber , x_GAPDH_fiber , x_G3PDH_fiber , x_PGK_fiber , x_PGYM_fiber , x_ENO_fiber , x_PYK_fiber , x_ATPASE_fiber , x_CK_fiber , x_AK_fiber , x_LDH_fiber , Kop_ADP_fiber , Kop_Pi_fiber , x_OxPhosO2_fiber , nH_op_fiber , Kop_PYR_fiber , x_MCT_fiber_to_extracellular , Kmct_lac_fiber_to_extracellular , x_CO2Diff_fiber_to_extracellular , x_CO2weg_extracellular_to_capillary , x_HCO3weg_extracellular_to_capillary , x_LacWeg_extracellular_to_capillary ]

function [f,J,Kapp] = dXdTMuscleMetabolism_OxPhos_FT(t,x,BX,K_BX,Xcp0_NADH,Kn_NADH,par,clamp_idx)
%% GLOBAL VARIABLES
% temperature 37 
% ionic_strength 0.1 
T =310.15; 
I =0.1; 
%% LIST OF STATE VARIABLES
% 1 glucose1phos_fiber
% 2 glucose6phos_fiber
% 3 fructose6phos_fiber
% 4 ATP_fiber
% 5 fructose16phos_fiber
% 6 ADP_fiber
% 7 glyceraldehydephos_fiber
% 8 dihydroxyacetonephos_fiber
% 9 Pi_fiber
% 10 NAD_fiber
% 11 bpg_fiber
% 12 NADH_fiber
% 13 glycerol3phos_fiber
% 14 pg3_fiber
% 15 pg2_fiber
% 16 pep_fiber
% 17 pyruvate_fiber
% 18 phosphocreatine_fiber
% 19 creatine_fiber
% 20 AMP_fiber
% 21 lactate_fiber
% 22 O2aq_fiber
% 23 co2g_fiber
% 24 coaq_fiber
% 25 co2g_extracellular
% 26 coaq_extracellular
% 27 lactate_extracellular
% 28 co2g_capillary
% 29 coaq_capillary
% 30 lactate_capillary
% 31 h_fiber
% 32 m_fiber
% 33 k_fiber
% 34 h_extracellular
% 35 m_extracellular
% 36 k_extracellular
% 37 h_capillary
% 38 m_capillary
% 39 k_capillary
% 40 DPsi_fiber_to_extracellular
% 41 DPsi_extracellular_to_capillary

% PARTIAL VOLUME FRACTIONS
VWater_fiber = 0.754; % [=] l water (l region)^{-1}, Vinnakota and Bassingthwaighte, AJP, 2004 
VRegion_fiber = 0.875; % [=] l region (l tissue)^{-1}, Vinnakota and Bassingthwaighte, AJP, 2004 
VWater_extracellular = 1; % [=] l water (l region)^{-1}, Vinnakota and Bassingthwaighte, AJP, 2004 
VRegion_extracellular = 0.085; % [=] l region (l tissue)^{-1}, Vinnakota and Bassingthwaighte, AJP, 2004 
VWater_capillary = 1; % [=] l water (l region)^{-1}, Vinnakota and Bassingthwaighte, AJP, 2004 
VRegion_capillary = 0.04; % [=] l region (l tissue)^{-1}, Vinnakota and Bassingthwaighte, AJP, 2004 

%% THERMODYNAMIC DATA
RT = 8.314*T/1e3; % kJ  mol^{-1}
F = 0.096484; % kJ mol^{-1} mV^{-1}

%% STATE VARIABLES

% Concentrations of Reference Species
glucose1phos_fiber = x(1);
glucose6phos_fiber = x(2);
fructose6phos_fiber = x(3);
ATP_fiber = x(4);
fructose16phos_fiber = x(5);
ADP_fiber = x(6);
glyceraldehydephos_fiber = x(7);
dihydroxyacetonephos_fiber = x(8);
Pi_fiber = x(9);
NAD_fiber = x(10);
bpg_fiber = x(11);
NADH_fiber = x(12);
glycerol3phos_fiber = x(13);
pg3_fiber = x(14);
pg2_fiber = x(15);
pep_fiber = x(16);
pyruvate_fiber = x(17);
phosphocreatine_fiber = x(18);
creatine_fiber = x(19);
AMP_fiber = x(20);
lactate_fiber = x(21);
O2aq_fiber = x(22);
co2g_fiber = x(23);
coaq_fiber = x(24);
co2g_extracellular = x(25);
coaq_extracellular = x(26);
lactate_extracellular = x(27);
co2g_capillary = x(28);
coaq_capillary = x(29);
lactate_capillary = x(30);
% Concentrations of H, Mg, and K
h_fiber = x(31);
m_fiber = x(32);
k_fiber = x(33);
h_extracellular = x(34);
m_extracellular = x(35);
k_extracellular = x(36);
h_capillary = x(37);
m_capillary = x(38);
k_capillary = x(39);

% Membrane potentials
DPsi_fiber_to_extracellular = x(40);
DPsi_extracellular_to_capillary = x(41);

%% DISSOCIATION CONSTANTS
% glucose1phos_fiber
Kh(1) = 8.0708547586062666766027742445976933538531739031896e-07;
Km(1) = 0.0028545804190859376189837171011731697944924235343933;
Kk(1) = Inf;
% glucose6phos_fiber
Kh(2) = 7.7624711662869114663060960898621765124971716431901e-07;
Km(2) = Inf;
Kk(2) = Inf;
% fructose6phos_fiber
Kh(3) = 1.2882495516931349724339432236130953413066890789196e-06;
Km(3) = Inf;
Kk(3) = Inf;
% ATP_fiber
Kh(4) = 3.1841473060411534226195827933236781603909548721276e-07;
Km(4) = 5.2693462599814525461686121055038256599800661206245e-05;
Kk(4) = 0.069200292958055950598428296416386729106307029724121;
% fructose16phos_fiber
Kh(5) = 3.9810717055349692247584740718846507689931968343444e-07;
Km(5) = 0.0019952623149688789201683380980512083624489605426788;
Kk(5) = Inf;
% ADP_fiber
Kh(6) = 4.0957070424497166278660293184388230258718976983801e-07;
Km(6) = 0.00047168222467162355124661865524160475615644827485085;
Kk(6) = 0.1000000000000000055511151231257827021181583404541;
% glyceraldehydephos_fiber
Kh(7) = 3.5481338923357531451701621025285326993525814032182e-07;
Km(7) = Inf;
Kk(7) = Inf;
% dihydroxyacetonephos_fiber
Kh(8) = 1.2589254117941661142835753361968187391539686359465e-06;
Km(8) = 0.026915348039269152563557341295563674066215753555298;
Kk(8) = Inf;
% Pi_fiber
Kh(9) = 1.9001128754107379615305569928174200811099581187591e-07;
Km(9) = 0.022244795193035578340090552273977664299309253692627;
Kk(9) = 0.31622776601683794117647607890830840915441513061523;
% NAD_fiber
Kh(10) = Inf;
Km(10) = Inf;
Kk(10) = Inf;
% bpg_fiber
Kh(11) = 3.1622776601683791911296186882829317710275063291192e-08;
Km(11) = Inf;
Kk(11) = Inf;
% NADH_fiber
Kh(12) = Inf;
Km(12) = Inf;
Kk(12) = Inf;
% glycerol3phos_fiber
Kh(13) = 5.8536851325533857077701308316286521460369840497151e-07;
Km(13) = 0.023442288153199226236056418315456539858132600784302;
Kk(13) = Inf;
% pg3_fiber
Kh(14) = 6.1659500186148219147094264036557120789439068175852e-07;
Km(14) = Inf;
Kk(14) = Inf;
% pg2_fiber
Kh(15) = 9.9999999999999995474811182588625868561393872369081e-08;
Km(15) = 0.003548133892335753221403127355415563215501606464386;
Kk(15) = 0.066069344800759613467455722002341644838452339172363;
% pep_fiber
Kh(16) = 4.4668359215096348637208086354566383135988871799782e-07;
Km(16) = 0.0054954087385762481407502910712992161279544234275818;
Kk(16) = 0.083176377110267082914951686234417138621211051940918;
% pyruvate_fiber
Kh(17) = 0.0032359365692962811442145998341857193736359477043152;
Km(17) = Inf;
Kk(17) = Inf;
% phosphocreatine_fiber
Kh(18) = 3.3610474547182032044583682717231454262218903750181e-05;
Km(18) = 0.029675916796848916812123775343934539705514907836914;
Kk(18) = 0.48977881936844619437110281978675629943609237670898;
% creatine_fiber
Kh(19) = 0.0050118723362727246248282675367136107524856925010681;
Km(19) = Inf;
Kk(19) = Inf;
% AMP_fiber
Kh(20) = 4.9900767775049894178177479348024192518096242565662e-07;
Km(20) = 0.01111850542291117731330540863154965336434543132782;
Kk(20) = Inf;
% lactate_fiber
Kh(21) = 0.00021477547471276953053656577630192714423174038529396;
Km(21) = 0.104712854805089961018893518485128879547119140625;
Kk(21) = Inf;
% O2aq_fiber
Kh(22) = Inf;
Km(22) = Inf;
Kk(22) = Inf;
% co2g_fiber
Kh(23) = Inf;
Km(23) = Inf;
Kk(23) = Inf;
% coaq_fiber
Kh(24) = Inf;
Km(24) = Inf;
Kk(24) = Inf;
% co2g_extracellular
Kh(25) = Inf;
Km(25) = Inf;
Kk(25) = Inf;
% coaq_extracellular
Kh(26) = Inf;
Km(26) = Inf;
Kk(26) = Inf;
% lactate_extracellular
Kh(27) = 0.00021477547471276953053656577630192714423174038529396;
Km(27) = 0.104712854805089961018893518485128879547119140625;
Kk(27) = Inf;
% co2g_capillary
Kh(28) = Inf;
Km(28) = Inf;
Kk(28) = Inf;
% coaq_capillary
Kh(29) = Inf;
Km(29) = Inf;
Kk(29) = Inf;
% lactate_capillary
Kh(30) = 0.00021477547471276953053656577630192714423174038529396;
Km(30) = 0.104712854805089961018893518485128879547119140625;
Kk(30) = Inf;

%% BINDING POLYNOMIALS
P( 1 ) = 1  + h_fiber/Kh(1) + m_fiber/Km(1) + k_fiber/Kk(1);
P( 2 ) = 1  + h_fiber/Kh(2) + m_fiber/Km(2) + k_fiber/Kk(2);
P( 3 ) = 1  + h_fiber/Kh(3) + m_fiber/Km(3) + k_fiber/Kk(3);
P( 4 ) = 1  + h_fiber/Kh(4) + m_fiber/Km(4) + k_fiber/Kk(4);
P( 5 ) = 1  + h_fiber/Kh(5) + m_fiber/Km(5) + k_fiber/Kk(5);
P( 6 ) = 1  + h_fiber/Kh(6) + m_fiber/Km(6) + k_fiber/Kk(6);
P( 7 ) = 1  + h_fiber/Kh(7) + m_fiber/Km(7) + k_fiber/Kk(7);
P( 8 ) = 1  + h_fiber/Kh(8) + m_fiber/Km(8) + k_fiber/Kk(8);
P( 9 ) = 1  + h_fiber/Kh(9) + m_fiber/Km(9) + k_fiber/Kk(9);
P( 10 ) = 1  + h_fiber/Kh(10) + m_fiber/Km(10) + k_fiber/Kk(10);
P( 11 ) = 1  + h_fiber/Kh(11) + m_fiber/Km(11) + k_fiber/Kk(11);
P( 12 ) = 1  + h_fiber/Kh(12) + m_fiber/Km(12) + k_fiber/Kk(12);
P( 13 ) = 1  + h_fiber/Kh(13) + m_fiber/Km(13) + k_fiber/Kk(13);
P( 14 ) = 1  + h_fiber/Kh(14) + m_fiber/Km(14) + k_fiber/Kk(14);
P( 15 ) = 1  + h_fiber/Kh(15) + m_fiber/Km(15) + k_fiber/Kk(15);
P( 16 ) = 1  + h_fiber/Kh(16) + m_fiber/Km(16) + k_fiber/Kk(16);
P( 17 ) = 1  + h_fiber/Kh(17) + m_fiber/Km(17) + k_fiber/Kk(17);
P( 18 ) = 1  + h_fiber/Kh(18) + m_fiber/Km(18) + k_fiber/Kk(18);
P( 19 ) = 1  + h_fiber/Kh(19) + m_fiber/Km(19) + k_fiber/Kk(19);
P( 20 ) = 1  + h_fiber/Kh(20) + m_fiber/Km(20) + k_fiber/Kk(20);
P( 21 ) = 1  + h_fiber/Kh(21) + m_fiber/Km(21) + k_fiber/Kk(21);
P( 22 ) = 1  + h_fiber/Kh(22) + m_fiber/Km(22) + k_fiber/Kk(22);
P( 23 ) = 1  + h_fiber/Kh(23) + m_fiber/Km(23) + k_fiber/Kk(23);
P( 24 ) = 1  + h_fiber/Kh(24) + m_fiber/Km(24) + k_fiber/Kk(24);
P( 25 ) = 1  + h_extracellular/Kh(25) + m_extracellular/Km(25) + k_extracellular/Kk(25);
P( 26 ) = 1  + h_extracellular/Kh(26) + m_extracellular/Km(26) + k_extracellular/Kk(26);
P( 27 ) = 1  + h_extracellular/Kh(27) + m_extracellular/Km(27) + k_extracellular/Kk(27);
P( 28 ) = 1  + h_capillary/Kh(28) + m_capillary/Km(28) + k_capillary/Kk(28);
P( 29 ) = 1  + h_capillary/Kh(29) + m_capillary/Km(29) + k_capillary/Kk(29);
P( 30 ) = 1  + h_capillary/Kh(30) + m_capillary/Km(30) + k_capillary/Kk(30);

%% THERMODYNAMIC EQUATIONS
DGro_PGLM =-7.07;
DGro_PGI =3.14;
DGro_PFKa =17.6975;
DGro_FBA =21.4366;
DGro_TPI =-7.66;
DGro_GAPDH =43.8651;
DGro_G3PDH =63.1917;
DGro_PGK =-8.37;
DGro_PGYM =6.16;
DGro_ENO =-4.46;
DGro_PYK =-69.3658;
DGro_ATPASE =4.2842;
DGro_CK =-568.1942;
DGro_AK =-2.4858;
DGro_LDH =-64.6517;
DGro_OxPhosO2 =-1559.104;
DGro_CarbAnhyd =511.0379;
DGro_CO2Hyd =511.0379;
DGro_MCT =0;
DGro_CO2Diff =0;
DGro_CO2weg =0;
DGro_HCO3weg =0;
DGro_LacWeg =0;

Keq_PGLM_fiber = exp(-DGro_PGLM/RT)/P(1)*P(2);
Keq_PGI_fiber = exp(-DGro_PGI/RT)/P(2)*P(3);
Keq_PFKa_fiber = exp(-DGro_PFKa/RT)/P(3)/P(4)*P(5)*P(6)/h_fiber;
Keq_FBA_fiber = exp(-DGro_FBA/RT)/P(5)*P(7)*P(8);
Keq_TPI_fiber = exp(-DGro_TPI/RT)/P(7)*P(8);
Keq_GAPDH_fiber = exp(-DGro_GAPDH/RT)/P(7)/P(9)/P(10)*P(11)*P(12)/h_fiber;
Keq_G3PDH_fiber = exp(-DGro_G3PDH/RT)*P(8)/P(10)*P(12)/P(13)/h_fiber;
Keq_PGK_fiber = exp(-DGro_PGK/RT)*P(4)/P(6)/P(11)*P(14);
Keq_PGYM_fiber = exp(-DGro_PGYM/RT)/P(14)*P(15);
Keq_ENO_fiber = exp(-DGro_ENO/RT)/P(15)*P(16);
Keq_PYK_fiber = exp(-DGro_PYK/RT)*P(4)/P(6)/P(16)*P(17)*h_fiber;
Keq_ATPASE_fiber = exp(-DGro_ATPASE/RT)/P(4)*P(6)*P(9)/h_fiber;
Keq_CK_fiber = exp(-DGro_CK/RT)*P(4)/P(6)/P(18)*P(19)*h_fiber;
Keq_AK_fiber = exp(-DGro_AK/RT)/P(4)*P(6)^2/P(20);
Keq_LDH_fiber = exp(-DGro_LDH/RT)*P(10)/P(12)/P(17)*P(21)*h_fiber;
Keq_OxPhosO2_fiber = exp(-DGro_OxPhosO2/RT)*P(4)^16/P(6)^16/P(9)^16*P(10)/P(12)/P(17)/P(22)^3*P(23)^3*h_fiber^-18;
Keq_CarbAnhyd_fiber = exp(-DGro_CarbAnhyd/RT)/P(23)*P(24)/h_fiber;
Keq_CO2Hyd_extracellular = exp(-DGro_CO2Hyd/RT)/P(25)*P(26)/h_extracellular;
Keq_MCT = exp(-DGro_MCT/RT)/P(21)*P(27)*h_fiber^1/h_extracellular^1;
Keq_CO2Diff = exp(-DGro_CO2Diff/RT)/P(23)*P(25);
Keq_CO2weg = exp(-DGro_CO2weg/RT)/P(25)*P(28);
Keq_HCO3weg = exp(-DGro_HCO3weg/RT)/P(26)*P(29);
Keq_LacWeg = exp(-DGro_LacWeg/RT)/P(27)*P(30);

Kapp = [Keq_PGLM_fiber Keq_PGI_fiber Keq_PFKa_fiber Keq_FBA_fiber Keq_TPI_fiber Keq_GAPDH_fiber Keq_G3PDH_fiber Keq_PGK_fiber Keq_PGYM_fiber Keq_ENO_fiber Keq_PYK_fiber Keq_ATPASE_fiber Keq_CK_fiber Keq_AK_fiber Keq_LDH_fiber Keq_OxPhosO2_fiber Keq_CarbAnhyd_fiber Keq_CO2Hyd_extracellular ];
%% FLUX EQUATIONS
 
% PGLM_fiber: E.PGLM.0
 
%glucose1phos_fiber=glucose6phos_fiber
x_PGLM = par(1);
a=glucose1phos_fiber;
p=glucose6phos_fiber;
Kpglm_g1p=6.300000e-005;
Kpglm_g6p=3.000000e-005;
Vffpglm=x_PGLM;
pH_cy=-log10(h_fiber);
Vfpglm=Vffpglm.*(1.329)./(1+10.^(-pH_cy+6.64)+10.^(pH_cy-8.36));
Vbpglm=Vfpglm.*Kpglm_g6p./(Kpglm_g1p.*Keq_PGLM_fiber);
J_PGLM_fiber=(Vfpglm.*a./Kpglm_g1p-Vbpglm.*p./Kpglm_g6p)./...
(1+a./Kpglm_g1p+p./Kpglm_g6p);
 
% PGI_fiber: E.PGI.2
 
%glucose6phos_fiber=fructose6phos_fiber
a=glucose6phos_fiber;
p=fructose6phos_fiber;
x_PGI = par(2);
Vbbpgi=x_PGI;
pH_cy=-log10(h_fiber);
Kpgi_g6p=4.8e-4;
Kpgi_f6p=1.19e-4;
Vbpgi=Vbbpgi*(0.4203*pH_cy-2.3731);
Vfpgi=Vbpgi.*Kpgi_g6p.*Keq_PGI_fiber./Kpgi_f6p;
J_PGI_fiber=(Vfpgi.*a./Kpgi_g6p-Vbpgi.*p./Kpgi_f6p)./(1+p./Kpgi_f6p+a./Kpgi_g6p);
 
% PFKa_fiber: E.PFKa.1
 
Km1_f6p=68.49e-6;
%68.49muM
Km2_f6p=58.7e-6;
%M
Km1_MgATP=26.54e-6;
%M
Km2_MgATP=37.81e-6;
%M
Ki_f6p=4.33e-6;
%M
Ki_MgATP=124.6e-6;
%M
Ki_ATPH=0.649e-6;
%M
Ka=0.0812e-6;
%M
k_PFK=0.990;
c1=19.09;
c2=2.63;
Kadp_act=0.7e-6;
%adjusted
Kamp_act=1.23e-8;
%adjusted
x_PFKa = par(3);
%Vmax
%computebindingconstantsforATP
alpha_K=1.10708-1.54508e-3*298.15+5.95584e-6*298.15^2;
alpha_H=-1.28466e-5*298.15^2+9.90399e-8*298.15^3;
beta_K=(alpha_K/2.303)*((0.1^0.5/(1+1.6*0.1^0.5)-I^0.5/(1+1.6*I^0.5)));
beta_H=alpha_H*I^0.5/(1+1.6*I^0.5);
zm=2^2+(-4)^2-(-4+2)^2;
pKm=4.19+beta_K*zm;
drH_m=-18+beta_H*zm;
pKm=pKm+(1/T-1/298.15)*drH_m/(2.3026*8.314e-3);
KmATP=10^(-pKm);
zh=1^2+(-4)^2-(-4+1)^2;
pKa=6.48+beta_K*zh;
drH_h=-5+beta_H*zh;
pKa=pKa+(1/T-1/298.15)*drH_h/(2.3026*8.314e-3);
KhATP=10^(-pKa);
MgATP=m_fiber.*ATP_fiber./KmATP;
ATPH=ATP_fiber.*h_fiber./KhATP;
deinhibition=1+ADP_fiber./Kadp_act+AMP_fiber./Kamp_act;
Qi_MgATP=Ki_MgATP.*deinhibition;
Qi_ATPH=Ki_ATPH.*deinhibition;
QA=fructose6phos_fiber./Km1_f6p;
QB=MgATP./Km1_MgATP;
Q7=h_fiber./Ka.*(1+(ATPH./Qi_ATPH).^4);
Q1=QA.*(1+Q7./c1);
Q2=QB.*(1+Q7./c2);
Q=Q7.*(1+QA./c1).^3.*(1+QB./c2).^3./((1+QA).^3.*(1+QB).^3);
Q3=(1+Q./c1).*(Km2_MgATP./MgATP);
Q4=(1+Q/c2).*(Km2_f6p./fructose6phos_fiber);
Q5=(Ki_f6p./fructose6phos_fiber).*(1+Q).*Q3.*QA;
Q6=(Qi_MgATP./MgATP).*(1+Q).*Q4.*QB;
E=(1+Q./(c1.*c2)+Q5).*(1+Q./(c1.*c2)+Q6).*(1+Q7);
EA=Q1.*(1+Q./(c1.*c2)+Q6).*Q3.*(1+Q./c1);
EB=Q2.*(1+Q./(c1.*c2)+Q5).*Q4.*(1+Q./c2);
NUM=Q1.*(1+Q./(c1.*c2)+Q6)+Q2.*(1+Q./(c1.*c2)+Q5).*(1+Q./(c1.*c2));
EABBA=NUM.*(1+Q./(c1.*c2));
J_PFKa_fiber=x_PFKa.*k_PFK.*NUM./(E+EA+EABBA);
 
% FBA_fiber: E.FBA.2
 
%fructose16phos_fiber=glyceraldehydephos_fiber+dihydroxyacetonephos_fiber
a=fructose16phos_fiber;
p=glyceraldehydephos_fiber;
q=dihydroxyacetonephos_fiber;
x_FBA = par(4);
Vffald=x_FBA;
Kald_fbp=5.000000e-005;
Kald_dhap=2.000000e-003;
Kald_gap=1.000000e-003;
pH_cy=-log10(h_fiber);
Vfald=Vffald.*(0.5+0.5.*sin(pi./3.2.*(pH_cy+0.25)));
Vbald=Vfald.*Kald_gap.*Kald_dhap./(Kald_fbp.*Keq_FBA_fiber);
J_FBA_fiber=(Vfald.*a./Kald_fbp-Vbald.*p.*q./(Kald_gap.*Kald_dhap))./...
(1+a./Kald_fbp+p./Kald_gap+q./Kald_dhap);
 
% TPI_fiber: E.TPI.1
 
%dihydroxyacetonephos_fiber=glyceraldehydephos_fiber
a=glyceraldehydephos_fiber;
p=dihydroxyacetonephos_fiber;
x_TPI = par(5);
Vfftpi=x_TPI;
Ktpi_gap=3.200000e-004;
Ktpi_dhap=6.100000e-004;
Vftpi=Vfftpi;
Vbtpi=Vftpi.*Ktpi_dhap./(Ktpi_gap.*Keq_TPI_fiber);
J_TPI_fiber=(Vftpi.*a./Ktpi_gap-Vbtpi.*p./Ktpi_dhap)./...
(1+a./Ktpi_gap+p./Ktpi_dhap);
 
% GAPDH_fiber: E.GAPDH.1
 
%glyceraldehydephosphate+Pi_fiber+NAD_fiber=1,3-bisphospholycerate+NADH_fiber+H_fiber
a=glyceraldehydephos_fiber;
b=Pi_fiber;
c=NAD_fiber;
p=bpg_fiber;
q=NADH_fiber;
x_GAPDH = par(6);
Kgapdh_gap=2.500000e-006;
Kgapdh_nad=9.000000e-005;
Kgapdh_pi=2.900000e-004;
Kgapdh_bpg=8.000000e-007;
Kgapdh_nadh=3.300000e-006;
Vffgad=x_GAPDH;
pH_cy=-log10(h_fiber);
Dgap=1+b./Kgapdh_pi+a./Kgapdh_gap+c./Kgapdh_nad+...
a.*c./(Kgapdh_gap.*Kgapdh_nad)+a.*c.*b./(Kgapdh_gap.*...
Kgapdh_nad.*Kgapdh_pi)+p./Kgapdh_bpg+q./Kgapdh_nadh...
+p.*q./(Kgapdh_nadh.*Kgapdh_bpg);
Vfgad=Vffgad.*(0.0007.*exp(pH_cy.*0.8979));
Vbgad=Vfgad.*Kgapdh_bpg.*Kgapdh_nadh./...
(Kgapdh_gap.*Kgapdh_pi.*Kgapdh_nad.*Keq_GAPDH_fiber);
J_GAPDH_fiber=(Vfgad.*a.*c.*b./(Kgapdh_nad.*Kgapdh_gap.*Kgapdh_pi)-...
Vbgad.*p.*q./(Kgapdh_bpg.*Kgapdh_nadh))./Dgap;
 
% G3PDH_fiber: E.G3PDH.0
 
%glycerol3phos_fiber+NAD_fiber=dihydroxyacetonephos_fiber+NADH_fiber+H_fiber
x_G3PDH = par(7);
a=glycerol3phos_fiber;
b=NAD_fiber;
p=dihydroxyacetonephos_fiber;
q=NADH_fiber;
Kg3pdh_g3p=1.800000e-004;
Kg3pdh_nad=1.200000e-005;
Kg3pdh_dhap=2.200000e-004;
Kg3pdh_nadh=8.000000e-006;
Vbbg3pdh=x_G3PDH;
Dg3pdh=(1+a./Kg3pdh_g3p+q./Kg3pdh_nadh).*...
(1+p./Kg3pdh_dhap+b./Kg3pdh_nad);
Vbg3pdh=Vbbg3pdh;
Vfg3pdh=Vbg3pdh.*Kg3pdh_g3p.*Kg3pdh_nad.*Keq_G3PDH_fiber./...
(Kg3pdh_dhap.*Kg3pdh_nadh);
J_G3PDH_fiber=(Vfg3pdh.*a.*b./(Kg3pdh_g3p.*Kg3pdh_nad)-...
Vbg3pdh.*p.*q./(Kg3pdh_dhap.*Kg3pdh_nadh))./Dg3pdh;
 
% PGK_fiber: E.PGK.1
 
%bpg_fiber+ADP_fiber=pg3_fiber+ATP_fiber
a=bpg_fiber;
b=ADP_fiber;
p=ATP_fiber;
q=pg3_fiber;
x_PGK = par(8);
Kpgk_bpg=2.000000e-003;
Kpgk_adp=8.000000e-006;
Kpgk_3pg=1.200000e-003;
Kpgk_atp=3.500000e-004;
Vbbpgk=x_PGK;
Vbpgk=Vbbpgk;
Vfpgk=Vbpgk.*Kpgk_bpg.*Kpgk_adp.*Keq_PGK_fiber./(Kpgk_3pg.*Kpgk_atp);
D_PGK=(1+b./Kpgk_adp+a./Kpgk_bpg+a.*b./(Kpgk_bpg.*Kpgk_adp)+...
q./Kpgk_3pg+p./Kpgk_atp+q.*p./(Kpgk_3pg.*Kpgk_atp));
J_PGK_fiber=(Vfpgk.*a.*b./(Kpgk_adp.*Kpgk_bpg)-...
Vbpgk.*p.*q./(Kpgk_atp.*Kpgk_3pg))./D_PGK;
 
% PGYM_fiber: E.PGYM.2
 
%pg3_fiber=pg2_fiber
a=pg3_fiber;
p=pg2_fiber;
x_PGYM = par(9);
Kpgm_3pg=2.000000e-004;
Kpgm_2pg=1.400000e-005;
pH_cy=-log10(h_fiber);
Vffpgm=x_PGYM;
Vfpgm=Vffpgm.*(0.4444.*pH_cy-2.2407);
Vbpgm=Vfpgm.*Kpgm_2pg./(Kpgm_3pg.*Keq_PGYM_fiber);
J_PGYM_fiber=(Vfpgm.*a./Kpgm_3pg-Vbpgm.*p./Kpgm_2pg)./...
(1+a./Kpgm_3pg+p./Kpgm_2pg);
 
% ENO_fiber: E.ENO.1
 
%pg2_fiber=pep_fiber+H2O_fiber
a=pg2_fiber;
p=pep_fiber;
x_ENO = par(10);
Ken_2pg=1.000000e-004;
Ken_pep=3.700000e-004;
Vffen=x_ENO;
Vfen=Vffen;
Vben=Vfen.*Ken_pep./(Ken_2pg.*Keq_ENO_fiber);
J_ENO_fiber=(Vfen.*a./Ken_2pg-Vben.*p./Ken_pep)./...
(1+p./Ken_pep+a./Ken_2pg);
 
% PYK_fiber: E.PYK.2
 
%pep_fiber+ADP_fiber+H_fiber=pyruvate_fiber+ATP_fiber
a=pep_fiber;
b=ADP_fiber;
p=pyruvate_fiber;
q=ATP_fiber;
x_PYK = par(11);
Kpk_pep=8.000000e-005;
Kpk_adp=3.000000e-004;
Kpk_pyr=7.050000e-003;
Kpk_atp=1.130000e-003;
Vffpk=x_PYK;
pH_cy=-log10(h_fiber);
Vfpk=Vffpk./(1+10.^(pH_cy-8.2));
Vbpk=Vfpk.*Kpk_pyr.*Kpk_atp./(Kpk_pep.*Kpk_adp.*Keq_PYK_fiber);
J_PYK_fiber=(Vfpk.*a.*b./(Kpk_pep.*Kpk_adp)-...
Vbpk.*p.*q./(Kpk_pyr.*Kpk_atp))./...
(1+a./Kpk_pep+b./Kpk_adp+a.*b./(Kpk_pep.*Kpk_adp)...
+q./Kpk_atp+p./Kpk_pyr+p.*q./(Kpk_pyr.*Kpk_atp));
 
% ATPASE_fiber: E.ATPASE.3
 
x_ATPASE = par(12);
J_ATPASE_fiber=0.480*1e-3+x_ATPASE;
 
% CK_fiber: E.CK.5
 
%phosphocreatine_fiber+ADP_fiber+H_fiber=creatine_fiber+ATP_fiber
a=phosphocreatine_fiber;
b=ADP_fiber;
p=ATP_fiber;
q=creatine_fiber;
x_CK = par(13);
Kck_pcr=1.110000e-003;
Kck_iatp=3.500000e-003;
Kck_iadp=1.350000e-004;
Kck_ipcr=3.900000e-003;
Kck_cr=3.800000e-003;
VrevCK=x_CK;
%EquilibriumsconstantforCK:K_CK=1.66e9;
K_CK=1.66e9;
Keq_CK_fiber=K_CK*h_fiber;
VforCK=(VrevCK.*Keq_CK_fiber)./(Kck_iatp.*Kck_cr./(Kck_iadp.*Kck_pcr));
J_CK_fiber=-(VrevCK.*p.*q./(Kck_iatp.*Kck_cr)-(VforCK.*b.*a./...
(Kck_iadp.*Kck_pcr)))./(1+b./Kck_iadp+a./Kck_ipcr+a.*b./...
(Kck_iadp.*Kck_pcr)+p./Kck_iatp+q.*p./(Kck_cr.*Kck_iatp));
 
% AK_fiber: E.AK.1
 
%ATP_fiber+AMP_fiber=2ADP
a=ADP_fiber;
p=ATP_fiber;
q=AMP_fiber;
Kadk_amp=3.200000e-004;
Kadk_atp=2.700000e-004;
Kadk_adp=3.500000e-004;
x_AK = par(14);
Vfadk=x_AK;
Vbadk=Vfadk.*Kadk_adp.^2./(Kadk_amp.*Kadk_atp.*Keq_AK_fiber);
J_AK_fiber=(Vfadk.*p.*q./(Kadk_atp.*Kadk_amp)-Vbadk.*(a./...
Kadk_adp).^2)./(1+p./Kadk_atp+q./Kadk_amp+p.*q./...
(Kadk_atp.*Kadk_amp)+2.*a./Kadk_adp+a.^2./Kadk_adp.^2);
 
% LDH_fiber: E.LDH.1
 
%pyruvate_fiber+NADH_fiber+H_fiber=lactate_fiber+NAD_fiber
x_LDH = par(15);
a=pyruvate_fiber;
b=NADH_fiber;
p=lactate_fiber;
q=NAD_fiber;
Kldh_pyr=3.350000e-004;
Kldh_nadh=2.000000e-006;
Kldh_lac=1.700000e-002;
Kldh_nad=8.490000e-004;
Vffldh=x_LDH;
pH_cy=-log10(h_fiber);
Vfldh=Vffldh.*(-0.1376.*pH_cy+1.745);
Vbldh=Vfldh.*Kldh_lac.*Kldh_nad./(Kldh_pyr.*Kldh_nadh.*Keq_LDH_fiber);
J_LDH_fiber=(Vfldh.*a.*b./(Kldh_pyr.*Kldh_nadh)-...
Vbldh.*p.*q./(Kldh_lac.*Kldh_nad))./...
(1+a./Kldh_pyr+b./Kldh_nadh+...
a.*b./(Kldh_pyr.*Kldh_nadh)+p./Kldh_lac+q./Kldh_nad+...
p.*q./(Kldh_lac.*Kldh_nad));
 
% OxPhosO2_fiber: E.OxPhosO2.2
 
Kop_ADP = par(16);
Kop_Pi = par(17);
x_OxPhosO2 = par(18);
nH_op = par(19);
Kop_PYR = par(20);
J_OxPhosO2_fiber=x_OxPhosO2.*(ADP_fiber./Kop_ADP)^nH_op./(1+(ADP_fiber./Kop_ADP)^nH_op).*1/(1+Kop_PYR/pyruvate_fiber).*1/(1+Kop_Pi/Pi_fiber);
 
% CarbAnhyd_fiber: E.CarbAnhyd.0
 
X_CA=950;
kf_CA=0.0469;
kb_CA=9.89e4;
J_CarbAnhyd_fiber=X_CA*(kf_CA*co2g_fiber-kb_CA*h_fiber*coaq_fiber);
 
% CO2Hyd_extracellular: E.CO2Hyd.0
 
kf=0.0469;
kb=9.89e4;
J_CO2Hyd_extracellular=(kf*co2g_extracellular-kb*h_extracellular*coaq_extracellular);
%HCO3
 
% MCT : fiber_to_extracellular
 
a = lactate_fiber;
b = h_fiber;
c = lactate_extracellular;
d = h_extracellular;
x_MCT  = par(21);
Kmct_h = 10^(-8.89);
Kmct_lac  = par(22);
VmaxMCT = x_MCT;
num = a.*b./(Kmct_lac.*Kmct_h) - c.*d./(Kmct_lac.*Kmct_h);
den = 2 + a.*b./(Kmct_lac.*Kmct_h) + a./Kmct_lac + b./Kmct_h ...
+ c.*d./(Kmct_lac.*Kmct_h) + c./Kmct_lac + d./Kmct_h;
J_MCT_fiber_to_extracellular = VmaxMCT.*(num./den);
 
% CO2Diff : fiber_to_extracellular
 
x_CO2Diff  = par(23);
J_CO2Diff_fiber_to_extracellular = x_CO2Diff*(co2g_fiber - co2g_extracellular);
% 2.3084
 
% CO2weg : extracellular_to_capillary
 
x_CO2weg  = par(24);
J_CO2weg_extracellular_to_capillary = - x_CO2weg*(1.5e-3 -co2g_extracellular);
 
% HCO3weg : extracellular_to_capillary
 
x_HCO3weg  = par(25);
J_HCO3weg_extracellular_to_capillary = - x_HCO3weg*(30e-3 -coaq_extracellular);
 
% LacWeg : extracellular_to_capillary
 
x_LacWeg  = par(26);
J_LacWeg_extracellular_to_capillary = x_LacWeg*lactate_extracellular*0.6517;

%% Variable stoichiometry for OxPhos
Ka_NADH = 0.1e-6;
s_NADH = (1/(1+(Ka_NADH/NADH_fiber)^4));
s_NAD = s_NADH;
s_O2 = (5+s_NADH)/2;
s_CO2 = 3;
s_ATP = (46+9*s_NADH)/(11/3)+1;
s_ADP = s_ATP;
s_Pi = s_ATP;
s_H = 1+s_ATP+s_NADH;
%% REACTANT TIME DERIVATIVES
f(1,:) = ( 0  ); % [clamped] % glucose1phos_fiber
f(2,:) = ( 0  + 1*J_PGLM_fiber - 1*J_PGI_fiber ) / VWater_fiber; % glucose6phos_fiber
f(3,:) = ( 0  + 1*J_PGI_fiber - 1*J_PFKa_fiber ) / VWater_fiber; % fructose6phos_fiber
f(4,:) = ( 0  - 1*J_PFKa_fiber + 1*J_PGK_fiber + 1*J_PYK_fiber - 1*J_ATPASE_fiber + 1*J_CK_fiber - 1*J_AK_fiber + s_ATP*J_OxPhosO2_fiber ) / VWater_fiber; % ATP_fiber
f(5,:) = ( 0  + 1*J_PFKa_fiber - 1*J_FBA_fiber ) / VWater_fiber; % fructose16phos_fiber
f(6,:) = ( 0  + 1*J_PFKa_fiber - 1*J_PGK_fiber - 1*J_PYK_fiber + 1*J_ATPASE_fiber - 1*J_CK_fiber + 2*J_AK_fiber - s_ADP*J_OxPhosO2_fiber ) / VWater_fiber; % ADP_fiber
f(7,:) = ( 0  + 1*J_FBA_fiber - 1*J_TPI_fiber - 1*J_GAPDH_fiber ) / VWater_fiber; % glyceraldehydephos_fiber
f(8,:) = ( 0  + 1*J_FBA_fiber + 1*J_TPI_fiber + 1*J_G3PDH_fiber ) / VWater_fiber; % dihydroxyacetonephos_fiber
f(9,:) = ( 0  - 1*J_PGLM_fiber - 1*J_GAPDH_fiber + 1*J_ATPASE_fiber - s_Pi*J_OxPhosO2_fiber ) / VWater_fiber; % Pi_fiber
f(10,:) = ( 0  - 1*J_GAPDH_fiber - 1*J_G3PDH_fiber + 1*J_LDH_fiber + s_NAD*J_OxPhosO2_fiber ) / VWater_fiber; % NAD_fiber
f(11,:) = ( 0  + 1*J_GAPDH_fiber - 1*J_PGK_fiber ) / VWater_fiber; % bpg_fiber
f(12,:) = ( 0  + 1*J_GAPDH_fiber + 1*J_G3PDH_fiber - 1*J_LDH_fiber - s_NADH*J_OxPhosO2_fiber ) / VWater_fiber; % NADH_fiber
f(13,:) = ( 0  - 1*J_G3PDH_fiber ) / VWater_fiber; % glycerol3phos_fiber
f(14,:) = ( 0  + 1*J_PGK_fiber - 1*J_PGYM_fiber ) / VWater_fiber; % pg3_fiber
f(15,:) = ( 0  + 1*J_PGYM_fiber - 1*J_ENO_fiber ) / VWater_fiber; % pg2_fiber
f(16,:) = ( 0  + 1*J_ENO_fiber - 1*J_PYK_fiber ) / VWater_fiber; % pep_fiber
f(17,:) = ( 0  + 1*J_PYK_fiber - 1*J_LDH_fiber - 1*J_OxPhosO2_fiber ) / VWater_fiber; % pyruvate_fiber
f(18,:) = ( 0  - 1*J_CK_fiber ) / VWater_fiber; % phosphocreatine_fiber
f(19,:) = ( 0  + 1*J_CK_fiber ) / VWater_fiber; % creatine_fiber
f(20,:) = ( 0  - 1*J_AK_fiber ) / VWater_fiber; % AMP_fiber
f(21,:) = ( 0  + 1*J_LDH_fiber - 1*J_MCT_fiber_to_extracellular/VRegion_fiber*VRegion_extracellular ) / VWater_fiber; % lactate_fiber
f(22,:) = ( 0  ); % [clamped] % O2aq_fiber
f(23,:) = ( 0  + s_CO2*J_OxPhosO2_fiber - 1*J_CarbAnhyd_fiber - 1*J_CO2Diff_fiber_to_extracellular/VRegion_fiber*VRegion_extracellular ) / VWater_fiber; % co2g_fiber
f(24,:) = ( 0  + 1*J_CarbAnhyd_fiber ) / VWater_fiber; % coaq_fiber
f(25,:) = ( 0  - 1*J_CO2Hyd_extracellular + 1*J_CO2Diff_fiber_to_extracellular - 1*J_CO2weg_extracellular_to_capillary/VRegion_extracellular*VRegion_capillary ) / VWater_extracellular; % co2g_extracellular
f(26,:) = ( 0  + 1*J_CO2Hyd_extracellular - 1*J_HCO3weg_extracellular_to_capillary/VRegion_extracellular*VRegion_capillary ) / VWater_extracellular; % coaq_extracellular
f(27,:) = ( 0  + 1*J_MCT_fiber_to_extracellular - 1*J_LacWeg_extracellular_to_capillary/VRegion_extracellular*VRegion_capillary ) / VWater_extracellular; % lactate_extracellular
f(28,:) = ( 0  + 1*J_CO2weg_extracellular_to_capillary ) / VWater_capillary; % co2g_capillary
f(29,:) = ( 0  + 1*J_HCO3weg_extracellular_to_capillary ) / VWater_capillary; % coaq_capillary
f(30,:) = ( 0  + 1*J_LacWeg_extracellular_to_capillary ) / VWater_capillary; % lactate_capillary

%% ION EQUATIONS
% COMPARTMENT fiber:
ii = [1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24]; % Indices of SVs in compartment fiber
% PARTIAL DERIVATIVES
pHBpM = -sum( (h_fiber*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpK = -sum( (h_fiber*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_fiber./Km(ii)+k_fiber./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_fiber*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_fiber*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_fiber./Kh(ii)+k_fiber./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_fiber*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_fiber*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_fiber./Kh(ii)+m_fiber./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0 + 0*J_PGLM_fiber + 0*J_PGI_fiber + 1*J_PFKa_fiber + 0*J_FBA_fiber + 0*J_TPI_fiber + 1*J_GAPDH_fiber + 1*J_G3PDH_fiber + 0*J_PGK_fiber + 0*J_PGYM_fiber + 0*J_ENO_fiber - 1*J_PYK_fiber + 1*J_ATPASE_fiber - 1*J_CK_fiber + 0*J_AK_fiber - 1*J_LDH_fiber - s_H*J_OxPhosO2_fiber + 1*J_CarbAnhyd_fiber - 1*J_MCT_fiber_to_extracellular/VRegion_fiber*VRegion_extracellular) / VWater_fiber;
J_M = (0) / VWater_fiber;
J_K = (0) / VWater_fiber;
Phi_H = J_H - sum( h_fiber*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = -sum( m_fiber*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K -sum( k_fiber*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(1)/K_BX(1)/(1+h_fiber/K_BX(1))^2 + 1e-14/h_fiber^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(31,:) =  ( (pKBpM*pMBpK - aM*aK)*Phi_H + ...
            (aK*pHBpM - pHBpK*pKBpM)*Phi_M + ...
            (aM*pHBpK - pHBpM*pMBpK)*Phi_K ) / D;
f(32,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(33,:) =  0;
% COMPARTMENT extracellular:
ii = [25  26  27]; % Indices of SVs in compartment extracellular
% PARTIAL DERIVATIVES
pHBpM = -sum( (h_extracellular*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpK = -sum( (h_extracellular*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_extracellular./Km(ii)+k_extracellular./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_extracellular*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_extracellular*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_extracellular./Kh(ii)+k_extracellular./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_extracellular*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_extracellular*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_extracellular./Kh(ii)+m_extracellular./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0 + 1*J_CO2Hyd_extracellular + 1*J_MCT_fiber_to_extracellular) / VWater_extracellular;
J_M = (0) / VWater_extracellular;
J_K = (0) / VWater_extracellular;
Phi_H = J_H - sum( h_extracellular*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = -sum( m_extracellular*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K -sum( k_extracellular*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(2)/K_BX(2)/(1+h_extracellular/K_BX(2))^2 + 1e-14/h_extracellular^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(34,:) =  0;
f(35,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(36,:) =  ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;
% COMPARTMENT capillary:
ii = [28  29  30]; % Indices of SVs in compartment capillary
% PARTIAL DERIVATIVES
pHBpM = -sum( (h_capillary*x(ii)'./Kh(ii))./(Km(ii).*P(ii).^2) );
pHBpK = -sum( (h_capillary*x(ii)'./Kh(ii))./(Kk(ii).*P(ii).^2) );
pHBpH = +sum( (1+m_capillary./Km(ii)+k_capillary./Kk(ii)).*x(ii)'./(Kh(ii).*P(ii).^2) );
pMBpH = -sum( (m_capillary*x(ii)'./Km(ii))./(Kh(ii).*P(ii).^2) );
pMBpK = -sum( (m_capillary*x(ii)'./Km(ii))./(Kk(ii).*P(ii).^2) );
pMBpM = +sum( (1+h_capillary./Kh(ii)+k_capillary./Kk(ii)).*x(ii)'./(Km(ii).*P(ii).^2) );
pKBpH = -sum( (k_capillary*x(ii)'./Kk(ii))./(Kh(ii).*P(ii).^2) );
pKBpM = -sum( (k_capillary*x(ii)'./Kk(ii))./(Km(ii).*P(ii).^2) );
pKBpK = +sum( (1+h_capillary./Kh(ii)+m_capillary./Km(ii)).*x(ii)'./(Kk(ii).*P(ii).^2) );
% PHIs
J_H = (0) / VWater_capillary;
J_M = (0) / VWater_capillary;
J_K = (0) / VWater_capillary;
Phi_H = J_H - sum( h_capillary*f(ii)'./(Kh(ii).*P(ii)) );
Phi_M = -sum( m_capillary*f(ii)'./(Km(ii).*P(ii)) );
Phi_K = J_K -sum( k_capillary*f(ii)'./(Kk(ii).*P(ii)) );
% ALPHAs
aH = 1 + pHBpH;
aM = 1 + pMBpM;
aK = 1 + pKBpK;
% ADDITIONAL BUFFER for [H+]
aH = 1 + pHBpH + BX(3)/K_BX(3)/(1+h_capillary/K_BX(3))^2 + 1e-14/h_capillary^2; % M
% DENOMINATOR
D = aH*pKBpM*pMBpK + aK*pHBpM*pMBpH + aM*pHBpK*pKBpH - ...
    aM*aK*aH - pHBpK*pKBpM*pMBpH - pHBpM*pMBpK*pKBpH;
% DERIVATIVES for H,Mg,K
f(37,:) =  0;
f(38,:) =  ( (aK*pMBpH - pKBpH*pMBpK)*Phi_H + ...
            (pKBpH*pHBpK - aH*aK)*Phi_M + ...
            (aH*pMBpK - pHBpK*pMBpH)*Phi_K ) / D;
f(39,:) =  ( (aM*pKBpH - pKBpM*pMBpH)*Phi_H + ...
            (aH*pKBpM - pKBpH*pHBpM)*Phi_M + ...
            (pMBpH*pHBpM - aH*aM)*Phi_K ) / D;

%% ELECTROPHYS EQUATIONS
C_fiber_to_extracellular = par(27);
f(40) = ( 0) / C_fiber_to_extracellular;
C_extracellular_to_capillary = par(28);
f(41) = ( 0) / C_extracellular_to_capillary;
%% FLUX VECTOR:
J = [ J_PGLM_fiber J_PGI_fiber J_PFKa_fiber J_FBA_fiber J_TPI_fiber J_GAPDH_fiber J_G3PDH_fiber J_PGK_fiber J_PGYM_fiber J_ENO_fiber J_PYK_fiber J_ATPASE_fiber J_CK_fiber J_AK_fiber J_LDH_fiber J_OxPhosO2_fiber J_CarbAnhyd_fiber J_CO2Hyd_extracellular J_MCT_fiber_to_extracellular J_CO2Diff_fiber_to_extracellular J_CO2weg_extracellular_to_capillary J_HCO3weg_extracellular_to_capillary J_LacWeg_extracellular_to_capillary];

f(clamp_idx,:) = 0;


