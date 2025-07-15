format long
clear
close all

%% 'Delete previously existing files'
delete brate.csv;
delete frate.csv;
delete coverage.csv;
delete rate.csv;
delete TOF;
delete combined_coverages.xls
delete combined_rate.xls
delete combined_DKRC_E.csv
delete combined_DTRC_E.csv
delete combined_DKRC_M.csv
delete combined_DTRC_M.csv
delete combined_sigma_E.csv
delete combined_sigma_M.csv
delete dlnTOF_E_dGts_thetaconst.csv
delete dlnTOF_E_dGads_thetaconst.csv
delete dlnTOF_E_dtheta_Gconst.csv
delete dlnTOF_M_dGts_thetaconst.csv
delete dlnTOF_M_dGads_thetaconst.csv
delete dlnTOF_M_dtheta_Gconst.csv
delete dthetadGts.csv
delete dthetadGads.csv
delete dthetadtheta0.csv
delete jacobian.csv
delete sites
delete time_evolution.csv
delete sum_DTRC_DKRC_M.csv
delete sum_DTRC_DKRC_E.csv
delete true_coverages.csv
delete sum_DTRC_M.csv
delete sum_DKRC_M.csv
delete sum_DTRC_E.csv
delete sum_DKRC_E.csv

%% 'Add CasADi to path'
addpath('/home/bamidele/per_bin/casadi-3.7.0-linux64-matlab2018b') ;

%% 'Set run time here in seconds' 
run_time = 2.40e-0 ;

%% 'Constants'
kb = 8.617333262e-5 ; % 'Boltzmanns constant'
h = 4.135667696e-15 ; % 'Plancks constant'

%% 'Operating conditions'
T = 573.15 ;

%% 'DFT Energy file names'
species_file = 'species_ads_111.csv' ;
ts_file = 'ts_energy_111_new_2.xlsx' ;

%% 'Read DFT energies and parameters'
[G_ads_111, G_ts_111, G_ads_act_111, G_ads_rxn_111, EZPEc_gas_111, lnQ_gas_111, LI_params_111] = read_Energies(T, species_file, ts_file, 111) ;

%% 'Add adsorption-desorption steps activation energy to G_ts_111'
G_ts_111(27:31) = G_ads_act_111 ;

%% 'Read Lateral Interactions Parameters'
LI_111_H = LI_params_111.('H_LI') ;

LI_111_H = ActivateLI(LI_111_H, "no") ; % 'Make all Lateral interactions zero'

H2_pressure = 0.1317225 ; % '0.13 atm to Bar'
Ethane_pressure = 0.03343725 ; % '0.033 atm to Bar'
coverages = zeros(length(G_ads_111), 1) ;

PCH3CH3 = Ethane_pressure ;
PCH2CH2 = 0;
PCHCH = 0;
PCH4 = 0;
PH2 = H2_pressure ;

Gases_0 = [PCH3CH3; PCH2CH2; PCHCH; PCH4; PH2] ;

Ptot = sum(Gases_0) ;

% 'Initialize surface coverages'
x0=zeros(length(G_ads_111),1);

% 'From clean surface: 100% free site'
x0(17) = 1;

% 'Gas mole fractions'
no_of_gases = length(EZPEc_gas_111) ;
x = zeros(no_of_gases, 1) ;
x(:) = Gases_0 / Ptot ;

%% 'Initialize relational variables and columnize'
Q_Gts0 = zeros(length(G_ts_111)*length(x0), 1) ;
Q_Gads0 = zeros(length(G_ads_111)*length(x0), 1) ;
R0 = reshape(eye(length(x0)), [length(x0)*length(x0), 1]) ;

%% 'Add Qs and R inital to y'
y0 = [x0 ; Q_Gts0; Q_Gads0; R0] ;

%% 'Mass matrix'
M=eye(length(y0)) ;
M(length(x0), length(x0)) = 0 ; % 'Compute free site balance as Algebraic eqn'

optode1 = odeset('NonNegative',1:17, 'Mass',M,'Abstol',1E-15,'RelTol',1E-15) ;
optode = odeset(optode1,'Abstol',1E-10,'RelTol',1E-10) ;
optlsq = optimset('TolFun',1E-13,'TolX',1E-13) ;      

%% 'Solve the composite ODE/DAE equation'
[t,y] = ode15s(@(t,y)C2_Pt_111_function_inline(t, y, x, Ptot, T, G_ts_111, G_ads_111, EZPEc_gas_111, lnQ_gas_111, LI_111_H),[0, run_time],y0,optode) ;

%% 'Write coverages and time solution to files'
writematrix(y, 'coverages_evolution.xlsx') ;
writematrix(t, 'time_evolution.xlsx') ;

%% 'Final coverages'
coverages(:) = y(end,1:length(G_ads_111))' ;

%% 'Make DL arrays'
dl_G_ts_111 = dlarray(G_ts_111) ;
dl_G_ads_111 = dlarray(G_ads_111) ;

dl_G_ads_act_111 = dlarray(G_ads_act_111) ;
dl_G_ads_rxn_111 = dlarray(G_ads_rxn_111) ;

dl_EZPEc_gas_111 = dlarray(EZPEc_gas_111) ;
dl_lnQ_gas_111 = dlarray(lnQ_gas_111) ;

% 'Initialize reaction energies and activation barriers'
dG_act_111 = zeros(length(G_ts_111), 1) ;
dG_rxn_111 = zeros(length(G_ts_111), 1) ;
dG_act_111_LI = zeros(length(G_ts_111), 1) ;
dG_rxn_111_LI = zeros(length(G_ts_111), 1) ;
dl_dG_act_111 = dlarray(dG_act_111) ;
dl_dG_rxn_111 = dlarray(dG_rxn_111) ;
dl_dG_act_111_LI = dlarray(dG_act_111_LI) ;
dl_dG_rxn_111_LI = dlarray(dG_rxn_111_LI) ;
dl_LI_111_H = dlarray(LI_111_H) ;

r111 = zeros(length(G_ts_111),1) ;
Keq_111 = zeros(length(G_ts_111),1) ;
kf_111 = zeros(length(G_ts_111),1) ;
kb_111 = zeros(length(G_ts_111),1) ;

dl_r111 = dlarray(r111) ;
dl_Keq_111 = dlarray(Keq_111) ;
dl_kf_111 = dlarray(kf_111) ;
dl_kb_111 = dlarray(kb_111) ;

no_of_gases = length(EZPEc_gas_111) ;
y_len = length(G_ads_111) ;
no_of_rxns = length(G_ts_111) ;

% 'Initialize rate arrays'
r111_27_plot = zeros(length(y(:,1)),1) ;
r111_23_plot = zeros(length(y(:,1)),1) ;
r111_eth_cons = zeros(length(y(:,1)),1) ;
r111_meth_ads_deso = zeros(length(y(:,1)),1) ;
size(r111_meth_ads_deso)

% 'Initialize dkrc array'
true_coverages_combined = zeros(length(G_ads_111), length(t)) ;
DKRC_E_combined = zeros(length(G_ts_111), length(t)) ;
DTRC_E_combined = zeros(length(G_ads_111), length(t)) ;
DKRC_M_combined = zeros(length(G_ts_111), length(t)) ;
DTRC_M_combined = zeros(length(G_ads_111), length(t)) ;

sigma_E_combined = zeros(length(G_ads_111), length(t)) ;
sigma_M_combined = zeros(length(G_ads_111), length(t)) ;

Sum_DKRC_E = zeros(1, length(t)) ;
Sum_DTRC_E = zeros(1, length(t)) ;
Sum_DKRC_M = zeros(1, length(t)) ;
Sum_DTRC_M = zeros(1, length(t)) ;

Sum_DTRC_DKRC_E = zeros(1, length(t)) ;
Sum_DTRC_DKRC_M = zeros(1, length(t)) ;

%% 'Use ODE/DAE solution to compute DRCs'
for jk = 1:length(t)
    
    current_y_reshape = reshape(y(jk,:), [length(y(end,:)), 1]) ;

    %% 'Add LI to G_ads'
    G_ads_111_LI = zeros(length(G_ads_111), 1) ;
    dl_G_ads_111_LI = dlarray(G_ads_111_LI) ;
    for idx = 1:length(G_ads_111)
        G_ads_111_LI(idx) = G_ads_111(idx) + LI_111_H(idx) * current_y_reshape(16) ;
        dl_G_ads_111_LI(idx) = dl_G_ads_111(idx) + dl_LI_111_H(idx) * current_y_reshape(16) ;
    end
    
    %% 'Get the Reaction energies and Activation Energies'
    %% '111 Reaction energies without LI'
    dG_rxn_111(1) = G_ads_111(12) + G_ads_111(12) - G_ads_111(1) - G_ads_111(17) ;
    dG_rxn_111(2) = G_ads_111(2) + G_ads_111(16) - G_ads_111(1) - G_ads_111(17) ;
    dG_rxn_111(3) = G_ads_111(12) + G_ads_111(13) - G_ads_111(2) - G_ads_111(17) ;
    dG_rxn_111(4) = G_ads_111(3) + G_ads_111(16) - G_ads_111(2) - G_ads_111(17) ;
    dG_rxn_111(5) = G_ads_111(5) + G_ads_111(16) - G_ads_111(2) - G_ads_111(17) ;
    dG_rxn_111(6) = G_ads_111(12) + G_ads_111(14) - G_ads_111(3) - G_ads_111(17) ;
    dG_rxn_111(7) = G_ads_111(4) + G_ads_111(16) - G_ads_111(3) - G_ads_111(17) ;
    dG_rxn_111(8) = G_ads_111(6) + G_ads_111(16) - G_ads_111(3) - G_ads_111(17) ;
    dG_rxn_111(9) = G_ads_111(12) + G_ads_111(15) - G_ads_111(4) - G_ads_111(17) ;
    dG_rxn_111(10) = G_ads_111(7) + G_ads_111(16) - G_ads_111(4) - G_ads_111(17) ;
    dG_rxn_111(11) = G_ads_111(13) + G_ads_111(13) - G_ads_111(5) - G_ads_111(17) ;
    dG_rxn_111(12) = G_ads_111(6) + G_ads_111(16) - G_ads_111(5) - G_ads_111(17) ;
    dG_rxn_111(13) = G_ads_111(13) + G_ads_111(14) - G_ads_111(6) - G_ads_111(17) ;
    dG_rxn_111(14) = G_ads_111(7) + G_ads_111(16) - G_ads_111(6) - G_ads_111(17) ;
    dG_rxn_111(15) = G_ads_111(8) + G_ads_111(16) - G_ads_111(6) - G_ads_111(17) ;
    dG_rxn_111(16) = G_ads_111(13) + G_ads_111(15) - G_ads_111(7) - G_ads_111(17) ;
    dG_rxn_111(17) = G_ads_111(9) + G_ads_111(16) - G_ads_111(7) - G_ads_111(17) ;
    dG_rxn_111(18) = G_ads_111(14) + G_ads_111(14) - G_ads_111(8) - G_ads_111(17) ;
    dG_rxn_111(19) = G_ads_111(9) + G_ads_111(16) - G_ads_111(8) - G_ads_111(17) ;
    dG_rxn_111(20) = G_ads_111(14) + G_ads_111(15) - G_ads_111(9) - G_ads_111(17) ;
    dG_rxn_111(21) = G_ads_111(10) + G_ads_111(16) - G_ads_111(9) - G_ads_111(17) ;
    dG_rxn_111(22) = G_ads_111(15) + G_ads_111(15) - G_ads_111(10) - G_ads_111(17) ;
    dG_rxn_111(23) = G_ads_111(12) + G_ads_111(16) - G_ads_111(11) - G_ads_111(17) ;
    dG_rxn_111(24) = G_ads_111(13) + G_ads_111(16) - G_ads_111(12) - G_ads_111(17) ;
    dG_rxn_111(25) = G_ads_111(14) + G_ads_111(16) - G_ads_111(13) - G_ads_111(17) ;
    dG_rxn_111(26) = G_ads_111(15) + G_ads_111(16) - G_ads_111(14) - G_ads_111(17) ;
    
    %% '111 Reaction energies with LI'
    dG_rxn_111_LI(1) = G_ads_111_LI(12) + G_ads_111_LI(12) - G_ads_111_LI(1) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(2) = G_ads_111_LI(2) + G_ads_111_LI(16) - G_ads_111_LI(1) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(3) = G_ads_111_LI(12) + G_ads_111_LI(13) - G_ads_111_LI(2) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(4) = G_ads_111_LI(3) + G_ads_111_LI(16) - G_ads_111_LI(2) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(5) = G_ads_111_LI(5) + G_ads_111_LI(16) - G_ads_111_LI(2) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(6) = G_ads_111_LI(12) + G_ads_111_LI(14) - G_ads_111_LI(3) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(7) = G_ads_111_LI(4) + G_ads_111_LI(16) - G_ads_111_LI(3) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(8) = G_ads_111_LI(6) + G_ads_111_LI(16) - G_ads_111_LI(3) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(9) = G_ads_111_LI(12) + G_ads_111_LI(15) - G_ads_111_LI(4) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(10) = G_ads_111_LI(7) + G_ads_111_LI(16) - G_ads_111_LI(4) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(11) = G_ads_111_LI(13) + G_ads_111_LI(13) - G_ads_111_LI(5) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(12) = G_ads_111_LI(6) + G_ads_111_LI(16) - G_ads_111_LI(5) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(13) = G_ads_111_LI(13) + G_ads_111_LI(14) - G_ads_111_LI(6) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(14) = G_ads_111_LI(7) + G_ads_111_LI(16) - G_ads_111_LI(6) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(15) = G_ads_111_LI(8) + G_ads_111_LI(16) - G_ads_111_LI(6) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(16) = G_ads_111_LI(13) + G_ads_111_LI(15) - G_ads_111_LI(7) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(17) = G_ads_111_LI(9) + G_ads_111_LI(16) - G_ads_111_LI(7) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(18) = G_ads_111_LI(14) + G_ads_111_LI(14) - G_ads_111_LI(8) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(19) = G_ads_111_LI(9) + G_ads_111_LI(16) - G_ads_111_LI(8) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(20) = G_ads_111_LI(14) + G_ads_111_LI(15) - G_ads_111_LI(9) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(21) = G_ads_111_LI(10) + G_ads_111_LI(16) - G_ads_111_LI(9) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(22) = G_ads_111_LI(15) + G_ads_111_LI(15) - G_ads_111_LI(10) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(23) = G_ads_111_LI(12) + G_ads_111_LI(16) - G_ads_111_LI(11) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(24) = G_ads_111_LI(13) + G_ads_111_LI(16) - G_ads_111_LI(12) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(25) = G_ads_111_LI(14) + G_ads_111_LI(16) - G_ads_111_LI(13) - G_ads_111_LI(17) ;
    dG_rxn_111_LI(26) = G_ads_111_LI(15) + G_ads_111_LI(16) - G_ads_111_LI(14) - G_ads_111_LI(17) ;
    
    %% '111 Activation energies without LI'
    dG_act_111(1) = G_ts_111(1) + G_ads_111(17) - G_ads_111(1) - G_ads_111(17) ;
    dG_act_111(2) = G_ts_111(2) + G_ads_111(17) - G_ads_111(1) - G_ads_111(17) ;
    dG_act_111(3) = G_ts_111(3) + G_ads_111(17) - G_ads_111(2) - G_ads_111(17) ;
    dG_act_111(4) = G_ts_111(4) + G_ads_111(17) - G_ads_111(2) - G_ads_111(17) ;
    dG_act_111(5) = G_ts_111(5) + G_ads_111(17) - G_ads_111(2) - G_ads_111(17) ;
    dG_act_111(6) = G_ts_111(6) + G_ads_111(17) - G_ads_111(3) - G_ads_111(17) ;
    dG_act_111(7) = G_ts_111(7) + G_ads_111(17) - G_ads_111(3) - G_ads_111(17) ;
    dG_act_111(8) = G_ts_111(8) + G_ads_111(17) - G_ads_111(3) - G_ads_111(17) ;
    dG_act_111(9) = G_ts_111(9) + G_ads_111(17) - G_ads_111(4) - G_ads_111(17) ;
    dG_act_111(10) = G_ts_111(10) + G_ads_111(17) - G_ads_111(4) - G_ads_111(17) ;
    dG_act_111(11) = G_ts_111(11) + G_ads_111(17) - G_ads_111(5) - G_ads_111(17) ;
    dG_act_111(12) = G_ts_111(12) + G_ads_111(17) - G_ads_111(5) - G_ads_111(17) ;
    dG_act_111(13) = G_ts_111(13) + G_ads_111(17) - G_ads_111(6) - G_ads_111(17) ;
    dG_act_111(14) = G_ts_111(14) + G_ads_111(17) - G_ads_111(6) - G_ads_111(17) ;
    dG_act_111(15) = G_ts_111(15) + G_ads_111(17) - G_ads_111(6) - G_ads_111(17) ;
    dG_act_111(16) = G_ts_111(16) + G_ads_111(17) - G_ads_111(7) - G_ads_111(17) ;
    dG_act_111(17) = G_ts_111(17) + G_ads_111(17) - G_ads_111(7) - G_ads_111(17) ;
    dG_act_111(18) = G_ts_111(18) + G_ads_111(17) - G_ads_111(8) - G_ads_111(17) ;
    dG_act_111(19) = G_ts_111(19) + G_ads_111(17) - G_ads_111(8) - G_ads_111(17) ;
    dG_act_111(20) = G_ts_111(20) + G_ads_111(17) - G_ads_111(9) - G_ads_111(17) ;
    dG_act_111(21) = G_ts_111(21) + G_ads_111(17) - G_ads_111(9) - G_ads_111(17) ;
    dG_act_111(22) = G_ts_111(22) + G_ads_111(17) - G_ads_111(10) - G_ads_111(17) ;
    dG_act_111(23) = G_ts_111(23) + G_ads_111(17) - G_ads_111(11) - G_ads_111(17) ;
    dG_act_111(24) = G_ts_111(24) + G_ads_111(17) - G_ads_111(12) - G_ads_111(17) ;
    dG_act_111(25) = G_ts_111(25) + G_ads_111(17) - G_ads_111(13) - G_ads_111(17) ;
    dG_act_111(26) = G_ts_111(26) + G_ads_111(17) - G_ads_111(14) - G_ads_111(17) ;
    
    %% '111 Activation energies with LI'
    for i = 1:26
        dG_act_111_LI(i) = dG_act_111(i) + 0.5 * (dG_rxn_111_LI(i) - dG_rxn_111(i)) ;
    end
    
    %% 'dl deltaG reaction without LI'
    dl_dG_rxn_111(1) = dl_G_ads_111(12) + dl_G_ads_111(12) - dl_G_ads_111(1) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(2) = dl_G_ads_111(2) + dl_G_ads_111(16) - dl_G_ads_111(1) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(3) = dl_G_ads_111(12) + dl_G_ads_111(13) - dl_G_ads_111(2) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(4) = dl_G_ads_111(3) + dl_G_ads_111(16) - dl_G_ads_111(2) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(5) = dl_G_ads_111(5) + dl_G_ads_111(16) - dl_G_ads_111(2) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(6) = dl_G_ads_111(12) + dl_G_ads_111(14) - dl_G_ads_111(3) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(7) = dl_G_ads_111(4) + dl_G_ads_111(16) - dl_G_ads_111(3) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(8) = dl_G_ads_111(6) + dl_G_ads_111(16) - dl_G_ads_111(3) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(9) = dl_G_ads_111(12) + dl_G_ads_111(15) - dl_G_ads_111(4) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(10) = dl_G_ads_111(7) + dl_G_ads_111(16) - dl_G_ads_111(4) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(11) = dl_G_ads_111(13) + dl_G_ads_111(13) - dl_G_ads_111(5) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(12) = dl_G_ads_111(6) + dl_G_ads_111(16) - dl_G_ads_111(5) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(13) = dl_G_ads_111(13) + dl_G_ads_111(14) - dl_G_ads_111(6) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(14) = dl_G_ads_111(7) + dl_G_ads_111(16) - dl_G_ads_111(6) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(15) = dl_G_ads_111(8) + dl_G_ads_111(16) - dl_G_ads_111(6) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(16) = dl_G_ads_111(13) + dl_G_ads_111(15) - dl_G_ads_111(7) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(17) = dl_G_ads_111(9) + dl_G_ads_111(16) - dl_G_ads_111(7) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(18) = dl_G_ads_111(14) + dl_G_ads_111(14) - dl_G_ads_111(8) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(19) = dl_G_ads_111(9) + dl_G_ads_111(16) - dl_G_ads_111(8) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(20) = dl_G_ads_111(14) + dl_G_ads_111(15) - dl_G_ads_111(9) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(21) = dl_G_ads_111(10) + dl_G_ads_111(16) - dl_G_ads_111(9) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(22) = dl_G_ads_111(15) + dl_G_ads_111(15) - dl_G_ads_111(10) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(23) = dl_G_ads_111(12) + dl_G_ads_111(16) - dl_G_ads_111(11) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(24) = dl_G_ads_111(13) + dl_G_ads_111(16) - dl_G_ads_111(12) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(25) = dl_G_ads_111(14) + dl_G_ads_111(16) - dl_G_ads_111(13) - dl_G_ads_111(17) ;
    dl_dG_rxn_111(26) = dl_G_ads_111(15) + dl_G_ads_111(16) - dl_G_ads_111(14) - dl_G_ads_111(17) ;
    
    %% 'dl deltaG reaction without LI'
    dl_dG_rxn_111_LI(1) = dl_G_ads_111_LI(12) + dl_G_ads_111_LI(12) - dl_G_ads_111_LI(1) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(2) = dl_G_ads_111_LI(2) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(1) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(3) = dl_G_ads_111_LI(12) + dl_G_ads_111_LI(13) - dl_G_ads_111_LI(2) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(4) = dl_G_ads_111_LI(3) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(2) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(5) = dl_G_ads_111_LI(5) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(2) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(6) = dl_G_ads_111_LI(12) + dl_G_ads_111_LI(14) - dl_G_ads_111_LI(3) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(7) = dl_G_ads_111_LI(4) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(3) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(8) = dl_G_ads_111_LI(6) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(3) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(9) = dl_G_ads_111_LI(12) + dl_G_ads_111_LI(15) - dl_G_ads_111_LI(4) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(10) = dl_G_ads_111_LI(7) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(4) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(11) = dl_G_ads_111_LI(13) + dl_G_ads_111_LI(13) - dl_G_ads_111_LI(5) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(12) = dl_G_ads_111_LI(6) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(5) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(13) = dl_G_ads_111_LI(13) + dl_G_ads_111_LI(14) - dl_G_ads_111_LI(6) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(14) = dl_G_ads_111_LI(7) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(6) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(15) = dl_G_ads_111_LI(8) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(6) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(16) = dl_G_ads_111_LI(13) + dl_G_ads_111_LI(15) - dl_G_ads_111_LI(7) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(17) = dl_G_ads_111_LI(9) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(7) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(18) = dl_G_ads_111_LI(14) + dl_G_ads_111_LI(14) - dl_G_ads_111_LI(8) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(19) = dl_G_ads_111_LI(9) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(8) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(20) = dl_G_ads_111_LI(14) + dl_G_ads_111_LI(15) - dl_G_ads_111_LI(9) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(21) = dl_G_ads_111_LI(10) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(9) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(22) = dl_G_ads_111_LI(15) + dl_G_ads_111_LI(15) - dl_G_ads_111_LI(10) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(23) = dl_G_ads_111_LI(12) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(11) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(24) = dl_G_ads_111_LI(13) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(12) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(25) = dl_G_ads_111_LI(14) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(13) - dl_G_ads_111_LI(17) ;
    dl_dG_rxn_111_LI(26) = dl_G_ads_111_LI(15) + dl_G_ads_111_LI(16) - dl_G_ads_111_LI(14) - dl_G_ads_111_LI(17) ;
    
    %% 'dlarray activation energy'
    dl_dG_act_111(1) = dl_G_ts_111(1) + dl_G_ads_111(17) - dl_G_ads_111(1) - dl_G_ads_111(17) ;
    dl_dG_act_111(2) = dl_G_ts_111(2) + dl_G_ads_111(17) - dl_G_ads_111(1) - dl_G_ads_111(17) ;
    dl_dG_act_111(3) = dl_G_ts_111(3) + dl_G_ads_111(17) - dl_G_ads_111(2) - dl_G_ads_111(17) ;
    dl_dG_act_111(4) = dl_G_ts_111(4) + dl_G_ads_111(17) - dl_G_ads_111(2) - dl_G_ads_111(17) ;
    dl_dG_act_111(5) = dl_G_ts_111(5) + dl_G_ads_111(17) - dl_G_ads_111(2) - dl_G_ads_111(17) ;
    dl_dG_act_111(6) = dl_G_ts_111(6) + dl_G_ads_111(17) - dl_G_ads_111(3) - dl_G_ads_111(17) ;
    dl_dG_act_111(7) = dl_G_ts_111(7) + dl_G_ads_111(17) - dl_G_ads_111(3) - dl_G_ads_111(17) ;
    dl_dG_act_111(8) = dl_G_ts_111(8) + dl_G_ads_111(17) - dl_G_ads_111(3) - dl_G_ads_111(17) ;
    dl_dG_act_111(9) = dl_G_ts_111(9) + dl_G_ads_111(17) - dl_G_ads_111(4) - dl_G_ads_111(17) ;
    dl_dG_act_111(10) = dl_G_ts_111(10) + dl_G_ads_111(17) - dl_G_ads_111(4) - dl_G_ads_111(17) ;
    dl_dG_act_111(11) = dl_G_ts_111(11) + dl_G_ads_111(17) - dl_G_ads_111(5) - dl_G_ads_111(17) ;
    dl_dG_act_111(12) = dl_G_ts_111(12) + dl_G_ads_111(17) - dl_G_ads_111(5) - dl_G_ads_111(17) ;
    dl_dG_act_111(13) = dl_G_ts_111(13) + dl_G_ads_111(17) - dl_G_ads_111(6) - dl_G_ads_111(17) ;
    dl_dG_act_111(14) = dl_G_ts_111(14) + dl_G_ads_111(17) - dl_G_ads_111(6) - dl_G_ads_111(17) ;
    dl_dG_act_111(15) = dl_G_ts_111(15) + dl_G_ads_111(17) - dl_G_ads_111(6) - dl_G_ads_111(17) ;
    dl_dG_act_111(16) = dl_G_ts_111(16) + dl_G_ads_111(17) - dl_G_ads_111(7) - dl_G_ads_111(17) ;
    dl_dG_act_111(17) = dl_G_ts_111(17) + dl_G_ads_111(17) - dl_G_ads_111(7) - dl_G_ads_111(17) ;
    dl_dG_act_111(18) = dl_G_ts_111(18) + dl_G_ads_111(17) - dl_G_ads_111(8) - dl_G_ads_111(17) ;
    dl_dG_act_111(19) = dl_G_ts_111(19) + dl_G_ads_111(17) - dl_G_ads_111(8) - dl_G_ads_111(17) ;
    dl_dG_act_111(20) = dl_G_ts_111(20) + dl_G_ads_111(17) - dl_G_ads_111(9) - dl_G_ads_111(17) ;
    dl_dG_act_111(21) = dl_G_ts_111(21) + dl_G_ads_111(17) - dl_G_ads_111(9) - dl_G_ads_111(17) ;
    dl_dG_act_111(22) = dl_G_ts_111(22) + dl_G_ads_111(17) - dl_G_ads_111(10) - dl_G_ads_111(17) ;
    dl_dG_act_111(23) = dl_G_ts_111(23) + dl_G_ads_111(17) - dl_G_ads_111(11) - dl_G_ads_111(17) ;
    dl_dG_act_111(24) = dl_G_ts_111(24) + dl_G_ads_111(17) - dl_G_ads_111(12) - dl_G_ads_111(17) ;
    dl_dG_act_111(25) = dl_G_ts_111(25) + dl_G_ads_111(17) - dl_G_ads_111(13) - dl_G_ads_111(17) ;
    dl_dG_act_111(26) = dl_G_ts_111(26) + dl_G_ads_111(17) - dl_G_ads_111(14) - dl_G_ads_111(17) ;
    
    %% '111 Activation energies with LI'
    for i = 1:26
        dl_dG_act_111_LI(i) = dl_dG_act_111(i) + 0.5 * (dl_dG_rxn_111_LI(i) - dl_dG_rxn_111(i)) ;
    end
    
    %% 'Gasphase adsorption-desorption reactions: for H2, remember that slab and product energies are squared'
    prod_idx = [1, 5, 8, 11, 16] ;
    for i = 1:length(prod_idx)
        if i == length(prod_idx)
            dG_rxn_111(26 + i) = 2 * G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - 2 * G_ads_111(17) ;
            dl_dG_rxn_111(26 + i) = 2 * dl_G_ads_111(prod_idx(i)) - (dl_EZPEc_gas_111(i) - kb * T * dl_lnQ_gas_111(i)) - 2 * dl_G_ads_111(17) ;
            dG_rxn_111_LI(26 + i) = dG_rxn_111(26 + i) + 2 * ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * current_y_reshape(16)) ;
            dl_dG_rxn_111_LI(26 + i) = dl_dG_rxn_111(26 + i) + 2 * ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * current_y_reshape(16)) ;
            dG_act_111(26 + i) = G_ts_111(26 + i) + G_ads_111(17) - 2 * G_ads_111(17) ;
            dl_dG_act_111(26 + i) = dl_G_ts_111(26 + i) + dl_G_ads_111(17) - 2 * dl_G_ads_111(17) ;
            dG_act_111_LI(26 + i) = dG_act_111(26 + i) + 0.5 * (dG_rxn_111_LI(26 + i) - dG_rxn_111(26 + i)) ;
            dl_dG_act_111(26 + i) = dl_dG_act_111(26 + i) + 0.5 * (dl_dG_rxn_111_LI(26 + i) - dl_dG_rxn_111(26 + i)) ;
        else
            dG_rxn_111(26 + i) = G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - G_ads_111(17) ;
            dl_dG_rxn_111(26 + i) = dl_G_ads_111(prod_idx(i)) - (dl_EZPEc_gas_111(i) - kb * T * dl_lnQ_gas_111(i)) - dl_G_ads_111(17) ;
            dG_rxn_111_LI(26 + i) = dG_rxn_111(26 + i) + ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * current_y_reshape(16)) ;
            dl_dG_rxn_111_LI(26 + i) = dl_dG_rxn_111(26 + i) + ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * current_y_reshape(16)) ;
            dG_act_111(26 + i) = G_ts_111(26 + i) - G_ads_111(17) ;
            dl_dG_act_111(26 + i) = dl_G_ts_111(26 + i) - dl_G_ads_111(17) ;
            dG_act_111_LI(26 + i) = dG_act_111(26 + i) + 0.5 * (dG_rxn_111_LI(26 + i) - dG_rxn_111(26 + i)) ;
            dl_dG_act_111_LI(26 + i) = dl_dG_act_111(26 + i) + 0.5 * (dl_dG_rxn_111_LI(26 + i) - dl_dG_rxn_111(26 + i)) ;
        end
    end
    
    %% 'Get Reaction constants'
    for idx = 1:length(G_ts_111)
        if (dG_act_111_LI(idx) < 0)
            dG_act_111(idx) = 0.001 ;
        end
    
        if (dG_act_111_LI(idx) < dG_rxn_111_LI(idx))
		    dG_act_111_LI(idx) = dG_rxn_111_LI(idx) ; 
        end
    
        Keq_111(idx) = exp(-(dG_rxn_111_LI(idx)) / (kb * T));
        kf_111(idx) = (kb * T / h) * exp(-(dG_act_111_LI(idx)) / (kb * T)) ;
        kb_111(idx) = kf_111(idx) / Keq_111(idx) ;
    
        dl_Keq_111(idx) = exp(-(dl_dG_rxn_111_LI(idx)) / (kb * T)) ;
        dl_kf_111(idx) = (kb * T / h) * exp(-(dl_dG_act_111_LI(idx)) / (kb * T)) ;
        dl_kb_111(idx) = dl_kf_111(idx) / dl_Keq_111(idx) ;
    end
    
    % 'Add the evolution of a reaction step to the plots'
    r111_27_plot(jk) = kf_111(27)*x(1)*Ptot*current_y_reshape(17) - kb_111(27)*current_y_reshape(1) ;
    r111_23_plot(jk) = -(kf_111(23)*current_y_reshape(11)*current_y_reshape(17) - kb_111(23)*current_y_reshape(12)*current_y_reshape(16)) ;
    r111_eth_cons(jk) = kf_111(1)*current_y_reshape(1)*current_y_reshape(17) - kb_111(1)*current_y_reshape(12)*current_y_reshape(12) + kf_111(2)*current_y_reshape(1)*current_y_reshape(17) - kb_111(2)*current_y_reshape(2)*current_y_reshape(16) ;
    r111_meth_ads_deso(jk) = -(kf_111(30)*x(4)*Ptot*current_y_reshape(17) - kb_111(30)*current_y_reshape(11)) ;

    %% 'Extract Q_Gts, Q_Gads and R matrices from y solution and reshape'
    Q_Gts = reshape(current_y_reshape(length(x0)+1:length(x0)+length(x0)*length(G_ts_111)), [length(x0), length(G_ts_111)]) ;
    Q_Gads = reshape(current_y_reshape(length(x0)+1+length(x0)*length(G_ts_111):length(x0)+length(x0)*length(G_ts_111)+length(x0)*length(G_ads_111)), [length(x0), length(G_ads_111)]) ;
    R = reshape(current_y_reshape(length(x0)+1+length(x0)*length(G_ts_111)+length(x0)*length(G_ads_111):end), [length(x0),length(x0)]) ;

    dl_y = dlarray(current_y_reshape) ;
    dl_x = dlarray(x) ;

    site_counts = [1, 1, 2, 3, 2, 3, 3, 3, 3, 4, 1, 1, 2, 3, 3, 1, 1] ;
    y_site_balance = site_counts * current_y_reshape(1:length(x0)) ;
    true_coverages = site_counts' .* current_y_reshape(1:length(x0)) ;
    true_coverages_combined(:, jk) = true_coverages ;

    sum(x(1:5)) % 'Sum of Gas mole fractions'
    
    PCH3CH3_f = Ptot * x(1) ;
    PCH2CH2_f = Ptot * x(2) ;
    PCHCH_f = Ptot * x(3) ;
    PCH4_f = Ptot * x(4) ;
    PH2_f = Ptot * x(5) ;

    Ptot_f = PCH3CH3_f + PCH2CH2_f + PCHCH_f + PCH4_f + PH2_f ;

    %% 'Compute instantaneous differentials of lnTOF w.r.t G_ts, G_ads and coverages at current time-coverages'
    [lnTOF_Ethane, lnTOF_Methane, dlnTOF_E_dGts_thetaconst, dlnTOF_M_dGts_thetaconst, dlnTOF_E_dGads_thetaconst, dlnTOF_M_dGads_thetaconst, dlnTOF_E_dtheta_Gconst,  dlnTOF_M_dtheta_Gconst] = dlfeval(@compute_DKRC_DTRC_ThetaConstant_DL, current_y_reshape(1:length(x0)), dl_y(1:length(x0)), x, dl_x, Ptot_f, G_ts_111, dl_G_ts_111, G_ads_111, dl_G_ads_111, EZPEc_gas_111, dl_EZPEc_gas_111, lnQ_gas_111, dl_lnQ_gas_111, LI_111_H, dl_LI_111_H) ;

    %% 'Compute DTRC and DKRC'
    % 'Ethane'
    n_dkrc_E = -1 * kb * T * (dlnTOF_E_dGts_thetaconst + Q_Gts.' * dlnTOF_E_dtheta_Gconst) ;
    n_dtrc_E = -1 * kb * T * (dlnTOF_E_dGads_thetaconst + Q_Gads.' * dlnTOF_E_dtheta_Gconst) ;
    sigma_E = -1 * n_dtrc_E ./ true_coverages ;

    % 'Methane'
    n_dkrc_M = -1 * kb * T * (dlnTOF_M_dGts_thetaconst + Q_Gts.' * dlnTOF_M_dtheta_Gconst) ;
    n_dtrc_M = -1 * kb * T * (dlnTOF_M_dGads_thetaconst + Q_Gads.' * dlnTOF_M_dtheta_Gconst) ;
    sigma_M = -1 * n_dtrc_M ./ true_coverages ;

    DKRC_E_combined(:, jk) = extractdata(n_dkrc_E) ;
    DTRC_E_combined(:, jk) = extractdata(n_dtrc_E) ;
    DKRC_M_combined(:, jk) = extractdata(n_dkrc_M) ;
    DTRC_M_combined(:, jk) = extractdata(n_dtrc_M) ;
    sigma_E_combined(:, jk) = extractdata(sigma_E) ;
    sigma_M_combined(:, jk) = extractdata(sigma_M) ;

    Sum_DTRC_E(jk) = sum(DTRC_E_combined(:, jk)) ;
    Sum_DKRC_E(jk) = sum(DKRC_E_combined(:, jk)) ;
    Sum_DTRC_M(jk) = sum(DTRC_M_combined(:, jk)) ;
    Sum_DKRC_M(jk) = sum(DKRC_M_combined(:, jk)) ;

    Sum_DTRC_DKRC_E(jk) = Sum_DTRC_E(jk) + Sum_DKRC_E(jk) ;
    Sum_DTRC_DKRC_M(jk) = Sum_DTRC_M(jk) + Sum_DKRC_M(jk) ;
end

%% 'Plot of coverage solution'
loglog(t,r111_27_plot(:,1),'r', t,r111_23_plot(:,1),'k', t,r111_eth_cons(:,1),'m', t,r111_meth_ads_deso(:,1),'m',t, y(:,1),'b',t,y(:,2),'b',t,y(:,3),'b',t,y(:,4),'b',t,y(:,5),'b',t,y(:,6),'b',t,y(:,7),'b',t,y(:,8),'b',t,y(:,9),'b',t,y(:,10),'b',t,y(:,11),'b',t,y(:,12),'b',t,y(:,13),'b',t,y(:,14),'b',t,y(:,15),'b',t,y(:,16),'b',t,y(:,17),'k') ;
title('Solution of balance Equation');
xlabel('time t'); ylabel('solution y');
l = legend('r111_{27}','r111_{23}', 'r111_{ethcons}', 'r111_{meth_deso}','y 1','y 2','y 3','y 4','y 5','y 6','y 7','y 8','y 9','y 10','y 11','y 12','y 13','y 14','y 15','y 16','y 17') ;
set(l,'Edgecolor',[1 1 1]);

y_raw = y(end,:);
y = reshape(y_raw, [length(y_raw), 1]) ;

y ;
Ptot_f = PCH3CH3_f + PCH2CH2_f + PCHCH_f + PCH4_f + PH2_f ;

% 'Overall Rate Expressions'
r111(1) = kf_111(1)*y(1)*y(17) - kb_111(1)*y(12)*y(12) ;
r111(2) = kf_111(2)*y(1)*y(17) - kb_111(2)*y(2)*y(16) ;
r111(3) = kf_111(3)*y(2)*y(17)^2 - kb_111(3)*y(12)*y(13) ;
r111(4) = kf_111(4)*y(2)*y(17)^2 - kb_111(4)*y(3)*y(16) ;
r111(5) = kf_111(5)*y(2)*y(17)^2 - kb_111(5)*y(5)*y(16) ;
r111(6) = kf_111(6)*y(3)*y(17)^2 - kb_111(6)*y(12)*y(14) ;
r111(7) = kf_111(7)*y(3)*y(17)^2 - kb_111(7)*y(4)*y(16) ;
r111(8) = kf_111(8)*y(3)*y(17)^2 - kb_111(8)*y(6)*y(16) ;
r111(9) = kf_111(9)*y(4)*y(17) - kb_111(9)*y(12)*y(15) ;
r111(10) = kf_111(10)*y(4)*y(17) - kb_111(10)*y(7)*y(16) ;
r111(11) = kf_111(11)*y(5)*y(17)^2 - kb_111(11)*y(13)*y(13) ;
r111(12) = kf_111(12)*y(5)*y(17)^2 - kb_111(12)*y(6)*y(16) ;
r111(13) = kf_111(13)*y(6)*y(17)^2 - kb_111(13)*y(13)*y(14) ;
r111(14) = kf_111(14)*y(6)*y(17) - kb_111(14)*y(7)*y(16) ;
r111(15) = kf_111(15)*y(6)*y(17) - kb_111(15)*y(8)*y(16) ;
r111(16) = kf_111(16)*y(7)*y(17)^2 - kb_111(16)*y(13)*y(15) ;
r111(17) = kf_111(17)*y(7)*y(17) - kb_111(17)*y(9)*y(16) ;
r111(18) = kf_111(18)*y(8)*y(17)^3 - kb_111(18)*y(14)*y(14) ;
r111(19) = kf_111(19)*y(8)*y(17) - kb_111(19)*y(9)*y(16) ;
r111(20) = kf_111(20)*y(9)*y(17)^3 - kb_111(20)*y(14)*y(15) ;
r111(21) = kf_111(21)*y(9)*y(17)^2 - kb_111(21)*y(10)*y(16) ;
r111(22) = kf_111(22)*y(10)*y(17)^2 - kb_111(22)*y(15)*y(15) ;
r111(23) = kf_111(23)*y(11)*y(17) - kb_111(23)*y(12)*y(16) ;
r111(24) = kf_111(24)*y(12)*y(17)^2 - kb_111(24)*y(13)*y(16) ;
r111(25) = kf_111(25)*y(13)*y(17)^2 - kb_111(25)*y(14)*y(16) ;
r111(26) = kf_111(26)*y(14)*y(17) - kb_111(26)*y(15)*y(16) ;
r111(27) = kf_111(27)*x(1)*Ptot*y(17) - kb_111(27)*y(1) ;
r111(28) = kf_111(28)*x(2)*Ptot*y(17)^2 - kb_111(28)*y(5) ;
r111(29) = kf_111(29)*x(3)*Ptot*y(17)^3 - kb_111(29)*y(8) ;
r111(30) = kf_111(30)*x(4)*Ptot*y(17) - kb_111(30)*y(11) ;
r111(31) = kf_111(31)*x(5)*Ptot*y(17)^2 - kb_111(31)*y(16)^2 ;

% 'Forward Rates Expressions'
f111(1) = kf_111(1)*y(1)*y(17) ;
f111(2) = kf_111(2)*y(1)*y(17) ;
f111(3) = kf_111(3)*y(2)*y(17)^2 ;
f111(4) = kf_111(4)*y(2)*y(17)^2 ;
f111(5) = kf_111(5)*y(2)*y(17)^2 ;
f111(6) = kf_111(6)*y(3)*y(17)^2 ;
f111(7) = kf_111(7)*y(3)*y(17)^2 ;
f111(8) = kf_111(8)*y(3)*y(17)^2 ;
f111(9) = kf_111(9)*y(4)*y(17) ;
f111(10) = kf_111(10)*y(4)*y(17) ;
f111(11) = kf_111(11)*y(5)*y(17)^2 ;
f111(12) = kf_111(12)*y(5)*y(17)^2 ;
f111(13) = kf_111(13)*y(6)*y(17)^2 ;
f111(14) = kf_111(14)*y(6)*y(17) ;
f111(15) = kf_111(15)*y(6)*y(17) ;
f111(16) = kf_111(16)*y(7)*y(17)^2 ;
f111(17) = kf_111(17)*y(7)*y(17) ;
f111(18) = kf_111(18)*y(8)*y(17)^3 ;
f111(19) = kf_111(19)*y(8)*y(17) ;
f111(20) = kf_111(20)*y(9)*y(17)^3 ;
f111(21) = kf_111(21)*y(9)*y(17)^2 ;
f111(22) = kf_111(22)*y(10)*y(17)^2 ;
f111(23) = kf_111(23)*y(11)*y(17) ;
f111(24) = kf_111(24)*y(12)*y(17)^2 ;
f111(25) = kf_111(25)*y(13)*y(17)^2 ;
f111(26) = kf_111(26)*y(14)*y(17) ;
f111(27) = kf_111(27)*x(1)*Ptot*y(17) ;
f111(28) = kf_111(28)*x(2)*Ptot*y(17)^2 ;
f111(29) = kf_111(29)*x(3)*Ptot*y(17)^3 ;
f111(30) = kf_111(30)*x(4)*Ptot*y(17) ;
f111(31) = kf_111(31)*x(5)*Ptot*y(17)^2 ;

% 'Backward Rate Expressions'
b111(1) = kb_111(1)*y(12)*y(12) ;
b111(2) = kb_111(2)*y(2)*y(16) ;
b111(3) = kb_111(3)*y(12)*y(13) ;
b111(4) = kb_111(4)*y(3)*y(16) ;
b111(5) = kb_111(5)*y(5)*y(16) ;
b111(6) = kb_111(6)*y(12)*y(14) ;
b111(7) = kb_111(7)*y(4)*y(16) ;
b111(8) = kb_111(8)*y(6)*y(16) ;
b111(9) = kb_111(9)*y(12)*y(15) ;
b111(10) = kb_111(10)*y(7)*y(16) ;
b111(11) = kb_111(11)*y(13)*y(13) ;
b111(12) = kb_111(12)*y(6)*y(16) ;
b111(13) = kb_111(13)*y(13)*y(14) ;
b111(14) = kb_111(14)*y(7)*y(16) ;
b111(15) = kb_111(15)*y(8)*y(16) ;
b111(16) = kb_111(16)*y(13)*y(15) ;
b111(17) = kb_111(17)*y(9)*y(16) ;
b111(18) = kb_111(18)*y(14)*y(14) ;
b111(19) = kb_111(19)*y(9)*y(16) ;
b111(20) = kb_111(20)*y(14)*y(15) ;
b111(21) = kb_111(21)*y(10)*y(16) ;
b111(22) = kb_111(22)*y(15)*y(15) ;
b111(23) = kb_111(23)*y(12)*y(16) ;
b111(24) = kb_111(24)*y(13)*y(16) ;
b111(25) = kb_111(25)*y(14)*y(16) ;
b111(26) = kb_111(26)*y(15)*y(16) ;
b111(27) = kb_111(27)*y(1) ;
b111(28) = kb_111(28)*y(5) ;
b111(29) = kb_111(29)*y(8) ;
b111(30) = kb_111(30)*y(11) ;
b111(31) = kb_111(31)*y(16)^2 ;

r=r111(1:31);
f=f111(1:31);
b=b111(1:31);

%% 'Write results to files'

% 'Rates'
rate_combined = r111 ;
fileID = fopen('coverage.csv','a');
fmt = '%5d\r\n';
fprintf(fileID,fmt,y);
fclose(fileID);
rate = fopen('rate.csv','a');
frate= fopen('frate.csv','a');
brate= fopen('brate.csv','a');
fmt = '%5d\r\n';
fprintf(rate,fmt,r);
fprintf(frate,fmt,f);
fprintf(brate,fmt,b);
fclose(rate);
fclose(frate);
fclose(brate);

Ethane_conv_rate=fopen('TOF','a');
Ethane_conv_TOF = r111(1) + r111(2) ;
Methane_prod_TOF = -r111(23) ;
fprintf(Ethane_conv_rate,'%22s %5d %20s %5d %12s %6.3f %12s %6.3f\r\n','Ethane conversion= ', Ethane_conv_TOF,'Methane production= ', Methane_prod_TOF, 'PCH3CH3= ',PCH3CH3,'PH2= ',PH2);
fclose(Ethane_conv_rate);

% 'Write site balance'
sites_bal=fopen('sites', 'a');
fprintf(sites_bal, '%6s %6.3f\r\n', 'Site Balance', y_site_balance);
fclose(sites_bal);

writematrix(r111_27_plot, 'ethane_ads_deso.xlsx') ;
writematrix(r111_23_plot, 'methane_production.xlsx') ;
writematrix(r111_eth_cons, 'ethane_consumption.xlsx') ;
writematrix(r111_meth_ads_deso, 'methane_ads_deso.xlsx')
writematrix(coverages, 'combined_coverages.xls');
writematrix(true_coverages_combined, 'true_coverages.csv');
writematrix(rate_combined, 'combined_rate.xls');
writematrix(DKRC_E_combined, 'combined_DKRC_E.csv') ;
writematrix(DTRC_E_combined, 'combined_DTRC_E.csv') ;
writematrix(DKRC_M_combined, 'combined_DKRC_M.csv') ;
writematrix(DTRC_M_combined, 'combined_DTRC_M.csv') ;
writematrix(sigma_E_combined, 'combined_sigma_E.csv') ;
writematrix(sigma_M_combined, 'combined_sigma_M.csv') ;
writematrix(J, 'jacobian.csv') ;
writematrix(extractdata(dlnTOF_E_dtheta_Gconst), 'dlnTOF_E_dtheta_Gconst.csv') ;
writematrix(extractdata(dlnTOF_E_dGts_thetaconst), 'dlnTOF_E_dGts_thetaconst.csv') ;
writematrix(extractdata(dlnTOF_E_dGads_thetaconst), 'dlnTOF_E_dGads_thetaconst.csv') ;
writematrix(extractdata(dlnTOF_M_dtheta_Gconst), 'dlnTOF_M_dtheta_Gconst.csv') ;
writematrix(extractdata(dlnTOF_M_dGts_thetaconst), 'dlnTOF_M_dGts_thetaconst.csv') ;
writematrix(extractdata(dlnTOF_M_dGads_thetaconst), 'dlnTOF_M_dGads_thetaconst.csv') ;
writematrix(Q_Gts, 'dthetadGts.csv') ;
writematrix(Q_Gads, 'dthetadGads.csv') ;
writematrix(R, 'dthetadtheta0.csv') ;
writematrix(t, 'time_evolution.csv')
writematrix(Sum_DTRC_DKRC_M, 'sum_DTRC_DKRC_M.csv')
writematrix(Sum_DTRC_DKRC_E, 'sum_DTRC_DKRC_E.csv')
writematrix(Sum_DTRC_M, 'sum_DTRC_M.csv')
writematrix(Sum_DKRC_M, 'sum_DKRC_M.csv')
writematrix(Sum_DTRC_E, 'sum_DTRC_E.csv')
writematrix(Sum_DKRC_E, 'sum_DKRC_E.csv')

%%%%%%%%%%%% Inline Functions %%%%%%%%%%%%%%%%
%% 'Read DFT Energies and Parameters supplied in Excel Worksheet'
function [G_ads, G_ts, G_ads_act, G_ads_rxn, EZPEc_gas, lnQ_gas, LI_params] = read_Energies(T, species_file, ts_file, facet)

    % 'TS energy sheet name in excel file'
    ts_sheet = 'ts_energy_' + string(facet) ;
    
    G_ads_raw = readtable(species_file,'ReadRowNames',true) ;
    G_ts_raw = readtable(ts_file,'Sheet',ts_sheet,'ReadRowNames',true) ;
    G_ads_act_raw = readtable(ts_file,'Sheet','ads_act','ReadRowNames',true) ;
    G_ads_rxn_raw = readtable(ts_file,'Sheet','ads_rxn','ReadRowNames',true) ;
    EZPEc_gas_raw = readtable(ts_file,'Sheet','EZPEc_gas','ReadRowNames',true) ;
    lnQ_gas_raw = readtable(ts_file,'Sheet','lnQ_gas','ReadRowNames',true) ;
    LI_params_raw = readtable(ts_file,'Sheet','LI_params','ReadRowNames',true) ;
    
    T_range = linspace(473.15, 883.15, 42) ;
    if ismember(T, T_range)
        % 'Get the G(s) at that temperature'
        T_modified = T * 100 ;
        T_readable = 'T' + string(T_modified) + '_K' ;
        G_ads = G_ads_raw.(T_readable) ;
        G_ts = G_ts_raw.(T_readable) ;
        G_ads_act = G_ads_act_raw.(T_readable) ;
        G_ads_rxn = G_ads_rxn_raw.(T_readable) ;
        lnQ_gas = lnQ_gas_raw.(T_readable) ;
    else
        % 'Get the temp bounds where operating temperature falls in the array of T-range'
        lower_T_bound = max(T_range(T >= T_range));
        upper_T_bound = min(T_range(T <= T_range));
        lower_T_readable = 'T' + string(lower_T_bound * 100) + '_K';
        upper_T_readable = 'T' + string(upper_T_bound * 100) + '_K';
        
        % 'Get intermediate energies'
        G_ads_lower = G_ads_raw.(lower_T_readable) ;
        G_ads_upper = G_ads_raw.(upper_T_readable) ;
        
        G_ts_lower = G_ts_raw.(lower_T_readable) ;
        G_ts_upper = G_ts_raw.(upper_T_readable) ;
        
        G_ads_act_lower = G_ads_act_raw.(lower_T_readable) ;
        G_ads_act_upper = G_ads_act_raw.(upper_T_readable) ;
        
        G_ads_rxn_lower = G_ads_rxn_raw.(lower_T_readable) ;
        G_ads_rxn_upper = G_ads_rxn_raw.(upper_T_readable) ;
    
        lnQ_gas_lower = lnQ_gas_raw.(lower_T_readable) ;
        lnQ_gas_upper = lnQ_gas_raw.(upper_T_readable) ;
        
        delta_H_ads = (((G_ads_upper/upper_T_bound) - (G_ads_lower/lower_T_bound))/((1/upper_T_bound) - (1/lower_T_bound))) ;
        delta_H_ts = (((G_ts_upper/upper_T_bound) - (G_ts_lower/lower_T_bound))/((1/upper_T_bound) - (1/lower_T_bound))) ;
        delta_H_ads_act = (((G_ads_act_upper/upper_T_bound) - (G_ads_act_lower/lower_T_bound))/((1/upper_T_bound) - (1/lower_T_bound))) ;
        delta_H_ads_rxn = (((G_ads_rxn_upper/upper_T_bound) - (G_ads_rxn_lower/lower_T_bound))/((1/upper_T_bound) - (1/lower_T_bound))) ;
        delta_H_lnQ_gas = (((lnQ_gas_upper/upper_T_bound) - (lnQ_gas_lower/lower_T_bound))/((1/upper_T_bound) - (1/lower_T_bound))) ;
        
        G_ads = T * ((delta_H_ads * ((1/T) - (1/lower_T_bound))) + (G_ads_lower/lower_T_bound));
        G_ts = T * ((delta_H_ts * ((1/T) - (1/lower_T_bound))) + (G_ts_lower/lower_T_bound));
        G_ads_act = T * ((delta_H_ads_act * ((1/T) - (1/lower_T_bound))) + (G_ads_act_lower/lower_T_bound));
        G_ads_rxn = T * ((delta_H_ads_rxn * ((1/T) - (1/lower_T_bound))) + (G_ads_rxn_lower/lower_T_bound));
        lnQ_gas = T * ((delta_H_lnQ_gas * ((1/T) - (1/lower_T_bound))) + (lnQ_gas_lower/lower_T_bound));
    end
    
    EZPEc_gas = EZPEc_gas_raw.('EZPEc_gas') ;
    LI_params = LI_params_raw ;

end

%% 'Activate or Deactivate Lateral Interactions'
function LI_arr_new = ActivateLI(LI_arr, turn_on)
    if turn_on == "yes"
        % 'Return original array'
        LI_arr_new = LI_arr ;
    elseif turn_on == "no"
        % 'Return an array of zeros'
        LI_arr_new = zeros(length(LI_arr), 1) ;
    end
end
