clear
close all
format long

%% 'Delete previously existing result files if they exist'
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
delete sum_DTRC_M.csv
delete sum_DKRC_M.csv
delete sum_DTRC_E.csv
delete sum_DKRC_E.csv

%% 'Add casadi to path'
addpath('/home/x-bamidele/per_bin/casadi-3.6.7-linux64-matlab2018b/')

%% 'run_time'
% 'set time limit here'
run_time = 1.60e-0 ;

%% 'Constants'
kb=8.617333262e-5 ; % 'Boltzmanns constant'
h=4.135667696e-15 ; % 'Plancks constant'

% 'Operating conditions'
T = 573.15 ;

%% 'File names'
species_file = 'species_ads_111.csv' ;
ts_file = 'ts_energy_111.xlsx' ;

%% 'Read Adsorption energies'
G_ads_111_raw = readtable(species_file,'ReadRowNames',true) ;
G_ts_111_raw = readtable(ts_file,'Sheet','ts_energy_111','ReadRowNames',true) ;
G_ads_act_111_raw = readtable(ts_file,'Sheet','ads_act','ReadRowNames',true) ;
G_ads_rxn_111_raw = readtable(ts_file,'Sheet','ads_rxn','ReadRowNames',true) ;
EZPEc_gas_111_raw = readtable(ts_file,'Sheet','EZPEc_gas','ReadRowNames',true) ;
lnQ_gas_111_raw = readtable(ts_file,'Sheet','lnQ_gas','ReadRowNames',true) ;

% 'Get the G(s) at that temperature'
T_modified = T * 100 ;
T_readable = 'T' + string(T_modified) + '_K' ;
G_ads_111 = G_ads_111_raw.(T_readable) ;
G_ts_111 = G_ts_111_raw.(T_readable) ;
G_ads_act_111 = G_ads_act_111_raw.(T_readable) ;
G_ads_rxn_111 = G_ads_rxn_111_raw.(T_readable) ;
lnQ_gas_111 = lnQ_gas_111_raw.(T_readable) ;

EZPEc_gas_111 = EZPEc_gas_111_raw.('EZPEc_gas') ;

%% 'Add adsorption-desorption steps activation energy to G_ts_111'
G_ts_111(27:31) = G_ads_act_111 ;

%% 'Make DeepLearning arrays of read energies'
dl_G_ts_111 = dlarray(G_ts_111) ;
dl_G_ads_111 = dlarray(G_ads_111) ;
dl_G_ads_act_111 = dlarray(G_ads_act_111) ;
dl_G_ads_rxn_111 = dlarray(G_ads_rxn_111) ;
dl_EZPEc_gas_111 = dlarray(EZPEc_gas_111) ;
dl_lnQ_gas_111 = dlarray(lnQ_gas_111) ;

%% 'Initialize reaction energies, activation barriers rates and constants'
dG_act_111 = zeros(length(G_ts_111), 1) ;
dG_rxn_111 = zeros(length(G_ts_111), 1) ;
r111 = zeros(length(G_ts_111),1) ;
Keq_111 = zeros(length(G_ts_111),1) ;
kf_111 = zeros(length(G_ts_111),1) ;
kb_111 = zeros(length(G_ts_111),1) ;

%% 'Make the corresponding DeepLearning arrays'
dl_dG_act_111 = dlarray(dG_act_111) ;
dl_dG_rxn_111 = dlarray(dG_rxn_111) ;
dl_r111 = dlarray(r111) ;
dl_Keq_111 = dlarray(Keq_111) ;
dl_kf_111 = dlarray(kf_111) ;
dl_kb_111 = dlarray(kb_111) ;

%% 'Get the Reaction energies and Activation Energies'
% '111 Activation energies of TS steps'
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

%% 'dlarray activation energy of TS steps'
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

%% '111 Reaction energies of TS steps'
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

%% 'dlarray deltaG reaction of TS steps'
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


%% 'Get Reaction constants'
for idx = 1:26
    if (dG_act_111(idx) < 0)
        dl_dG_act_111(idx) = 0.001 ;
        dG_act_111(idx) = 0.001 ;
    end
    Keq_111(idx) = exp(-(dG_rxn_111(idx)) / (kb * T));
    kf_111(idx) = (kb * T / h) * exp(-(dG_act_111(idx)) / (kb * T)) ;
    kb_111(idx) = kf_111(idx) / Keq_111(idx) ;

    dl_Keq_111(idx) = exp(-(dl_dG_rxn_111(idx)) / (kb * T)) ;
    dl_kf_111(idx) = (kb * T / h) * exp(-(dl_dG_act_111(idx)) / (kb * T)) ;
    dl_kb_111(idx) = dl_kf_111(idx) / dl_Keq_111(idx) ;
end

% 'Activation and reaction energies of gas phase adsorption desorption steps'
% 'for H2 (index 16), remember that slab and product energy is squared'
prod_idx = [1, 5, 8, 11, 16] ;
for i = 1:length(prod_idx)
    if i == length(prod_idx)
	% "Dissociative adsorption of H2 to 2H(*)"
        dG_act_111(26 + i) = G_ts_111(26 + i) + G_ads_111(17) - 2 * G_ads_111(17) ;
        dl_dG_act_111(26 + i) = dl_G_ts_111(26 + i) + dl_G_ads_111(17) - 2 * dl_G_ads_111(17) ;

        dG_rxn_111(26 + i) = 2 * G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - 2 * G_ads_111(17) ;
        dl_dG_rxn_111(26 + i) = 2 * dl_G_ads_111(prod_idx(i)) - (dl_EZPEc_gas_111(i) - kb * T * dl_lnQ_gas_111(i)) - 2 * dl_G_ads_111(17) ;
    else
	% 'Other gas phase'
        dG_act_111(26 + i) = G_ts_111(26 + i) - G_ads_111(17) ;
        dl_dG_act_111(26 + i) = dl_G_ts_111(26 + i) - dl_G_ads_111(17) ;

        dG_rxn_111(26 + i) = G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - G_ads_111(17) ;
        dl_dG_rxn_111(26 + i) = dl_G_ads_111(prod_idx(i)) - (dl_EZPEc_gas_111(i) - kb * T * dl_lnQ_gas_111(i)) - dl_G_ads_111(17) ;
    end

    Keq_111(26 + i) = exp(-(dG_rxn_111(26 + i)) / (kb * T)) ;
    dl_Keq_111(26 + i) = exp(-(dl_dG_rxn_111(26 + i)) / (kb * T)) ;

    kf_111(26 + i) = (kb * T / h) * exp(-(dG_act_111(26 + i)) / (kb * T)) ;
    kb_111(26 + i) = kf_111(26 + i) / Keq_111(26 + i) ;

    dl_kf_111(26 + i) = (kb * T / h) * exp(-(dl_dG_act_111(26 + i)) / (kb * T)) ;
    dl_kb_111(26 + i) = dl_kf_111(26 + i) / dl_Keq_111(26 + i) ;
end

H2_pressure = 0.1317225 ; % "0.13 atm to bar"
Ethane_pressure = 0.03343725 ; % "0.033 atm to bar"
coverages = zeros(length(G_ads_111), 1) ;
rate_combined = zeros(length(G_ts_111), 1) ;

PCH3CH3 = Ethane_pressure;
PCH2CH2 = 0;
PCHCH = 0;
PCH4 = 0;
PH2 = H2_pressure;

Ptot = PCH3CH3 + PCH2CH2 + PCHCH + PCH4 + PH2 ;

% 'Initialize theta: state variables'
x0=zeros(length(G_ads_111),1);

% 'From clean surface, 100 percent free sites'
x0(17) = 1;

% 'Gas mole fraction'
x(1) = PCH3CH3 / Ptot ;
x(2) = PCH2CH2 / Ptot ;
x(3) = PCHCH / Ptot ;
x(4) = PCH4 / Ptot ;
x(5) = PH2 / Ptot ;

% 'Initialize Q_Gts0, Q_Gads0 and R0 and columnize them'
Q_Gts0 = zeros(length(G_ts_111)*length(x0), 1) ;
Q_Gads0 = zeros(length(G_ads_111)*length(x0), 1) ;
R0 = reshape(eye(length(x0)), [length(x0)*length(x0),1]) ;

%% 'Concatenate theta0, and columnized Q_Gts0, Q_Gads0 and R0'
y0 = [x0 ; Q_Gts0; Q_Gads0; R0] ;

% 'Mass matrix for implicit ODE setup'
M=eye(length(y0)) ;

% 'Algebraic equation for free sites balance'
M(length(x0),length(x0)) = 0 ;

optode1 = odeset('NonNegative',1:17, 'Mass',M,'Abstol',1E-15,'RelTol',1E-15) ;
optode = odeset(optode1,'Abstol',1E-10,'RelTol',1E-10) ;
optlsq = optimset('TolFun',1E-13,'TolX',1E-13) ;      

%% 'Solve the composite DAE'
[t,y]=ode15s(@(t,y)C2_Pt_111_Composite_ODE_function(t, y, x, T, Ptot, G_ts_111, G_ads_111, EZPEc_gas_111, lnQ_gas_111),[0, run_time], y0, optode);

% 'Final coverages, theta'
coverages(:) = y(end,1:17);

%% 'Initialize rate control arrays'
% 'Ethane'
DKRC_E_combined = zeros(length(G_ts_111), length(t)) ;
DTRC_E_combined = zeros(length(G_ads_111), length(t)) ;
% 'Methane'
DKRC_M_combined = zeros(length(G_ts_111), length(t)) ;
DTRC_M_combined = zeros(length(G_ads_111), length(t)) ;

%% 'Sum of DRCs'
% 'Ethane'
Sum_DKRC_E = zeros(1, length(t)) ;
Sum_DTRC_E = zeros(1, length(t)) ;
% 'Methane'
Sum_DKRC_M = zeros(1, length(t)) ;
Sum_DTRC_M = zeros(1, length(t)) ;

%% 'Sum DTRC + Sum DTRC'
Sum_DTRC_DKRC_E = zeros(1, length(t)) ;
Sum_DTRC_DKRC_M = zeros(1, length(t)) ;

%% 'Rates of ethane consumption and methane production'
r111_27_plot = zeros(length(t),1) ; % 'ethane consumption'
methane_TOF_plot = zeros(length(t),1) ;

%% 'Compute methane production and ethane consumption rates for every instance'
for jk = 1:length(t)
    current_y_reshape = reshape(y(jk,:), [length(y(end,:)), 1]) ;

    r111_27_plot(jk) = kf_111(27)*x(1)*Ptot*y(jk,17) - kb_111(27)*y(jk,1) ;
    methane_TOF_plot(jk) = -(kf_111(30)*x(4)*Ptot*y(jk,17) - kb_111(30)*y(jk,11)) ;

    % 'Retrieve Q_Gts, Q_ads and R in their original shape from the segmented column vectors'
    Q_Gts = reshape(current_y_reshape(length(x0)+1:length(x0)+length(x0)*length(G_ts_111)), [length(x0), length(G_ts_111)]) ;
    Q_Gads = reshape(current_y_reshape(length(x0)+1+length(x0)*length(G_ts_111):length(x0)+length(x0)*length(G_ts_111)+length(x0)*length(G_ads_111)), [length(x0), length(G_ads_111)]) ;
    R = reshape(current_y_reshape(length(x0)+1+length(x0)*length(G_ts_111)+length(x0)*length(G_ads_111):end), [length(x0),length(x0)]) ;

    % 'Create deep-learning arrays of state variables-theta and gas mole ratios'
    dl_y = dlarray(current_y_reshape) ;
    dl_x = dlarray(x) ;

	% 'Check that site conservation holds'
    y_site_balance = [1, 1, 2, 3, 2, 3, 3, 3, 3, 4, 1, 1, 2, 3, 3, 1, 1] * current_y_reshape(1:length(x0)) ;

	% 'Conservation of gas moles'
    sum(x(1:5)) % 'Sum of Gas partial fractions'
    
    % 'Instantaneous gas partial pressures'
    % 'Not expected to change since this evaluation is at 0 percent conversion'
    PCH3CH3_f = Ptot * x(1) ;
    PCH2CH2_f = Ptot * x(2) ;
    PCHCH_f = Ptot * x(3) ;
    PCH4_f = Ptot * x(4) ;
    PH2_f = Ptot * x(5) ;

    Ptot_f = PCH3CH3_f + PCH2CH2_f + PCHCH_f + PCH4_f + PH2_f ;

    % 'Compute d(ln(TOF))/dG at constant theta using MATLAB DL arrays and reverse propagation of Automatic Differentiation'
    [lnTOF_Ethane, lnTOF_Methane, dlnTOF_E_dGts_thetaconst, dlnTOF_M_dGts_thetaconst, dlnTOF_E_dGads_thetaconst, dlnTOF_M_dGads_thetaconst, dlnTOF_E_dtheta_Gconst,  dlnTOF_M_dtheta_Gconst] = dlfeval(@compute_DKRC_DTRC_ThetaConstant_DL, current_y_reshape(1:length(x0)), dl_y(1:length(x0)), x, dl_x, T, Ptot_f, dl_G_ts_111, G_ts_111, dl_G_ads_111, G_ads_111, dl_EZPEc_gas_111, EZPEc_gas_111, dl_lnQ_gas_111, lnQ_gas_111) ;

    %% 'compute DTRC and DKRC'
    % 'Based on Ethane Consumption Rate'
    n_dkrc_E = -1 * kb * T * (dlnTOF_E_dGts_thetaconst + Q_Gts.' * dlnTOF_E_dtheta_Gconst) ;
    n_dtrc_E = -1 * kb * T * (dlnTOF_E_dGads_thetaconst + Q_Gads.' * dlnTOF_E_dtheta_Gconst) ;

    % 'Based on Methane Production Rate'
    n_dkrc_M = -1 * kb * T * (dlnTOF_M_dGts_thetaconst + Q_Gts.' * dlnTOF_M_dtheta_Gconst) ;
    n_dtrc_M = -1 * kb * T * (dlnTOF_M_dGads_thetaconst + Q_Gads.' * dlnTOF_M_dtheta_Gconst) ;

	%% 'Collect DKRC and DTRC data'
    DKRC_E_combined(:, jk) = extractdata(n_dkrc_E) ;
    DTRC_E_combined(:, jk) = extractdata(n_dtrc_E) ;
    DKRC_M_combined(:, jk) = extractdata(n_dkrc_M) ;
    DTRC_M_combined(:, jk) = extractdata(n_dtrc_M) ;

    Sum_DTRC_E(jk) = sum(DTRC_E_combined(:, jk)) ;
    Sum_DKRC_E(jk) = sum(DKRC_E_combined(:, jk)) ;
    Sum_DTRC_M(jk) = sum(DTRC_M_combined(:, jk)) ;
    Sum_DKRC_M(jk) = sum(DKRC_M_combined(:, jk)) ;

    Sum_DTRC_DKRC_E(jk) = Sum_DTRC_E(jk) + Sum_DKRC_E(jk) ;
    Sum_DTRC_DKRC_M(jk) = Sum_DTRC_M(jk) + Sum_DKRC_M(jk) ;

end

%% 'Plot log-log figure of rates and pseudo-coverages'
figure
loglog(t,r111_27_plot(:,1),'r',t,methane_TOF_plot(:,1),'g',t, y(:,1),'k',t,y(:,2),'b',t,y(:,3),'b',t,y(:,4),'b',t,y(:,5),'b',t,y(:,6),'b',t,y(:,7),'b',t,y(:,8),'b',t,y(:,9),'b',t,y(:,10),'k',t,y(:,11),'b',t,y(:,12),'b',t,y(:,13),'b',t,y(:,14),'b',t,y(:,15),'b',t,y(:,16),'b',t,y(:,17),'b') ;
title('Solution of balance Equation');
xlabel('time t'); ylabel('solution y');
l = legend('C_{2}H_{6} consumption', 'CH_{4} production', 'y 1','y 2','y 3','y 4','y 5','y 6','y 7','y 8','y 9','y 10','y 11','y 12','y 13','y 14','y 15','y 16','y 17') ;
set(l,'Edgecolor',[1 1 1]);

y_raw = y(end,:);
y = reshape(y_raw, [length(y_raw), 1]) ;

%% 'Final total pressure'
Ptot_f = PCH3CH3_f + PCH2CH2_f + PCHCH_f + PCH4_f + PH2_f ;

%% 'Regular MKM'
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

% 'forward rates'
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

rate_combined(:) = r111 ;

% 'final composiite y'
fileID = fopen('coverage.csv','a');
fmt = '%5d\r\n';
fprintf(fileID,fmt,y);
fclose(fileID);

% 'Write out the final net rate, forward and reverse rate of all elementary step'
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

% 'Final TOF at the end'
Ethane_cons_rate=fopen('TOF','a');
Ethane_cons_TOF = r111(27) ;
Methane_prod_TOF = -r111(30) ;
fprintf(Ethane_cons_rate,'%22s %5d %20s %5d %12s %6.3f %12s %6.3f\r\n','Ethane consumption= ', Ethane_cons_TOF,'Methane production= ', Methane_prod_TOF, 'PCH3CH3= ',PCH3CH3,'PH2= ',PH2);
fclose(Ethane_cons_rate);
sites_bal=fopen('sites', 'a');
fprintf(sites_bal, '%6s %6.3f\r\n', 'Site Balance', y_site_balance);
fclose(sites_bal);

%%  'Write out Results'
writematrix(coverages, 'combined_coverages.xls');
writematrix(rate_combined, 'combined_rate.xls');
writematrix(DKRC_E_combined, 'combined_DKRC_E.csv') ;
writematrix(DTRC_E_combined, 'combined_DTRC_E.csv') ;
writematrix(DKRC_M_combined, 'combined_DKRC_M.csv') ;
writematrix(DTRC_M_combined, 'combined_DTRC_M.csv') ;
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
