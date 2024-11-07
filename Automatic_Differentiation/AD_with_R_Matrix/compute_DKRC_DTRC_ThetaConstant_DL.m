function [lnTOF_Ethane, lnTOF_Methane, dlnTOF_E_dGts_thetaconst, dlnTOF_M_dGts_thetaconst, dlnTOF_E_dGads_thetaconst, dlnTOF_M_dGads_thetaconst, dlnTOF_E_dtheta_Gconst,  dlnTOF_M_dtheta_Gconst] = compute_DKRC_DTRC_ThetaConstant_DL(y, dl_y, x, dl_x, T, Ptot, dl_G_ts_111, G_ts_111, dl_G_ads_111, G_ads_111, dl_EZPEc_gas_111, EZPEc_gas_111, dl_lnQ_gas_111, lnQ_gas_111)

kb = 8.617333262e-5; % 'Boltzmanns constant'
h = 4.135667696e-15; % 'Plancks constant'

%% 'Get the Reaction energies and Activation Energies'
% 'Initialize Activation and Reaction Energies'
dG_act_111 = zeros(length(G_ts_111), 1) ;
dG_rxn_111 = zeros(length(G_ts_111), 1) ;

% 'Initialize elementary reaction rates and rate constants'
r111 = zeros(length(G_ts_111), 1) ;
Keq_111 = zeros(length(G_ts_111), 1) ;
kf_111 = zeros(length(G_ts_111), 1) ;
kb_111 = zeros(length(G_ts_111), 1) ;

%% 'Corresponding DL-Arrays of initialized reaction parameters'
dl_dG_act_111 = dlarray(zeros(length(dl_G_ts_111), 1)) ;
dl_dG_rxn_111 = dlarray(zeros(length(dl_G_ts_111), 1)) ;

dl_r111 = dlarray(r111) ;
dl_Keq_111 = dlarray(Keq_111) ;
dl_kf_111 = dlarray(kf_111) ;
dl_kb_111 = dlarray(kb_111) ;

%% '111 Activation energies of TS steps'
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

%% '111 Reaction energies  of TS steps'
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

%% 'dl-arrays for deltaG reaction of TS steps'
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
        % 'Dissociative adsorption of H2 to 2H(*)'
        dG_act_111(26 + i) = G_ts_111(26 + i) + G_ads_111(17) - 2 * G_ads_111(17) ;
        dl_dG_act_111(26 + i) = dl_G_ts_111(26 + i) + dl_G_ads_111(17) - 2 * dl_G_ads_111(17) ;
        Keq_111(26 + i) = (exp(-2 * G_ads_111(prod_idx(i)) / (kb * T))) / (exp((kb * T * lnQ_gas_111(i) - EZPEc_gas_111(i) - 2 * G_ads_111(17)) / (kb * T))) ;
        dl_Keq_111(26 + i) = (exp(-2 * dl_G_ads_111(prod_idx(i)) / (kb * T))) / (exp((kb * T * dl_lnQ_gas_111(i) - dl_EZPEc_gas_111(i) - 2 * dl_G_ads_111(17)) / (kb * T))) ;
    else
        % 'Other gas phase'
        dG_act_111(26 + i) = G_ts_111(26 + i) - G_ads_111(17) ;
        dl_dG_act_111(26 + i) = dl_G_ts_111(26 + i) - dl_G_ads_111(17) ;
        Keq_111(26 + i) = (exp(-1 * G_ads_111(prod_idx(i)) / (kb * T))) / (exp((kb * T * lnQ_gas_111(i) - EZPEc_gas_111(i) - G_ads_111(17)) / (kb * T))) ;
        dl_Keq_111(26 + i) = (exp(-1 * dl_G_ads_111(prod_idx(i)) / (kb * T))) / (exp((kb * T * dl_lnQ_gas_111(i) - dl_EZPEc_gas_111(i) - dl_G_ads_111(17)) / (kb * T))) ;
    end
    dG_rxn_111(26 + i) = -kb * T * log(Keq_111(26 + i)) ;
    dl_dG_rxn_111(26 + i) = -kb * T * log(dl_Keq_111(26 + i)) ;

    kf_111(26 + i) = (kb * T / h) * exp(-(dG_act_111(26 + i)) / (kb * T)) ;
    kb_111(26 + i) = kf_111(26 + i) / Keq_111(26 + i) ;

    dl_kf_111(26 + i) = (kb * T / h) * exp(-(dl_dG_act_111(26 + i)) / (kb * T)) ;
    dl_kb_111(26 + i) = dl_kf_111(26 + i) / dl_Keq_111(26 + i) ;
end

%% 'Overall Rate Expressions'
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

dl_r111(1) = dl_kf_111(1)*dl_y(1)*dl_y(17) - dl_kb_111(1)*dl_y(12)*dl_y(12) ;
dl_r111(2) = dl_kf_111(2)*dl_y(1)*dl_y(17) - dl_kb_111(2)*dl_y(2)*dl_y(16) ;
dl_r111(3) = dl_kf_111(3)*dl_y(2)*dl_y(17)^2 - dl_kb_111(3)*dl_y(12)*dl_y(13) ;
dl_r111(4) = dl_kf_111(4)*dl_y(2)*dl_y(17)^2 - dl_kb_111(4)*dl_y(3)*dl_y(16) ;
dl_r111(5) = dl_kf_111(5)*dl_y(2)*dl_y(17)^2 - dl_kb_111(5)*dl_y(5)*dl_y(16) ;
dl_r111(6) = dl_kf_111(6)*dl_y(3)*dl_y(17)^2 - dl_kb_111(6)*dl_y(12)*dl_y(14) ;
dl_r111(7) = dl_kf_111(7)*dl_y(3)*dl_y(17)^2 - dl_kb_111(7)*dl_y(4)*dl_y(16) ;
dl_r111(8) = dl_kf_111(8)*dl_y(3)*dl_y(17)^2 - dl_kb_111(8)*dl_y(6)*dl_y(16) ;
dl_r111(9) = dl_kf_111(9)*dl_y(4)*dl_y(17) - dl_kb_111(9)*dl_y(12)*dl_y(15) ;
dl_r111(10) = dl_kf_111(10)*dl_y(4)*dl_y(17) - dl_kb_111(10)*dl_y(7)*dl_y(16) ;
dl_r111(11) = dl_kf_111(11)*dl_y(5)*dl_y(17)^2 - dl_kb_111(11)*dl_y(13)*dl_y(13) ;
dl_r111(12) = dl_kf_111(12)*dl_y(5)*dl_y(17)^2 - dl_kb_111(12)*dl_y(6)*dl_y(16) ;
dl_r111(13) = dl_kf_111(13)*dl_y(6)*dl_y(17)^2 - dl_kb_111(13)*dl_y(13)*dl_y(14) ;
dl_r111(14) = dl_kf_111(14)*dl_y(6)*dl_y(17) - dl_kb_111(14)*dl_y(7)*dl_y(16) ;
dl_r111(15) = dl_kf_111(15)*dl_y(6)*dl_y(17) - dl_kb_111(15)*dl_y(8)*dl_y(16) ;
dl_r111(16) = dl_kf_111(16)*dl_y(7)*dl_y(17)^2 - dl_kb_111(16)*dl_y(13)*dl_y(15) ;
dl_r111(17) = dl_kf_111(17)*dl_y(7)*dl_y(17) - dl_kb_111(17)*dl_y(9)*dl_y(16) ;
dl_r111(18) = dl_kf_111(18)*dl_y(8)*dl_y(17)^3 - dl_kb_111(18)*dl_y(14)*dl_y(14) ;
dl_r111(19) = dl_kf_111(19)*dl_y(8)*dl_y(17) - dl_kb_111(19)*dl_y(9)*dl_y(16) ;
dl_r111(20) = dl_kf_111(20)*dl_y(9)*dl_y(17)^3 - dl_kb_111(20)*dl_y(14)*dl_y(15) ;
dl_r111(21) = dl_kf_111(21)*dl_y(9)*dl_y(17)^2 - dl_kb_111(21)*dl_y(10)*dl_y(16) ;
dl_r111(22) = dl_kf_111(22)*dl_y(10)*dl_y(17)^2 - dl_kb_111(22)*dl_y(15)*dl_y(15) ;
dl_r111(23) = dl_kf_111(23)*dl_y(11)*dl_y(17) - dl_kb_111(23)*dl_y(12)*dl_y(16) ;
dl_r111(24) = dl_kf_111(24)*dl_y(12)*dl_y(17)^2 - dl_kb_111(24)*dl_y(13)*dl_y(16) ;
dl_r111(25) = dl_kf_111(25)*dl_y(13)*dl_y(17)^2 - dl_kb_111(25)*dl_y(14)*dl_y(16) ;
dl_r111(26) = dl_kf_111(26)*dl_y(14)*dl_y(17) - dl_kb_111(26)*dl_y(15)*dl_y(16) ;
dl_r111(27) = dl_kf_111(27)*dl_x(1)*Ptot*dl_y(17) - dl_kb_111(27)*dl_y(1) ;
dl_r111(28) = dl_kf_111(28)*dl_x(2)*Ptot*dl_y(17)^2 - dl_kb_111(28)*dl_y(5) ;
dl_r111(29) = dl_kf_111(29)*dl_x(3)*Ptot*dl_y(17)^3 - dl_kb_111(29)*dl_y(8) ;
dl_r111(30) = dl_kf_111(30)*dl_x(4)*Ptot*dl_y(17) - dl_kb_111(30)*dl_y(11) ;
dl_r111(31) = dl_kf_111(31)*dl_x(5)*Ptot*dl_y(17)^2 - dl_kb_111(31)*dl_y(16)^2 ;

Ethane_cons_TOF = abs(dl_r111(27)) ; % 'Absolute to ensure positive value always for log evaluation'
Methane_prod_TOF = -dl_r111(30) ;

% 'Ethane_cons_TOF'
lnTOF_Ethane = log(Ethane_cons_TOF) ;
lnTOF_Methane = log(Methane_prod_TOF) ;

% 'Ethane consumption'
dlnTOF_E_dGts_thetaconst = dlgradient(lnTOF_Ethane, dl_G_ts_111) ;
dlnTOF_E_dGads_thetaconst = dlgradient(lnTOF_Ethane, dl_G_ads_111) ;
dlnTOF_E_dtheta_Gconst = dlgradient(lnTOF_Ethane, dl_y) ;

% 'Methane production'
dlnTOF_M_dGts_thetaconst = dlgradient(lnTOF_Methane, dl_G_ts_111) ;
dlnTOF_M_dGads_thetaconst = dlgradient(lnTOF_Methane, dl_G_ads_111) ;
dlnTOF_M_dtheta_Gconst = dlgradient(lnTOF_Methane, dl_y) ;

end