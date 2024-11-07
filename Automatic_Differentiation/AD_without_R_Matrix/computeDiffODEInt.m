function [J_full, dfdGts_full, dfdGads_full] = computeDiffODEInt(y, x, T, Ptot, G_ts_111, G_ads_111, EZPEc_gas_111, lnQ_gas_111)
% 'import casadi functionalities from main script'
import casadi.*

kb = 8.617333262e-5; % 'Boltzmanns constant'
h = 4.135667696e-15; % 'Plancks constant'

%% 'Initialize arrays to hold reaction parameters'
dG_act_111 = zeros(length(G_ts_111),1) ;

%% 'Intialize fake-deeplearning arrays to hold constants'
dl_EZPEc_gas_111 = EZPEc_gas_111(:) ;
dl_lnQ_gas_111 = lnQ_gas_111(:) ;
dl_x = x(:) ;

%% 'Casadi matrices to hold calculated reaction constants'
dl_dG_act_111 = MX(1, length(G_ts_111)) ;
dl_dG_rxn_111 = MX(1, length(G_ts_111)) ;
dl_r111 = MX(1, length(G_ts_111)) ;
dl_Keq_111 = MX(1, length(G_ts_111)) ;
dl_kf_111 = MX(1, length(G_ts_111)) ;
dl_kb_111 = MX(1, length(G_ts_111)) ;

%% 'Create differentiable casadi variable matrices'
dl_F = MX.sym('dl_F', 1, length(y)) ;
dl_y = MX.sym('dl_y', 1, length(y)) ;
dl_G_ts_111 = MX.sym('dl_G_ts_111', 1, length(G_ts_111)) ;
dl_G_ads_111 = MX.sym('dl_G_ads_111', 1, length(G_ads_111)) ;

%% 'Get the Reaction energies and Activation Energies'
% '111 Activation energies of TS steps'
% 'I only used this in the for-loop that ensures activation energy is non-negative'
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

%% 'dl_deltaG reaction of TS steps'
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
    end

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
        dl_dG_act_111(26 + i) = dl_G_ts_111(26 + i) + dl_G_ads_111(17) - 2 * dl_G_ads_111(17) ;

        dl_dG_rxn_111(26 + i) = 2 * dl_G_ads_111(prod_idx(i)) - (dl_EZPEc_gas_111(i) - kb * T * dl_lnQ_gas_111(i)) - 2 * dl_G_ads_111(17) ;
    else
        % 'Other gas phase'
        dl_dG_act_111(26 + i) = dl_G_ts_111(26 + i) - dl_G_ads_111(17) ;

        dl_dG_rxn_111(26 + i) = dl_G_ads_111(prod_idx(i)) - (dl_EZPEc_gas_111(i) - kb * T * dl_lnQ_gas_111(i)) - dl_G_ads_111(17) ;
    end

    dl_Keq_111(26 + i) = exp(-(dl_dG_rxn_111(26 + i)) / (kb * T)) ;
    dl_kf_111(26 + i) = (kb * T / h) * exp(-(dl_dG_act_111(26 + i)) / (kb * T)) ;
    dl_kb_111(26 + i) = dl_kf_111(26 + i) / dl_Keq_111(26 + i) ;

end

%% 'DL-array of reaction rates'
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

%% 'ODE function'
dl_F(1) = dl_r111(27) - (dl_r111(1) + dl_r111(2)) ;
dl_F(2) = dl_r111(2) - (dl_r111(3) + dl_r111(4) + dl_r111(5)) ;
dl_F(3) = dl_r111(4) - (dl_r111(6) + dl_r111(7) + dl_r111(8)) ;
dl_F(4) = dl_r111(7) - (dl_r111(9) + dl_r111(10)) ;
dl_F(5) = dl_r111(5) + dl_r111(28) - (dl_r111(11) + dl_r111(12)) ;
dl_F(6) = dl_r111(8) + dl_r111(12) - (dl_r111(13) + dl_r111(14) + dl_r111(15)) ;
dl_F(7) = dl_r111(10) + dl_r111(14) - (dl_r111(16) + dl_r111(17)) ;
dl_F(8) = dl_r111(15) + dl_r111(29) - (dl_r111(18) + dl_r111(19)) ;
dl_F(9) = dl_r111(17) + dl_r111(19) - (dl_r111(20) + dl_r111(21)) ;
dl_F(10) = dl_r111(21) - (dl_r111(22)) ;
dl_F(11) = dl_r111(30) - (dl_r111(23)) ;
dl_F(12) = dl_r111(3) + dl_r111(6) + dl_r111(9) + dl_r111(23) + 2 * (dl_r111(1)) - (dl_r111(24)) ;
dl_F(13) = dl_r111(3) + dl_r111(13) + dl_r111(16) + dl_r111(24) + 2 * (dl_r111(11)) - (dl_r111(25)) ;
dl_F(14) = dl_r111(6) + dl_r111(13) + dl_r111(20) + dl_r111(25) + 2 * (dl_r111(18)) - (dl_r111(26)) ;
dl_F(15) = dl_r111(9) + dl_r111(16) + dl_r111(20) + dl_r111(26) + 2 * (dl_r111(22)) ;
dl_F(16) = dl_r111(2) + dl_r111(4) + dl_r111(5) + dl_r111(7) + dl_r111(8) + dl_r111(10) + dl_r111(12) + dl_r111(14) + dl_r111(15) + dl_r111(17) + dl_r111(19) + dl_r111(21) + dl_r111(23) + dl_r111(24) + dl_r111(25) + dl_r111(26) + 2 * (dl_r111(31)) ;
dl_F(17) =  - (dl_r111(1) + dl_r111(2) + dl_r111(9) + dl_r111(10) + dl_r111(14) + dl_r111(15) + dl_r111(17) + dl_r111(19) + dl_r111(23) + dl_r111(26) + dl_r111(27) + dl_r111(30)) - 2 * (dl_r111(3) + dl_r111(4) + dl_r111(5) + dl_r111(6) + dl_r111(7) + dl_r111(8) + dl_r111(11) + dl_r111(12) + dl_r111(13) + dl_r111(16) + dl_r111(21) + dl_r111(22) + dl_r111(24) + dl_r111(25) + dl_r111(28) + dl_r111(31)) - 3 * (dl_r111(18) + dl_r111(20) + dl_r111(29)) ;
% 'Note that dl_F(17) is evaluated explicitly, not implicitly as in main script'

%% 'Compute the desired 2-Dimenstional differentials symbolically'
% 'Jacobian, dfdy'
J_raw = jacobian(dl_F, dl_y) ;
% 'dfdGts'
dfdGts_raw = jacobian(dl_F, dl_G_ts_111) ;
% 'dfdGads'
dfdGads_raw = jacobian(dl_F, dl_G_ads_111) ;

%% 'Obtain their symbolic functions'
J_fn = Function('J_fn', {dl_y,dl_G_ts_111,dl_G_ads_111}, {J_raw}) ;
dfdGts_fn = Function('dfdGts_fn', {dl_y,dl_G_ts_111,dl_G_ads_111}, {dfdGts_raw}) ;
dfdGads_fn = Function('dfdGads_fn', {dl_y,dl_G_ts_111,dl_G_ads_111}, {dfdGads_raw}) ;

%% 'Evaluate the functions at current theta, Gts and Gads'
J_full = full(J_fn(y, G_ts_111, G_ads_111)) ;
dfdGts_full = full(dfdGts_fn(y, G_ts_111, G_ads_111)) ;
dfdGads_full = full(dfdGads_fn(y, G_ts_111, G_ads_111)) ;

end