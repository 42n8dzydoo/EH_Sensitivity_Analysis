function F = C2_Pt_111_ODE_function(t, y, x, T, Ptot, G_ts_111, G_ads_111, EZPEc_gas_111, lnQ_gas_111, LI_111_H)
    format long
    
    % 'Constants'
    kb = 8.617333262e-5 ; % 'Boltzmanns constant'
    h = 4.135667696e-15 ; % 'Plancks constant'
    
    %% 'Add LI to G_ads'
    G_ads_111_LI = zeros(length(G_ads_111), 1) ;
    for idx = 1:length(G_ads_111)
        G_ads_111_LI(idx) = G_ads_111(idx) + LI_111_H(idx) * y(16) ;
    end
    
    dG_act_111 = zeros(length(G_ts_111), 1) ;
    dG_rxn_111 = zeros(length(G_ts_111), 1) ;
    dG_act_111_LI = zeros(length(G_ts_111), 1) ; % 'Activation Energies with Lateral interactions'
    dG_rxn_111_LI = zeros(length(G_ts_111), 1) ; % 'Reaction Energies with Lateral interactions'
    
    r111 = zeros(length(G_ts_111), 1) ;
    Keq_111 = zeros(length(G_ts_111), 1) ;
    kf_111 = zeros(length(G_ts_111), 1) ;
    kb_111 = zeros(length(G_ts_111), 1) ;

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

    % 'Activation and reaction energies of gas phase adsorption desorption steps'
    % 'for H2 (index 16), remember that slab and product energy is squared'
    prod_idx = [1, 5, 8, 11, 16] ;
    for i = 1:length(prod_idx)
        if i == length(prod_idx)
            dG_rxn_111(26 + i) = 2 * G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - 2 * G_ads_111(17) ;
            dG_rxn_111_LI(26 + i) = dG_rxn_111(26 + i) + 2 * ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * y(16)) ;
            dG_act_111(26 + i) = G_ts_111(26 + i) + G_ads_111(17) - 2 * G_ads_111(17) ;
            dG_act_111_LI(26 + i) = dG_act_111(26 + i) + 0.5 * (dG_rxn_111_LI(26 + i) - dG_rxn_111(26 + i)) ;
        else
            dG_rxn_111(26 + i) = G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - G_ads_111(17) ;
            dG_rxn_111_LI(26 + i) = dG_rxn_111(26 + i) + ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * y(16)) ;
            dG_act_111(26 + i) = G_ts_111(26 + i) - G_ads_111(17) ;
            dG_act_111_LI(26 + i) = dG_act_111(26 + i) + 0.5 * (dG_rxn_111_LI(26 + i) - dG_rxn_111(26 + i)) ;
        end
    end

    %% 'Get Reaction constants'
    for idx = 1:length(G_ts_111)
        if (dG_act_111_LI(idx) < 0)
            dG_act_111_LI(idx) = 0.001 ;
        end
    
        if (dG_act_111_LI(idx) < dG_rxn_111_LI(idx))
		    dG_act_111_LI(idx) = dG_rxn_111_LI(idx) ;
        end
    
        Keq_111(idx) = exp(-(dG_rxn_111_LI(idx)) / (kb * T));
        kf_111(idx) = (kb * T / h) * exp(-(dG_act_111_LI(idx)) / (kb * T)) ;
        kb_111(idx) = kf_111(idx) / Keq_111(idx) ;
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

    %% 'Initialize ODE for coverages'
    F = zeros(17,1) ;
    
    %% 'ODE for coverages'
    F(1) = r111(27) - (r111(1) + r111(2)) ;
    F(2) = r111(2) - (r111(3) + r111(4) + r111(5)) ;
    F(3) = r111(4) - (r111(6) + r111(7) + r111(8)) ;
    F(4) = r111(7) - (r111(9) + r111(10)) ;
    F(5) = r111(5) + r111(28) - (r111(11) + r111(12)) ;
    F(6) = r111(8) + r111(12) - (r111(13) + r111(14) + r111(15)) ;
    F(7) = r111(10) + r111(14) - (r111(16) + r111(17)) ;
    F(8) = r111(15) + r111(29) - (r111(18) + r111(19)) ;
    F(9) = r111(17) + r111(19) - (r111(20) + r111(21)) ;
    F(10) = r111(21) - (r111(22)) ;
    F(11) = r111(30) - (r111(23)) ;
    F(12) = r111(3) + r111(6) + r111(9) + r111(23) + 2 * (r111(1)) - (r111(24)) ;
    F(13) = r111(3) + r111(13) + r111(16) + r111(24) + 2 * (r111(11)) - (r111(25)) ;
    F(14) = r111(6) + r111(13) + r111(20) + r111(25) + 2 * (r111(18)) - (r111(26)) ;
    F(15) = r111(9) + r111(16) + r111(20) + r111(26) + 2 * (r111(22)) ;
    F(16) = r111(2) + r111(4) + r111(5) + r111(7) + r111(8) + r111(10) + r111(12) + r111(14) + r111(15) + r111(17) + r111(19) + r111(21) + r111(23) + r111(24) + r111(25) + r111(26) + 2 * (r111(31)) ;
    
    %% 'Free site balances as ODE'
    % 'F(17) =  - (r111(1) + r111(2) + r111(9) + r111(10) + r111(14) + r111(15) + r111(17) + r111(19) + r111(23) + r111(26) + r111(27) + r111(30)) - 2 * (r111(3) + r111(4) + r111(5) + r111(6) + r111(7) + r111(8) + r111(11) + r111(12) + r111(13) + r111(16) + r111(21) + r111(22) + r111(24) + r111(25) + r111(28) + r111(31)) - 3 * (r111(18) + r111(20) + r111(29)) ;'
    
    %% 'Free site balances as Algebraic expression'
    F(17) = 1 - y(1) - y(2) - 2 * y(3) - 3 * y(4) - 2 * y(5) - 3 * y(6) - 3 * y(7) - 3 * y(8) - 3 * y(9) - 4 * y(10) - y(11) - y(12) - 2 * y(13) - 3 * y(14) - 3 * y(15) - y(16) - y(17) ;

end