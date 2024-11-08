format long
clear
close all
delete DTRC_methane.csv;
delete DTRC_ethane.csv;
delete DKRC_methane.csv;
delete DKRC_ethane.csv;
delete SUM_DTRC_methane.csv;
delete SUM_DTRC_ethane.csv;
delete SUM_DKRC_methane.csv;
delete SUM_DKRC_ethane.csv;
delete SUM_DTRC_DKRC_methane.csv ;
delete SUM_DTRC_DKRC_ethane.csv ;

%% 'Constants'
kb = 8.617333262e-5 ; % 'Boltzmanns constant'
h = 4.135667696e-15 ; % 'Plancks constant'

% 'Operating conditions'
T = 573.15 ;

%% 'file names'
species_file = 'species_ads_111.csv' ;
ts_file = 'ts_energy_111.xlsx' ;

%% 'Read Adsorption energies'
G_ads_111_raw = readtable(species_file,'ReadRowNames',true) ;
G_ts_111_raw = readtable(ts_file,'Sheet','ts_energy_111','ReadRowNames',true) ;
G_ads_act_111_raw = readtable(ts_file,'Sheet','ads_act','ReadRowNames',true) ;
G_ads_rxn_111_raw = readtable(ts_file,'Sheet','ads_rxn','ReadRowNames',true) ;
EZPEc_gas_111_raw = readtable(ts_file,'Sheet','EZPEc_gas','ReadRowNames',true) ;
lnQ_gas_111_raw = readtable(ts_file,'Sheet','lnQ_gas','ReadRowNames',true) ;

T_modified = T * 100 ;
T_readable = 'T' + string(T_modified) + '_K' ;
G_ads_111 = G_ads_111_raw.(T_readable) ;
G_ts_111 = G_ts_111_raw.(T_readable) ;
G_ads_act_111 = G_ads_act_111_raw.(T_readable) ;
G_ads_rxn_111 = G_ads_rxn_111_raw.(T_readable) ;
lnQ_gas_111 = lnQ_gas_111_raw.(T_readable) ;

EZPEc_gas_111 = EZPEc_gas_111_raw.('EZPEc_gas') ;

% 'time steps for finite differences evaluation specified as a list'
time_read = readmatrix('time_steps.csv') ;

%% 'Add activation energy of adsorption-desorption steps computed from collision theory to G_ts_111'
G_ts_111(27:31) = G_ads_act_111 ;

DTRC_methane = zeros(length(G_ads_111), length(time_read)) ;
DTRC_ethane = zeros(length(G_ads_111), length(time_read)) ;
DKRC_methane = zeros(length(G_ts_111), length(time_read)) ;
DKRC_ethane = zeros(length(G_ts_111), length(time_read)) ;

energy_perturb = 1e-4 ; % 'set delta G here'

%% 'DTRC loop'
for il = 1:length(G_ads_111)
    rTOF_ethane = zeros(3, length(time_read)) ;
    rTOF_methane = zeros(3, length(time_read)) ;
    ij = il ;
    G_ads_111_iter = [(G_ads_111(ij)-energy_perturb), G_ads_111(ij), (G_ads_111(ij)+energy_perturb)] ;
    
    % 'Get the rates at different energies of G_ads_111_iter'
    for ik = 1:3

        G_ads_111(ij) = G_ads_111_iter(ik) ;
        
        % 'Initialize reaction energies and activation barriers'
        dG_act_111 = zeros(length(G_ts_111), 1) ;
        dG_rxn_111 = zeros(length(G_ts_111), 1) ;
        
        % 'Initialize elementary rates, equilibrium constants and forward constants'
        r111 = zeros(length(G_ts_111),1) ;
        Keq_111 = zeros(length(G_ts_111), 1) ;
        kf_111 = zeros(length(G_ts_111), 1) ;
        kb_111 = zeros(length(G_ts_111), 1) ;
        
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
        
        % '111 Reaction energies of TS steps'
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
        
        %% 'Get Reaction constants'
        for idx = 1:26
            if (dG_act_111(idx) < 0)
                dG_act_111(idx) = 0.001 ;
            end
            Keq_111(idx) = exp(-(dG_rxn_111(idx)) / (kb * T));
            kf_111(idx) = (kb * T / h) * exp(-(dG_act_111(idx)) / (kb * T)) ;
            kb_111(idx) = kf_111(idx) / Keq_111(idx) ;
        end
        
        % 'Activation and reaction energies of gas phase adsorption desorption steps'
        % 'for H2 (index 16), remember that slab and product energy is squared'
        prod_idx = [1, 5, 8, 11, 16] ;
        for i = 1:length(prod_idx)
            if i == length(prod_idx)
                % 'Dissociative adsorption of H2 to 2H(*)'
                dG_act_111(26 + i) = G_ts_111(26 + i) + G_ads_111(17) - 2 * G_ads_111(17) ; % 'No need to subtract gas energy as it is embedded in the G_ts from collision theory'
        
                dG_rxn_111(26 + i) = 2 * G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - 2 * G_ads_111(17) ;
            else
                % 'Other gas phase'
                dG_act_111(26 + i) = G_ts_111(26 + i) - G_ads_111(17) ;
        
                dG_rxn_111(26 + i) = G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - G_ads_111(17) ;
            end
        
            Keq_111(26 + i) = exp(-(dG_rxn_111(26 + i)) / (kb * T)) ;
            kf_111(26 + i) = (kb * T / h) * exp(-(dG_act_111(26 + i)) / (kb * T)) ;
            kb_111(26 + i) = kf_111(26 + i) / Keq_111(26 + i) ;

        end
        
        H2_pressure = 0.1317225 ; % '0.13 atm to bar'
        Ethane_pressure = 0.03343725 ; % '0.033 atm to bar'
        
        PCH3CH3 = Ethane_pressure;
        PCH2CH2 = 0;
        PCHCH = 0;
        PCH4 = 0;
        PH2 = H2_pressure;

        Ptot = PCH3CH3 + PCH2CH2 + PCHCH + PCH4 + PH2 ;
        
	    % 'Initial conditions'
	    y0=zeros(length(G_ads_111),1);
	    y0(17) = 1 ; % '100 percent free sites at the start'

        x(1) = PCH3CH3 / Ptot ;
        x(2) = PCH2CH2 / Ptot ;
        x(3) = PCHCH / Ptot ;
        x(4) = PCH4 / Ptot ;
        x(5) = PH2 / Ptot ;

        % 'Mass matrix'
        M=eye(length(y0)) ;
        M(length(y0),length(y0)) = 0 ; % 'implicit ODE form for the free site balance'
        
        %% 'Solver options'
        optode = odeset('NonNegative',1:17,'Mass',M,'Abstol',1E-15,'RelTol',1E-15) ;
        optlsq = optimset('TolFun',1E-13,'TolX',1E-13) ;        
        
        [t,y]=ode15s(@(t,y)C2_Pt_111_ODE_function(t, y, x, T, Ptot, G_ts_111, G_ads_111, EZPEc_gas_111, lnQ_gas_111),time_read,y0,optode);

        for jk = 1:length(time_read)
            rTOF_ethane(ik, jk)= kf_111(27)*x(1)*Ptot*y(jk,17) - kb_111(27)*y(jk,1) ;
            rTOF_methane(ik, jk) = -(kf_111(30)*x(4)*Ptot*y(jk,17) - kb_111(30)*y(jk,11)) ;
        end
        
        y=y(end,:);
        
        y_site_balance = y(1:17).* [1, 1, 2, 3, 2, 3, 3, 3, 3, 4, 1, 1, 2, 3, 3, 1, 1] ;
        sum(y_site_balance) % 'Sum of site balance'
        sum(x(1:5)) % 'Sum of Gas partial fractions'

    end

    % 'Reassign the species adsorption energy to the original value'
    G_ads_111(ij) = G_ads_111_iter(2);
    G_ads_111(ij) ;

    Gdelta = [(G_ads_111(ij)-energy_perturb) ; G_ads_111(ij) ; (G_ads_111(ij)+energy_perturb)] ;
    X_TRC = [ones(length(Gdelta),1) Gdelta] ;

    for jm = 1:length(time_read)
        log_TOF_methane = log(rTOF_methane(:,jm)) ;
        log_TOF_ethane = log(rTOF_ethane(:,jm)) ;
        TRC_methane = X_TRC\log_TOF_methane ;
        TRC_ethane = X_TRC\log_TOF_ethane ;
        TRC_corr_methane = -kb * T * TRC_methane ;
        TRC_corr_ethane = -kb * T * TRC_ethane ;

        if (real(TRC_corr_methane(2)) < 0)
	        DTRC_methane(il,jm) = -1 * abs(TRC_corr_methane(2)) ;
        else
            DTRC_methane(il,jm) = abs(TRC_corr_methane(2)) ;
        end

        if (real(TRC_corr_ethane(2)) < 0)
	        DTRC_ethane(il,jm) = -1 * abs(TRC_corr_ethane(2)) ;
        else
            DTRC_ethane(il,jm) = abs(TRC_corr_ethane(2)) ;
        end
    end

end
writematrix(DTRC_methane, 'DTRC_methane.csv');
writematrix(DTRC_ethane, 'DTRC_ethane.csv');

%% 'DKRC loop'
for il = 1:length(G_ts_111)
    rTOF_ethane = zeros(3, length(time_read)) ;
    rTOF_methane = zeros(3, length(time_read)) ;
    ij = il ;

    G_ts_111_iter = [G_ts_111(ij)-energy_perturb, G_ts_111(ij), G_ts_111(ij)+energy_perturb] ;
    
    for ik = 1:3

        G_ts_111(ij) = G_ts_111_iter(ik) ;
        
        % 'Initialize reaction energies and activation barriers'
        dG_act_111 = zeros(length(G_ts_111), 1) ;
        dG_rxn_111 = zeros(length(G_ts_111), 1) ;
        
        % 'Initialize elementary rates, equilibrium constants and forward constants'
        r111 = zeros(length(G_ts_111), 1) ;
        Keq_111 = zeros(length(G_ts_111), 1) ;
        kf_111 = zeros(length(G_ts_111), 1) ;
        kb_111 = zeros(length(G_ts_111), 1) ;
        
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
        
        % '111 Reaction energies of TS steps'
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
        
        %% 'Get Reaction constants'
        for idx = 1:26
            if (dG_act_111(idx) < 0)
                dG_act_111(idx) = 0.001 ;
            end
            Keq_111(idx) = exp(-(dG_rxn_111(idx)) / (kb * T));
            kf_111(idx) = (kb * T / h) * exp(-(dG_act_111(idx)) / (kb * T)) ;
            kb_111(idx) = kf_111(idx) / Keq_111(idx) ;
        end
        
        % 'Activation and reaction energies of gas phase adsorption desorption steps'
        % 'for H2 (index 16), remember that slab and product energy is squared'
        prod_idx = [1, 5, 8, 11, 16] ;
        for i = 1:length(prod_idx)
            if i == length(prod_idx)
                % 'Dissociative adsorption of H2 to 2H(*)'
                dG_act_111(26 + i) = G_ts_111(26 + i) + G_ads_111(17) - 2 * G_ads_111(17) ; % 'No need to subtract gas energy as it is embedded in the G_ts from collision theory'
        
                dG_rxn_111(26 + i) = 2 * G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - 2 * G_ads_111(17) ;
            else
                % 'Other gas phase'
                dG_act_111(26 + i) = G_ts_111(26 + i) - G_ads_111(17) ;
        
                dG_rxn_111(26 + i) = G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - G_ads_111(17) ;
            end
        
            Keq_111(26 + i) = exp(-(dG_rxn_111(26 + i)) / (kb * T)) ;
            kf_111(26 + i) = (kb * T / h) * exp(-(dG_act_111(26 + i)) / (kb * T)) ;
            kb_111(26 + i) = kf_111(26 + i) / Keq_111(26 + i) ;

        end
        
        H2_pressure = 0.1317225 ; % '0.13 atm to bar'
        Ethane_pressure = 0.03343725 ; % '0.033 atm to bar'
        
        
        PCH3CH3 = Ethane_pressure;
        PCH2CH2 = 0;
        PCHCH = 0;
        PCH4 = 0;
        PH2 = H2_pressure;

        Ptot = PCH3CH3 + PCH2CH2 + PCHCH + PCH4 + PH2 ;
        
	    % 'Initial conditions'
	    y0=zeros(length(G_ads_111),1);
	    y0(17) = 1 ; % '100 percent free sites at the start'

        x(1) = PCH3CH3 / Ptot ;
        x(2) = PCH2CH2 / Ptot ;
        x(3) = PCHCH / Ptot ;
        x(4) = PCH4 / Ptot ;
        x(5) = PH2 / Ptot ;

        % 'Mass matrix'
        M=eye(length(y0)) ;
        M(length(y0),length(y0)) = 0 ; % 'implicit ODE form for the free site balance'
        
        %% 'Solver options'
        optode = odeset('NonNegative',1:17,'Mass',M,'Abstol',1E-15,'RelTol',1E-15) ;
        optlsq = optimset('TolFun',1E-13,'TolX',1E-13) ;        
        
        [t,y]=ode15s(@(t,y)C2_Pt_111_ODE_function(t, y, x, T, Ptot, G_ts_111, G_ads_111, EZPEc_gas_111, lnQ_gas_111),time_read,y0,optode);

        for jk = 1:length(time_read)
            rTOF_ethane(ik, jk)= kf_111(27)*x(1)*Ptot*y(jk,17) - kb_111(27)*y(jk,1) ;
            rTOF_methane(ik, jk) = -(kf_111(30)*x(4)*Ptot*y(jk,17) - kb_111(30)*y(jk,11)) ;
        end
        
        y=y(end,:);
        
        y_site_balance = y(1:17).* [1, 1, 2, 3, 2, 3, 3, 3, 3, 4, 1, 1, 2, 3, 3, 1, 1] ;
        sum(y_site_balance) % 'Sum of site balance'
        sum(x(1:5)) % 'Sum of Gas partial fractions'

    end

    % 'Reassign the steps transition state energy to the original value'
    G_ts_111(ij) = G_ts_111_iter(2) ;    
    G_ts_111(ij) ;

    Gdelta = [(G_ts_111(ij)-energy_perturb) ; G_ts_111(ij) ; (G_ts_111(ij)+energy_perturb)] ;
    X_KRC = [ones(length(Gdelta),1) Gdelta] ;

    for jm = 1:length(time_read)
        log_TOF_methane = log(rTOF_methane(:,jm)) ;
        log_TOF_ethane = log(rTOF_ethane(:,jm)) ;
        KRC_methane = X_KRC\log_TOF_methane ;
        KRC_ethane = X_KRC\log_TOF_ethane ;
        KRC_corr_methane = -kb * T * KRC_methane ;
        KRC_corr_ethane = -kb * T * KRC_ethane ;

        if (real(KRC_corr_methane(2)) < 0)
	        DKRC_methane(il,jm) = -1 * abs(KRC_corr_methane(2)) ;
        else
            DKRC_methane(il,jm) = abs(KRC_corr_methane(2)) ;
        end

        if (real(KRC_corr_ethane(2)) < 0)
	        DKRC_ethane(il,jm) = -1 * abs(KRC_corr_ethane(2)) ;
        else
            DKRC_ethane(il,jm) = abs(KRC_corr_ethane(2)) ;
        end
    end
end

SUM_DKRC_methane = zeros(length(time_read),1) ;
SUM_DTRC_methane = zeros(length(time_read),1) ;
SUM_DKRC_ethane = zeros(length(time_read),1) ;
SUM_DTRC_ethane = zeros(length(time_read),1) ;

SUM_DTRC_DKRC_methane = zeros(length(time_read),1) ;
SUM_DTRC_DKRC_ethane = zeros(length(time_read),1) ;

for jm = 1:length(time_read)
    SUM_DKRC_methane(jm) = sum(DKRC_methane(:,jm)) ;
    SUM_DKRC_ethane(jm) = sum(DKRC_ethane(:,jm)) ;
    SUM_DTRC_methane(jm) = sum(DTRC_methane(:,jm)) ;
    SUM_DTRC_ethane(jm) = sum(DTRC_ethane(:,jm)) ;

    SUM_DTRC_DKRC_methane(jm) = SUM_DTRC_methane(jm) + SUM_DKRC_methane(jm) ;
    SUM_DTRC_DKRC_ethane(jm) = SUM_DTRC_ethane(jm) + SUM_DKRC_ethane(jm) ;

end

writematrix(DKRC_methane, 'DKRC_methane.csv');
writematrix(DKRC_ethane, 'DKRC_ethane.csv');

writematrix(SUM_DKRC_methane, 'SUM_DKRC_methane.csv');
writematrix(SUM_DKRC_ethane, 'SUM_DKRC_ethane.csv');
writematrix(SUM_DTRC_methane, 'SUM_DTRC_methane.csv');
writematrix(SUM_DTRC_ethane, 'SUM_DTRC_ethane.csv');

writematrix(SUM_DTRC_DKRC_methane, 'SUM_DTRC_DKRC_methane.csv');
writematrix(SUM_DTRC_DKRC_ethane, 'SUM_DTRC_DKRC_ethane.csv');
writematrix(time_read, 'time_evolution.xlsx') ;
