format long
clear
close all

%% 'Constants'
kb = 8.617333262e-5 ; % 'Boltzmanns constant'
h = 4.135667696e-15 ; % 'Plancks constant'

% 'Operating conditions'
T = 573.15 ;

%% 'DFT Energy and Paramter file names'
species_file = 'species_ads_111.csv' ;
ts_file = 'ts_energy_111_new_2.xlsx' ;

%% 'Read Adsorption energies'
[G_ads_111, G_ts_111, G_ads_act_111, G_ads_rxn_111, EZPEc_gas_111, lnQ_gas_111, LI_params_111] = read_Energies(T, species_file, ts_file, 111) ;

%% 'Add adsorption-desorption steps activation energy to G_ts_111'
G_ts_111(27:31) = G_ads_act_111 ;

%% 'Read Lateral Interactions Parameters'
LI_111_H = LI_params_111.('H_LI') ;

LI_111_H = ActivateLI(LI_111_H, "no") ;

H2_pressure = 0.1317225 ; % '0.13 atm'
Ethane_pressure = 0.03343725 ; % '0.033 atm'
coverages = zeros(length(G_ads_111), 1) ;

% 'time steps for finite differences evaluation specified as a list'
time_read = readmatrix('time_steps.csv') ;

%% 'Add activation energy of adsorption-desorption steps computed from collision theory to G_ts_111'
G_ts_111(27:31) = G_ads_act_111 ;

DTRC_methane = zeros(length(G_ads_111), length(time_read)) ;
DTRC_ethane = zeros(length(G_ads_111), length(time_read)) ;
DKRC_methane = zeros(length(G_ts_111), length(time_read)) ;
DKRC_ethane = zeros(length(G_ts_111), length(time_read)) ;
sigma_methane = zeros(length(G_ads_111), length(time_read)) ;
sigma_ethane = zeros(length(G_ads_111), length(time_read)) ;
SUM_DKRC_methane = zeros(1,length(time_read)) ;
SUM_DTRC_methane = zeros(1,length(time_read)) ;
SUM_DKRC_ethane = zeros(1,length(time_read)) ;
SUM_DTRC_ethane = zeros(1,length(time_read)) ;

SUM_DTRC_DKRC_methane = zeros(1,length(time_read)) ;
SUM_DTRC_DKRC_ethane = zeros(1,length(time_read)) ;

%% %%%%%%%%%%%%%%%%%%%% 'Set delta G here' %%%%%%%%%%%%%%%%%%%%%%%%%
energy_perturb = 1e-8 ;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 'Delete pre-existing files to write new ones'
file_prefix = ["DTRC_methane" "DTRC_ethane" "DKRC_methane" "DKRC_ethane" "sigma_methane" "sigma_ethane" ...
    "SUM_DTRC_methane" "SUM_DTRC_ethane" "SUM_DKRC_methane" "SUM_DKRC_ethane" "SUM_DTRC_DKRC_methane" ...
    "SUM_DTRC_DKRC_ethane" "coverages" "true_coverages"] ;

for ki = 1 : length(file_prefix)
    baseFileName = file_prefix(ki) ;
    fullFileName = baseFileName + "_" + string(energy_perturb) + ".csv" ;
    
    delete(fullFileName); 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 'DTRC loop' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for il = 1:length(G_ads_111)
    rTOF_ethane = zeros(3, length(time_read)) ;
    rTOF_methane = zeros(3, length(time_read)) ;
    ij = il ;
    G_ads_111_iter = [(G_ads_111(ij)-energy_perturb), G_ads_111(ij), (G_ads_111(ij)+energy_perturb)] ;
    
    % 'Get the rates at different energies of G_ads_111_iter'
    for ik = 1:3

        G_ads_111(ij) = G_ads_111_iter(ik) ;
        
        PCH3CH3 = Ethane_pressure;
        PCH2CH2 = 0;
        PCHCH = 0;
        PCH4 = 0;
        PH2 = H2_pressure;

        Gases_0 = [PCH3CH3; PCH2CH2; PCHCH; PCH4; PH2] ;

        Ptot = sum(Gases_0) ;
        
	    %% 'Initial conditions'
	    y0=zeros(length(G_ads_111),1);
	    y0(17) = 1 ; % '100 percent free sites at the start'

        x = Gases_0 / Ptot ;

        %% 'Mass matrix'
        M=eye(length(y0)) ;
        M(length(y0), length(y0)) = 0 ; % 'Algebraic form for the free site balance'
        
        %% 'Set solver options'
        optode = odeset('NonNegative',1:17,'Mass',M,'Abstol',1E-15,'RelTol',1E-15) ;
        optlsq = optimset('TolFun',1E-13,'TolX',1E-13) ;        
        
        [t,y]=ode15s(@(t,y)C2_Pt_111_ODE_function(t, y, x, T, Ptot, G_ts_111, G_ads_111, EZPEc_gas_111, lnQ_gas_111, LI_111_H),time_read,y0,optode);

        for jk = 1:length(t)
            % 'Initialize reaction energies and activation barriers'
            dG_act_111 = zeros(length(G_ts_111), 1) ;
            dG_rxn_111 = zeros(length(G_ts_111), 1) ;
            dG_act_111_LI = zeros(length(G_ts_111), 1) ; % 'Activation Energies with Lateral interactions'
            dG_rxn_111_LI = zeros(length(G_ts_111), 1) ; % 'Reaction Energies with Lateral interactions'
            
            % 'Initialize elementary rates, equilibrium constants, forward rate constants, and backward rate constants'
            r111 = zeros(length(G_ts_111),1) ;
            Keq_111 = zeros(length(G_ts_111), 1) ;
            kf_111 = zeros(length(G_ts_111), 1) ;
            kb_111 = zeros(length(G_ts_111), 1) ;
    
            %% 'Add LI to G_ads'
            G_ads_111_LI = zeros(length(G_ads_111), 1) ;
            for idx = 1:length(G_ads_111)
                G_ads_111_LI(idx) = G_ads_111(idx) + LI_111_H(idx) * y(jk,16) ;
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
        
            % 'Activation and reaction energies of gas phase adsorption desorption steps'
            % 'for H2 (index 16), remember that slab and product energy is squared'
            prod_idx = [1, 5, 8, 11, 16] ;
            for i = 1:length(prod_idx)
                if i == length(prod_idx)
                    dG_rxn_111(26 + i) = 2 * G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - 2 * G_ads_111(17) ;
                    dG_rxn_111_LI(26 + i) = dG_rxn_111(26 + i) + 2 * ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * y(jk, 16)) ;
                    dG_act_111(26 + i) = G_ts_111(26 + i) + G_ads_111(17) - 2 * G_ads_111(17) ;
                    dG_act_111_LI(26 + i) = dG_act_111(26 + i) + 0.5 * (dG_rxn_111_LI(26 + i) - dG_rxn_111(26 + i)) ;
                else
                    dG_rxn_111(26 + i) = G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - G_ads_111(17) ;
                    dG_rxn_111_LI(26 + i) = dG_rxn_111(26 + i) + ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * y(jk, 16)) ;
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

            rTOF_ethane(ik, jk) = kf_111(27)*x(1)*Ptot*y(jk,17) - kb_111(27)*y(jk,1) ;
            rTOF_methane(ik, jk) = -(kf_111(30)*x(4)*Ptot*y(jk,17) - kb_111(30)*y(jk,11)) ;

        end
        
        if ik == 2
            time_coverages = y ;
        end
        sum(x(1:5)) % 'Sum of Gas partial fractions'

    end

    % 'Reassign the species adsorption energy to the original unperturbed value'
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
            DTRC_methane(il, jm) = -1 * abs(TRC_corr_methane(2)) ;
        else
            DTRC_methane(il, jm) = abs(TRC_corr_methane(2)) ;
        end
    
        if (real(TRC_corr_ethane(2)) < 0)
            DTRC_ethane(il, jm) = -1 * abs(TRC_corr_ethane(2)) ;
        else
            DTRC_ethane(il, jm) = abs(TRC_corr_ethane(2)) ;
        end
    
        current_y = reshape(time_coverages(end,:), [length(time_coverages(end,:)), 1]) ;
        site_counts = [1, 1, 2, 3, 2, 3, 3, 3, 3, 4, 1, 1, 2, 3, 3, 1, 1] ;
        y_site_balance = site_counts * current_y ;
        true_coverages = site_counts' .* current_y ; % 'The unpertubed case for reference'
        sigma_ethane(il, jm) = -1 * DTRC_ethane(il, jm) / true_coverages(il) ;
        sigma_methane(il, jm) = -1 * DTRC_methane(il, jm) / true_coverages(il) ;

    end
end

%% 'Write DTRC and sigma values'
writematrix(DTRC_methane, 'DTRC_methane_' + string(energy_perturb) + '.csv');
writematrix(DTRC_ethane, 'DTRC_ethane_' + string(energy_perturb) + '.csv');
writematrix(sigma_methane, 'sigma_methane_' + string(energy_perturb) + '.csv');
writematrix(sigma_ethane, 'sigma_ethane_' + string(energy_perturb) + '.csv');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 'DKRC loop' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for il = 1:length(G_ts_111)
    rTOF_ethane = zeros(3, length(time_read)) ;
    rTOF_methane = zeros(3, length(time_read)) ;
    ij = il ;
    G_ts_111_iter = [G_ts_111(ij)-energy_perturb, G_ts_111(ij), G_ts_111(ij)+energy_perturb] ;
    
    % 'Get the rates at different energies of G_ts_111_iter'
    for ik = 1:3

        G_ts_111(ij) = G_ts_111_iter(ik) ;
        
        PCH3CH3 = Ethane_pressure;
        PCH2CH2 = 0;
        PCHCH = 0;
        PCH4 = 0;
        PH2 = H2_pressure;

        Gases_0 = [PCH3CH3; PCH2CH2; PCHCH; PCH4; PH2] ;

        Ptot = sum(Gases_0) ;
        
	    %% 'Initial conditions'
	    y0=zeros(length(G_ads_111),1);
	    y0(17) = 1 ; % '100 percent free sites at the start'

        x = Gases_0 / Ptot ;

        %% 'Mass matrix'
        M=eye(length(y0)) ;
        M(length(y0), length(y0)) = 0 ; % 'Algebraic form for the free site balance'
        
        %% 'Set solver options'
        optode = odeset('NonNegative',1:17,'Mass',M,'Abstol',1E-15,'RelTol',1E-15) ;
        optlsq = optimset('TolFun',1E-13,'TolX',1E-13) ;        
        
        [t,y]=ode15s(@(t,y)C2_Pt_111_ODE_function(t, y, x, T, Ptot, G_ts_111, G_ads_111, EZPEc_gas_111, lnQ_gas_111, LI_111_H),time_read,y0,optode);

        for jk = 1:length(t)
            % 'Initialize reaction energies and activation barriers'
            dG_act_111 = zeros(length(G_ts_111), 1) ;
            dG_rxn_111 = zeros(length(G_ts_111), 1) ;
            dG_act_111_LI = zeros(length(G_ts_111), 1) ; % 'Activation Energy with Lateral interactions'
            dG_rxn_111_LI = zeros(length(G_ts_111), 1) ; % 'Reaction Energies with Lateral interactions'
            
            % 'Initialize elementary rates, equilibrium constants, forward rate constants, and backward rate constants'
            r111 = zeros(length(G_ts_111),1) ;
            Keq_111 = zeros(length(G_ts_111), 1) ;
            kf_111 = zeros(length(G_ts_111), 1) ;
            kb_111 = zeros(length(G_ts_111), 1) ;
    
            %% 'Add LI to G_ads'
            G_ads_111_LI = zeros(length(G_ads_111), 1) ;
            for idx = 1:length(G_ads_111)
                G_ads_111_LI(idx) = G_ads_111(idx) + LI_111_H(idx) * y(jk,16) ;
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
        
            % 'Activation and reaction energies of gas phase adsorption desorption steps'
            % 'for H2 (index 16), remember that slab and product energy is squared'
            prod_idx = [1, 5, 8, 11, 16] ;
            for i = 1:length(prod_idx)
                if i == length(prod_idx)
                    dG_rxn_111(26 + i) = 2 * G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - 2 * G_ads_111(17) ;
                    dG_rxn_111_LI(26 + i) = dG_rxn_111(26 + i) + 2 * ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * y(jk, 16)) ;
                    dG_act_111(26 + i) = G_ts_111(26 + i) + G_ads_111(17) - 2 * G_ads_111(17) ;
                    dG_act_111_LI(26 + i) = dG_act_111(26 + i) + 0.5 * (dG_rxn_111_LI(26 + i) - dG_rxn_111(26 + i)) ;
                else
                    dG_rxn_111(26 + i) = G_ads_111(prod_idx(i)) - (EZPEc_gas_111(i) - kb * T * lnQ_gas_111(i)) - G_ads_111(17) ;
                    dG_rxn_111_LI(26 + i) = dG_rxn_111(26 + i) + ((LI_111_H(prod_idx(i)) - LI_111_H(17)) * y(jk, 16)) ;
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

            rTOF_ethane(ik, jk) = kf_111(27)*x(1)*Ptot*y(jk,17) - kb_111(27)*y(jk,1) ;
            rTOF_methane(ik, jk) = -(kf_111(30)*x(4)*Ptot*y(jk,17) - kb_111(30)*y(jk,11)) ;

        end

        if ik == 2
            time_coverages = y ;
        end
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
            DKRC_methane(il, jm) = -1 * abs(KRC_corr_methane(2)) ;
        else
            DKRC_methane(il, jm) = abs(KRC_corr_methane(2)) ;
        end
    
        if (real(KRC_corr_ethane(2)) < 0)
            DKRC_ethane(il, jm) = -1 * abs(KRC_corr_ethane(2)) ;
        else
            DKRC_ethane(il, jm) = abs(KRC_corr_ethane(2)) ;
        end
    end
end

for xl = 1:length(time_read)
    SUM_DKRC_methane(:, xl) = sum(DKRC_methane(:, xl)) ;
    SUM_DKRC_ethane(:, xl) = sum(DKRC_ethane(:, xl)) ;
    SUM_DTRC_methane(:, xl) = sum(DTRC_methane(:, xl)) ;
    SUM_DTRC_ethane(:, xl) = sum(DTRC_ethane(:, xl)) ;
    
    SUM_DTRC_DKRC_methane(:, xl) = SUM_DTRC_methane(:, xl) + SUM_DKRC_methane(:, xl) ;
    SUM_DTRC_DKRC_ethane(:, xl) = SUM_DTRC_ethane(:, xl) + SUM_DKRC_ethane(:, xl) ;
end

%% 'Write DKRC and values of Sums'
writematrix(DKRC_methane, 'DKRC_methane_' + string(energy_perturb) + '.csv');
writematrix(DKRC_ethane, 'DKRC_ethane_' + string(energy_perturb) + '.csv');
writematrix(time_coverages(end,:)', 'coverages_' + string(energy_perturb) + '.csv')
writematrix(true_coverages, 'true_coverages_' + string(energy_perturb) + '.csv')

writematrix(SUM_DKRC_methane, 'SUM_DKRC_methane_' + string(energy_perturb) + '.csv');
writematrix(SUM_DKRC_ethane, 'SUM_DKRC_ethane_' + string(energy_perturb) + '.csv');
writematrix(SUM_DTRC_methane, 'SUM_DTRC_methane_' + string(energy_perturb) + '.csv');
writematrix(SUM_DTRC_ethane, 'SUM_DTRC_ethane_' + string(energy_perturb) + '.csv');

writematrix(SUM_DTRC_DKRC_methane, 'SUM_DTRC_DKRC_methane_' + string(energy_perturb) + '.csv');
writematrix(SUM_DTRC_DKRC_ethane, 'SUM_DTRC_DKRC_ethane_' + string(energy_perturb) + '.csv');

%% 'Plot values of Sums'
figure; 
% 'First subplot (position 1)'
subplot(2, 1, 1) ; 
semilogx(time_read', SUM_DKRC_methane, time_read', SUM_DTRC_methane, time_read', SUM_DTRC_DKRC_methane) ;
yline([1 0 -1],'--') ;
ylim([-8 8]) ;
xlim([1e-12 1]) ;
title('Based on Methane Production') ;

% 'Second subplot (position 2)'
subplot(2, 1, 2) ;
semilogx(time_read', SUM_DKRC_ethane, time_read', SUM_DTRC_ethane, time_read', SUM_DTRC_DKRC_ethane) ;
yline([1 0 -1],'--') ;
ylim([-8 8]) ;
xlim([1e-12 1]) ;
title('Based on Ethane Consumption') ;

%%%%%%%%%%%%%%%% 'Inline Functions' %%%%%%%%%%%%%%%%%%
%% 'Read DFT Energies and Parameters from Worksheet'
function [G_ads, G_ts, G_ads_act, G_ads_rxn, EZPEc_gas, lnQ_gas, LI_params] = read_Energies(T, species_file, ts_file, facet)

    %% 'Read Adsorption energies'
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

%% 'Activate or deactivate LI parameters'
function LI_arr_new = ActivateLI(LI_arr, turn_on)
    if turn_on == "yes"
        LI_arr_new = LI_arr ;
    elseif turn_on == "no"
        LI_arr_new = zeros(length(LI_arr), 1) ;
    end
end