clc
clear all
close all

%% DGP1
dir_case = 'Results_Baseline';
Figures_dgp1

dir_case = 'Results_p1';
Figures_dgp1

dir_case = 'Results_p12';
Figures_dgp1

dir_case = 'Results_FirstSkewed';
Figures_dgp1

dir_case = 'Results_SecondSkewed';
Figures_dgp1

dir_case = 'Results_ThirdSkewed';
Figures_dgp1

psi2 = [0; 0.25; 0.5];
psi3 = [0.8; 0.7; 0.5];
dir_case = 'Results_nonlin';
Figures_nonlin

psi2 = [0; 0.1; 0.25];
psi3 = [0.01; 0.05; 0.1];
dir_case = 'Results_nonlin2';
Figures_nonlin

psi2 = [0; 0.25; 0.5];
psi3 = [1; 0.1; 0];
dir_case = 'Results_nonlin3';
Figures_nonlin

dir_case = 'Results_twoproxies';
Figures_dgp1

dir_case = 'Results_wKurt_OneSynthProxy';
Figures_kurt

dir_case = 'Results_wKurt_TwoSynthProxy';
Figures_kurt

dir_case = 'Results_wKurt_DonaldTest';
Figures_kurt

dir_case = 'Results_wKurt_DonaldTestK';
Figures_kurt

dir_case = 'Results_irr_proxy';
Figures_dgp1

dir_case = 'Results_irr_proxy_cue';
Figures_dgp1

dir_case = 'Results_cueRob';
Figures_dgp1

%% DGP2

dir_case = 'Results_SimMR';
Figures_MR

dir_case = 'Results_SimMR_onestep';
Figures_MR
