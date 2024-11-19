clc
clear all
close all
addpath(genpath('./functions'));


seed = 123;
addpath('functions')

dir_case = 'ResultsMC\Results_nonlin2';

options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);

%%

significance_level = 0.1;

%test_type = 'one_step'

test_type = 'two_step'

K = 3;
N = 1;
K_1 = 1;
K_2 = K-K_1;

% p_true      = 0; assert(strcmp(test_type, 'two_step') )
p_true      = 0;

B = [1 0 1;...
    2 1 4;...
    4 6 6];

A1 = [0.79 0.00 0.25;...
    0.19 0.95 -0.46;...
    0.12 0.00 0.62];


eps_mean = 0;
eps_std = 1;

nrep            = 500;
T_vec           = [150; 300; 600; 1200;5000];
% T_vec           = [5000];
% rho_vec         = -[0;0.1;0.15;0.2;0.3];
corr_z_eps1_vec  = [1];
corr_z_eps2_vec  = -[1; 2; 3 ];
p_vec           = [0 ];
% p_vec           = [0];
eps_skew_vec    = [0.01 0.05 0.1 ]; % kaenzig proxy has E(z^3) = -1.38
% eps_skew_vec    = [2]; % kaenzig proxy has E(z^3) = -1.38

const   = 0;

n_sim = length(corr_z_eps1_vec)*length(corr_z_eps2_vec)*length(eps_skew_vec)*length(T_vec)*length(p_vec)
i_sim = 0;

eta_allsim          = NaN(length(corr_z_eps1_vec),length(corr_z_eps2_vec),length(eps_skew_vec),length(T_vec),length(p_vec));
rej_freq_allsim     = NaN(length(corr_z_eps1_vec),length(corr_z_eps2_vec),length(eps_skew_vec),length(T_vec),length(p_vec));

for i_z1 = 1:length(corr_z_eps1_vec)
    corr_z_eps1 = corr_z_eps1_vec(i_z1);
    for i_z2 = 1:length(corr_z_eps2_vec)
        
        if i_z2 == 1
            psi = [1 0 ];
        elseif i_z2 == 2
            psi = [1 0.1  ];
        elseif i_z2 == 3
            psi = [1 0.25  ];
        end
        
        for i_eps_skew = 1: length(eps_skew_vec)
            eps_skew = eps_skew_vec(i_eps_skew);
            eps_kurt = 3;
            
            
            psi = [psi eps_skew  ];
            
            for i_t = 1:length(T_vec)
                T = T_vec(i_t);
                for i_p = 1:length(p_vec)
                    p = p_vec(i_p);
                    rej_freq_all        = NaN(nrep,1);
                    J_statistic_all     = NaN(nrep,1);
                    parfor i_rep = 1:nrep
                        rng(seed + i_rep); 
                        % for i_rep = 1:nrep
                        % Generate data of length T+p (NOT T+p_true!)
                        [Y,z,eps] = zFC_GenData_06nonlin2(K,T,B,A1,p_true,...
                            psi,...
                            p,...
                            eps_mean,eps_std,eps_skew,eps_kurt);
                        
                        z_sq            = z.^2;
                        
                        z_used = z(:,p+1:end);
                        z_sq_used = z_sq(:,p+1:end);
                        
                        trend=0;
                        trend2=0;
                        
                        [A_hat,U_hat,~] = zFC_VarEstimation_02(Y,p,const, trend,trend2);
                        
                        step = U_hat*z_used'/(T-p);
                        step2 = step/step(1);
                        theta_z = step2(2:end);
                        
                        step = U_hat*z_sq_used'/(T-p);
                        step2 = step/step(1);
                        theta_z_sq = step2(2:end);
                        
                        if strcmp(test_type, 'one_step')
                            all_param_start = [A_hat(:); theta_z];
                            assert(length(all_param_start) == K*(K*p+const+trend+trend2) + (K-1));
                            [J_statistic,p_value,all_param_est] = EstAndJTest_04(Y,z_used,p,const,trend,trend2,all_param_start);
                        elseif strcmp(test_type, 'two_step')
                            thetastart = theta_z;
                            [J_statistic,p_value,theta_est] = EstAndJTest_03(U_hat,z_used,thetastart);
                        end
                        reject1 =  p_value < significance_level;
                        rej_freq_all(i_rep) = reject1;
                        
                        %                         thetastart = [theta_z' mean(U_hat(1,:).*z_used) mean(U_hat(1,:).*z_sq_used) ];
                        %                         [J_statistic,p_value,theta_est] = EstAndJTest(U_hat,z_used,thetastart);
                        %                         reject1 =  p_value < significance_level;
                        %                         rej_freq_all(i_rep) = reject1;
                        
                        J_statistic_all(i_rep) = J_statistic;
                        
                    end % end loop over nrep
                    
                    i_sim = i_sim + 1;
                    n_sim - i_sim
                    
                    
                    
                    
                    
                    %write to Excel
                    curr_dir = pwd;
                    results_dir = fullfile(curr_dir, dir_case);
                    if ~exist(results_dir, 'dir')
                        mkdir(results_dir);
                    end
                    cd(results_dir)
                    
                    sim_name = ['nrep_' num2str(nrep), ...
                        '_iz1_' num2str(i_z1), ...
                        '_psi2_' num2str(psi(2)), ...
                        '_psi3_' num2str(eps_skew),...
                        '_p_' num2str(p),...
                        '_T_' num2str(T)];
                    
 
                    results_file = fullfile(results_dir, 'Results.xlsx'); 
                    % Write simulation name and results
                    xlswrite(results_file, cellstr(sim_name), 'Rej_freq', ['A' num2str(i_sim + 1)]);       
                    results = squeeze(squeeze(mean(rej_freq_all,1,'omitnan')))  
                    xlswrite(results_file, results, 'Rej_freq', ['B' num2str(i_sim + 1)]); 
                    cd(curr_dir)
                            
                    
                    rej_freq_allsim(i_z1,i_z2,i_eps_skew,i_t,i_p)    = squeeze(squeeze(mean(rej_freq_all,1)));
                    
                end
            end
        end
    end
end

% save results to Matlab
curr_dir = pwd;
mkdir(dir_case)
cd(dir_case)

% Extract the base name from dir_case
[~, base_name, ~] = fileparts(dir_case);

% Create a descriptive filename
filename = [base_name,   '.mat'];

save(filename, 'rej_freq_allsim',...
    'significance_level',...
    'B', 'A1','eps_mean','eps_std','nrep','T_vec','corr_z_eps1_vec','corr_z_eps2_vec',...
    'p_vec','eps_skew_vec','const',...
    'test_type','p_true')

cd(curr_dir)

