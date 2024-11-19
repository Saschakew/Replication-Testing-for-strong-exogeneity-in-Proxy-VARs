clc
clear all
close all
addpath(genpath('./functions'));


seed = 123;
addpath('functions')

dir_case = 'ResultsMC\Results_SecondSkewed';

options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);

%%

significance_level = 0.1;

% test_type = 'one_step'
test_type = 'two_step'

% shock_NG = 'all';     % baseline
% shock_NG = 'one';
shock_NG = 'two';
% shock_NG = 'three';

K = 3;
N = 1;
K_1 = 1;
K_2 = K-K_1;

p_true      = 0; assert(strcmp(test_type, 'two_step') )
% p_true      = 1;

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
corr_z_eps1_vec  = [0.5;0.7;0.9];
corr_z_eps2_vec  = -[0.5 ];
p_vec           = [0];
eps_skew_vec    = [0; 1;2];

% const   = 1;
const   = 0;

if p_true == 0
    assert(const == 0)
    assert(p_vec == 0)
end

n_sim = length(corr_z_eps1_vec)*length(corr_z_eps2_vec)*length(eps_skew_vec)*length(T_vec)*length(p_vec)
i_sim = 0;

eta_allsim                  = NaN(length(corr_z_eps1_vec),length(corr_z_eps2_vec),length(eps_skew_vec),length(T_vec),length(p_vec));
rej_freq_allsim     = NaN(length(corr_z_eps1_vec),length(corr_z_eps2_vec),length(eps_skew_vec),length(T_vec),length(p_vec));

for i_z1 = 1:length(corr_z_eps1_vec)
    corr_z_eps1 = corr_z_eps1_vec(i_z1);
    for i_z2 = 1:length(corr_z_eps2_vec)
        corr_z_eps2 = corr_z_eps2_vec(i_z2);
        corr_target = [corr_z_eps1;corr_z_eps2];
        func_to_min = @(psi) zFC_get_param_from_corr(psi,corr_target)'*zFC_get_param_from_corr(psi,corr_target);
        [psi,~,~,~,~,~]= fminunc( func_to_min,[0,0],options);
        for i_eps_skew = 1: length(eps_skew_vec)
            eps_skew = eps_skew_vec(i_eps_skew);
            eps_kurt = max(eps_skew_vec)^2+2;
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
                        [Y,z,eps] = zFC_GenData_04(K,T,B,A1,p_true,...
                            psi,...
                            p,...
                            eps_mean,eps_std,eps_skew,eps_kurt, shock_NG);
                        z_sq            = z.^2;
                        
                        z_used = z(:,p+1:end);
                        z_sq_used = z_sq(:,p+1:end);
                        
                        trend=0;
                        trend2=0;
                        
                        [A_hat,U_hat,~] = zFC_VarEstimation_02(Y,p,const, trend,trend2);
                        
                        step = U_hat*z_used'/(T-p);
                        step2 = step/step(1);
                        theta_z = step2(2:end);
                        %   W = pinv(f_04(theta_z,U_hat, z_used)*f_04(theta_z,U_hat, z_used)'/size(U_hat,2));
                        
                        step = U_hat*z_sq_used'/(T-p);
                        step2 = step/step(1);
                        theta_z_sq = step2(2:end);
                        
                        if strcmp(test_type, 'one_step')
                            all_param_start = [A_hat(:); theta_z];
                            assert(length(all_param_start) == K*(K*p+const+trend+trend2) + (K-1));
                            [J_statistic,p_value,all_param_est] = EstAndJTest_04(Y,z_used,p,const,trend,trend2,all_param_start);
                        elseif strcmp(test_type, 'two_step')
                            thetastart = theta_z;
                            [J_statistic,p_value,theta_est] = EstAndJTest_03(U_hat,z_used,thetastart  );
                            %  [J_statistic,p_value,theta_est] = EstAndJTest_03(U_hat,z_used,thetastart );
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
                        '_corr_z_eps1_' num2str(corr_z_eps1), ...
                        '_corr_z_eps2_' num2str(corr_z_eps2), ...
                        '_eps_skew_' num2str(eps_skew),...
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

save(filename , 'rej_freq_allsim',...
    'significance_level',...
    'B', 'A1','eps_mean','eps_std','nrep','T_vec','corr_z_eps1_vec','corr_z_eps2_vec',...
    'p_vec','eps_skew_vec','const',...
    'test_type','p_true','shock_NG')

cd(curr_dir)