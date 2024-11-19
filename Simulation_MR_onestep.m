clc
clear all
close all
addpath(genpath('./functions'));


seed = 123;
addpath('functions')

dir_case = 'ResultsMC\Results_MR_onestep';

%% Set Parameters



load('MR_proxy_param_02.mat')

Pi_true = A_MR(:,const_MR+trend_MR+trend2_MR+1:end);
B_true  = B_MR;

K       = size(Pi_true,1);
p_true  = 1;

companion   = [Pi_true; kron(eye(p_true-1,p_true-1),eye(K)),zeros(K*(p_true-1),K)];
[~, C]      = eig(companion);
eig_all     = flip(sort(diag(abs(C))));
maxroot     = eig_all(1)

%% Run Simulation

T_vec           = [300; 600; 1200; 5000];
p_vec           = [1; 4; 8];
nrep            = 500;

const       = 1;
trend       = 0;
trend2      = 0;

significance_level = 0.1;

test_type = 'one_step';
%test_type = 'two_step'

n_sim = length(rho_vec)*length(T_vec)*length(p_vec)*length(eps_skew_vec)
i_sim = 0;

eta_allsim          = NaN(length(rho_vec),length(T_vec),length(p_vec),length(eps_skew_vec));
rej_freq_allsim     = NaN(length(rho_vec),length(T_vec),length(p_vec),length(eps_skew_vec));

for i_rho = 1:length(rho_vec)
    rho_true = rho_vec(i_rho);
    for i_t = 1:length(T_vec)
        T = T_vec(i_t);
        for i_p = 1:length(p_vec)
            p = p_vec(i_p);
            for i_eps = 1:length(eps_skew_vec)
                eps_skew = eps_skew_vec(i_eps);
                rej_freq_all        = NaN(nrep,1);
                J_statistic_all     = NaN(nrep,1);
                parfor i_rep = 1:nrep
                    %             for i_rep = 1:nrep
                    [Y,z] = zFC_GenMRdata(Pi_true,p,rho_true,sigma_nu_sq_true,T, eps_skew, B_true);
                    
                    [A_hat,U_hat,~] = zFC_VarEstimation_02(Y,p,const, trend,trend2);
                    z_used = z(:,p+1:end);
                    
                    step = U_hat*z_used'/(T-p);
                    step2 = step/step(1);
                    theta_z = step2(2:end);
                    
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
                    
                end
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
                    '_rho_' num2str(rho_true), ...
                    '_s_nu_sq_' num2str(sigma_nu_sq_true),...
                    '_eps_skew_' num2str(eps_skew),...
                    '_p_' num2str(p),...
                    '_T_' num2str(T)];
                    
 
                results_file = fullfile(results_dir, 'Results.xlsx'); 
                % Write simulation name and results
                xlswrite(results_file, cellstr(sim_name), 'Rej_freq', ['A' num2str(i_sim + 1)]);       
                results = squeeze(squeeze(mean(rej_freq_all,1,'omitnan')))  
                xlswrite(results_file, results, 'Rej_freq', ['B' num2str(i_sim + 1)]); 
                cd(curr_dir)
                        
                
                
                rej_freq_allsim(i_rho,i_t,i_p, i_eps)    = squeeze(squeeze(mean(rej_freq_all,1,'omitnan')));
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

save(filename ,  'rej_freq_allsim',...
    'significance_level',...
    'B_true', 'Pi_true','nrep','T_vec','rho_vec','sigma_nu_sq_true',...
    'p_vec','eps_skew_vec','const',...
    'test_type')

cd(curr_dir)

