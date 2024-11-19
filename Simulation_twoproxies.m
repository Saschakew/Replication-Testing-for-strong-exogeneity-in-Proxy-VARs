clc
clear all
close all
addpath(genpath('./functions'));


seed = 123;
addpath('functions')

dir_case = 'ResultsMC\Results_twoproxies';

options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-15,'StepTolerance', 1e-10,'MaxIterations', 1000);

%%

significance_level = 0.1;

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

A0 = inv(B);
A0sign = diag(sign(diag(A0)));
A0 = A0sign * A0   ;
A0scale = inv(diag(diag(A0)));
A0 = A0scale * A0 ;

A1 = [0.79 0.00 0.25;...
    0.19 0.95 -0.46;...
    0.12 0.00 0.62];


eps_mean = 0;
eps_std = 1;

nrep            = 500;
T_vec           = [150; 300; 600; 1200;5000];
corr_z_eps1_vec  = [0.4, 0.5, 0.6];
corr_z_eps2_vec  = -[0.6, 0.5,0]; % should be renamed to corr_z_eps3_vec here
p_vec           = [0];
eps_skew_vec    = [0;1;2];

const   = 1;

n_sim = length(corr_z_eps1_vec)*length(corr_z_eps2_vec)*length(eps_skew_vec)*length(T_vec)*length(p_vec)
i_sim = 0;

eta_allsim                  = NaN(length(corr_z_eps1_vec),length(corr_z_eps2_vec),length(eps_skew_vec),length(T_vec),length(p_vec));
rej_freq_allsim     = NaN(length(corr_z_eps1_vec),length(corr_z_eps2_vec),length(eps_skew_vec),length(T_vec),length(p_vec));

for i_z1 = 1:length(corr_z_eps1_vec)
    corr_z_eps1 = corr_z_eps1_vec(i_z1);
    for i_z2 = 1:length(corr_z_eps2_vec)
        corr_z_eps2 = corr_z_eps2_vec(i_z2);
        corr_target = [corr_z_eps1;corr_z_eps2];
        func_to_min = @(psi) zFC_get_param_from_corr_twoproxies(psi,corr_target)'*zFC_get_param_from_corr_twoproxies(psi,corr_target);
        [psi,~,~,~,~,~]= fminunc(func_to_min,[0,0],options);
        % % % ww = zFC_get_param_from_corr_twoproxies(psi,corr_target)*zFC_get_param_from_corr_twoproxies(psi,corr_target)';
        % % % func_to_min = @(psi) zFC_get_param_from_corr_twoproxies(psi,corr_target)'*inv(ww)*zFC_get_param_from_corr_twoproxies(psi,corr_target);
        % % % [psi,~,~,~,~,~]= fminunc(func_to_min,[0,0],options); ghjghj
        % % % [psi,~,~,~,~,~]= fminsearch(func_to_min,[0,0]); ghjghj
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
                        [Y,z1,z2,eps] = zFC_GenData_03twoproxies(K,T,B,A1,p_true,...
                            psi,...
                            p,...
                            eps_mean,eps_std,eps_skew,eps_kurt);
                        
                        % mean(z2 .*eps(1,:))
                        
                        % mean(z2 .* Y(3,:))
                        % -A0(3,1) * mean(z2 .* Y(1,:)) - A0(3,2) * mean(z2 .* Y(2,:))
                        %  mean(z2 .* Y(3,:)) +A0(3,1) * mean(z2 .* Y(1,:)) + A0(3,2) * mean(z2 .* Y(2,:))
                        
                        % A0signscale = A0sign * A0scale;
                        % A0(3,3)  * Y(3,:)
                        % -A0(3,1) *   Y(1,:)  - A0(3,2) *  Y(2,:) +   eps(3,:)
                        % -A0(3,1) *   Y(1,:)  - A0(3,2) *  Y(2,:) +  A0signscale(3,3) * eps(3,:)
                        % B(3,1) *   eps(1,:)  + B(3,2) *  eps(2,:) + B(3,3) * eps(3,:)
                        
                        z1_sq            = z1.^2;
                        z2_sq            = z2.^2;
                        
                        z1_used = z1(:,p+1:end);
                        z1_sq_used = z1_sq(:,p+1:end);
                        
                        z2_used = z2(:,p+1:end);
                        z2_sq_used = z2_sq(:,p+1:end);
                        
                        trend=0;
                        trend2=0;
                        
                        U_hat=Y;
                        
                        
                        theta_z = [A0(3:end,1)';A0(3:end,2)'];
                        options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);
                        thisloss_noW = @(theta) mean(f_04_twoproxies_nosynt(theta,U_hat, z1_used,z2_used),2)' * mean(f_04_twoproxies_nosynt(theta,U_hat, z1_used,z2_used),2);
                        [theta_est,minloss,~,~,~,~]= fminunc( thisloss_noW ,theta_z,options);
                        theta_z=theta_est;
                        %W = pinv(f_04_twoproxies(theta_est,U_hat, z1_used,z2_used)*f_04_twoproxies(theta_est,U_hat, z1_used,z2_used)'/size(U_hat,2));
                        
                        
                        if strcmp(test_type, 'one_step')
                            %  all_param_start = [A_hat(:); theta_z];
                            %  assert(length(all_param_start) == K*(K*p+const+trend+trend2) + (K-1));
                            % [J_statistic,p_value,all_param_est] = EstAndJTest_04(Y,z_used,p,const,trend,trend2,all_param_start);
                        elseif strcmp(test_type, 'two_step')
                            
                            thetastart = theta_z;
                            % [J_statistic,p_value,theta_est] = EstAndJTest_03_twoproxies(U_hat,z1_used,z2_used,thetastart,restrictions);
                            [J_statistic,p_value,theta_est] = EstAndJTest_04_twoproxies(U_hat,z1_used,z2_used,thetastart );
                        end
                        reject1 =  p_value < significance_level;
                        rej_freq_all(i_rep) = reject1;
                        
                        
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

save(filename, 'rej_freq_allsim',...
    'significance_level',...
    'B', 'A1','eps_mean','eps_std','nrep','T_vec','corr_z_eps1_vec','corr_z_eps2_vec',...
    'p_vec','eps_skew_vec','const',...
    'test_type','p_true')

cd(curr_dir)






