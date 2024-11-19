clc
clear all
close all


rng(123)

dir_case = 'Results_SimTax_2023_10_10';

curr_dir = pwd;
cd(dir_case)

load(['WS_', dir_case])

cd(curr_dir) 

%%

K = 3;

n_sim = length(rho_vec)*length(eps_skew_vec)*length(T_vec)*length(p_vec)

%% Excel output

generate_excel = 'no';

if strcmp(generate_excel, 'yes')
    
i_sim = 0;

for i_rho = 1:length(rho_vec)
    rho = rho_vec(i_rho);
    for i_s_nu = 1:length(sigma_nu_vec)
    sigma_nu_sq = sigma_nu_vec(i_s_nu);   
        corr_eps_z_analytical       = [1/sqrt(1 + sigma_nu_sq) rho/sqrt(rho^2 + sigma_nu_sq) zeros(1, K-2)];
        for i_eps_skew = 1: length(eps_skew_vec)
            eps_skew = eps_skew_vec(i_eps_skew);
            eps_kurt = max(eps_skew^2+2,3); 
            corr_eps_z_sq_analytical    = [eps_skew/sqrt(eps_kurt + rho^4*eps_kurt) ...
                                           0 zeros(1, K-2)];
            for i_t = 1:length(T_vec)
                T = T_vec(i_t);
                for i_p = 1:length(p_vec)
                    p = p_vec(i_p);
        
            i_sim = i_sim + 1;
            n_sim - i_sim
            
%             write to Excel
            curr_dir = pwd;
            mkdir(dir_case)      
            cd(dir_case)
            sim_name = ['nrep_' num2str(nrep), ...
                        '_nboot_' num2str(nboot), ...
                        '_rho_' num2str(rho), ...
                        '_s_nu_sq_' num2str(sigma_nu_sq),...
                        '_eps_skew_' num2str(eps_skew),...
                        '_p_' num2str(p),...
                        '_T_' num2str(T)];
            
            xlswrite('Results.xls',cellstr(sim_name),'Corr_z_eps_analytical',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',corr_eps_z_analytical,'Corr_z_eps_analytical',['B', num2str(i_sim + 1)])

            xlswrite('Results.xls',cellstr(sim_name),'Corr_z_eps_analytical',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',corr_eps_z_sq_analytical,'Corr_z_sq_eps_analytical',['B', num2str(i_sim + 1)])
            
            xlswrite('Results.xls',cellstr(sim_name),'Corr_z_eps_empirical',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',mean(corr_eps_z_empirical_all,1),'Corr_z_eps_empirical',['B', num2str(i_sim + 1)])

            xlswrite('Results.xls',cellstr(sim_name),'Corr_z_sq_eps_empirical',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',mean(corr_eps_z_sq_empirical_all,1),'Corr_z_sq_eps_empirical',['B', num2str(i_sim + 1)])
            
            xlswrite('Results.xls',cellstr(sim_name),'Rej_freq',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',rej_freq_allsim_Hausman(i_rho,i_s_nu,i_eps_skew,i_t,i_p),'Rej_freq',['B', num2str(i_sim + 1)])

            xlswrite('Results.xls',cellstr(sim_name),'Rej_freq',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',rej_freq_allsim_J1(i_rho,i_s_nu,i_eps_skew,i_t,i_p),'Rej_freq',['C', num2str(i_sim + 1)])

            xlswrite('Results.xls',cellstr(sim_name),'Rej_freq',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',rej_freq_allsim_J2(i_rho,i_s_nu,i_eps_skew,i_t,i_p),'Rej_freq',['D', num2str(i_sim + 1)])

            xlswrite('Results.xls',cellstr(sim_name),'Rej_freq',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',rej_freq_allsim_J3(i_rho,i_s_nu,i_eps_skew,i_t,i_p),'Rej_freq',['E', num2str(i_sim + 1)])

            xlswrite('Results.xls',cellstr(sim_name),'Rej_freq',['A', num2str(i_sim + 1)])
            xlswrite('Results.xls',rej_freq_allsim_J4(i_rho,i_s_nu,i_eps_skew,i_t,i_p),'Rej_freq',['F', num2str(i_sim + 1)])
            
            cd(curr_dir) 
                 
                end 
            end
        end
    end
end

end
%% Figures

generate_figures = 'yes';

if strcmp(generate_figures, 'yes')
    
for i_rho = 1:length(rho_vec)-1
    rho = rho_vec(i_rho);
%     for i_s_nu = 1:length(sigma_nu_vec)
%     sigma_nu_sq = sigma_nu_vec(i_s_nu);   
%         corr_eps_z_analytical       = [1/sqrt(1 + sigma_nu_sq) rho/sqrt(rho^2 + sigma_nu_sq) zeros(1, K-2)];
        for i_eps_skew = 1: length(eps_skew_vec)
            eps_skew = eps_skew_vec(i_eps_skew);
%             for i_t = 1:length(T_vec)
%                 T = T_vec(i_t);
%                 for i_p = 1:length(p_vec)
%                     p = p_vec(i_p);

fH = figure;
h = bar(categorical(p_vec),squeeze(rej_freq_allsim(i_rho,:,:,i_eps_skew))', 'FaceColor','flat');                     
% % %             rej_freq_allsim(i_rho,i_t,i_p, i_eps)    = squeeze(squeeze(mean(rej_freq_all,1,'omitnan')));

if rho == 0 
    ylim([0 0.5])
%     ylim([0 1])
else
    ylim([0 1])
end

step        = NaN(length(p_vec),3,5);
step(:,:,1) = kron(ones(length(p_vec),1),[1 0 0]);
step(:,:,2) = kron(ones(length(p_vec),1),[0 1 0]);
step(:,:,3) = kron(ones(length(p_vec),1),[0 0 1]);
% step(:,:,4) = kron(ones(length(eps_skew_vec),1),[1 0 1]);
% step(:,:,5) = kron(ones(length(eps_skew_vec),1),[0 1 1]);
h(1).CData = step(:,:,1);
h(2).CData = step(:,:,2);
h(3).CData = step(:,:,3);
% h(4).CData = step(:,:,4);
% h(5).CData = step(:,:,5);

ll = yline(significance_level,'k', 'Linewidth', 2); hold on

xlabel('$p$', 'interpreter','latex')

% legend({'Hausman', 'J1', 'J2'},'Orientation','horizontal','Location','NorthEast','NumColumns',2)
% legend({'Hausman', 'J1', 'J2', 'J3', 'J4'},'Orientation','horizontal','Location','NorthEast','NumColumns',2)
legend({['T = ', num2str(T_vec(1))], ['T = ', num2str(T_vec(2))], ['T = ', num2str(T_vec(3))], },'Orientation','horizontal','Location','NorthWest','NumColumns',2)
% legend({['T = ', num2str(T_vec(1))], ['T = ', num2str(T_vec(2))], ['T = ', num2str(T_vec(3))], ['T = ', num2str(T_vec(4))], ['T = ', num2str(T_vec(5))]},'Orientation','horizontal','Location','NorthEast','NumColumns',2)

box off
grid

% applyhatch(fH,'\-x./',[1 0 0;...
%                        0 1 0; ...
%                        0 0 1;...
%                        1 0 1])
% applyhatch(fH,'\-x./',[1 0 0;...
%                        0 1 0; ...
%                        0 0 1;...
%                        1 0 1; ...
%                        0 1 1])
% applyhatch(fH,'\-x.',[1 0 0; 0 1 0; 0 0 1; 1 0 1]), hold on
applyhatch(fH,'\-x',[1 0 0;...
                    0 1 0;...
                    0 0 1])

fig_name = ['rho' num2str(rho_vec(i_rho)) '_eps_skew_' num2str(eps_skew_vec(i_eps_skew))]

curr_dir = pwd;
mkdir(dir_case)      
cd(dir_case)
print([fig_name,'.pdf'],'-dpdf')
cd(curr_dir) 

close all
%                 end
%             end
%         end
    end
end
       
end
