% clc
% clear all
% close all

addpath('ResultsMC')
addpath('functions')

curr_dir = pwd;
cd('ResultsMC')
cd(dir_case)
load([dir_case])
cd(curr_dir) 

%%

K = 3;

n_sim = length(rho_vec)*length(eps_skew_vec)*length(T_vec)*length(p_vec)

%% Figures

generate_figures = 'yes';

if strcmp(generate_figures, 'yes')
    
for i_rho = 1:length(rho_vec)-1
    rho = rho_vec(i_rho);
        for i_eps_skew = 1: length(eps_skew_vec)
            eps_skew = eps_skew_vec(i_eps_skew);

fH = figure;
h = bar(categorical(p_vec),squeeze(rej_freq_allsim(i_rho,:,:,i_eps_skew))', 'FaceColor','flat');                     

if rho == 0 
    ylim([0 0.5])
else
    ylim([0 1])
end

step        = NaN(length(p_vec),3,5);
step(:,:,1) = kron(ones(length(p_vec),1),[1 0 0]);
step(:,:,2) = kron(ones(length(p_vec),1),[0 1 0]);
step(:,:,3) = kron(ones(length(p_vec),1),[0 0 1]);
step(:,:,4) = kron(ones(length(p_vec),1),[1 0 1]);
h(1).CData = step(:,:,1);
h(2).CData = step(:,:,2);
h(3).CData = step(:,:,3);
h(4).CData = step(:,:,4);

ll = yline(significance_level,'k', 'Linewidth', 2); hold on

xlabel('$p$', 'interpreter','latex')
legend({['T = ', num2str(T_vec(1))], ['T = ', num2str(T_vec(2))], ['T = ', num2str(T_vec(3))], ['T = ', num2str(T_vec(4))]},'Orientation','horizontal','Location','NorthEast','NumColumns',2)

box off
grid

applyhatch(fH,'\-x./',[1 0 0;...
                       0 1 0; ...
                       0 0 1;...
                       1 0 1])

fig_name = ['rho' num2str(rho_vec(i_rho)) '_eps_skew_' num2str(eps_skew_vec(i_eps_skew))]

curr_dir = pwd;
cd('ResultsMC')
cd(dir_case)
print([fig_name,'.pdf'],'-dpdf')
cd(curr_dir) 

close all

    end
end
       
end
