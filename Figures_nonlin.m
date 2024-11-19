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

n_sim = length(corr_z_eps1_vec)*length(corr_z_eps2_vec)*length(eps_skew_vec)*length(T_vec)*length(p_vec)


%% Figures

generate_figures = 'yes';

if strcmp(generate_figures, 'yes')

nfig = size(rej_freq_allsim,3)


    for i_z2 = 1:size(rej_freq_allsim,3)
                for i_p = 1:length(p_vec)
                    p = p_vec(i_p);

fH = figure;
h = bar(categorical(psi2),squeeze(rej_freq_allsim(1,:,i_z2,:))', 'FaceColor','flat');                     

    ylim([0 1])

step        = NaN(3,3,5);
step(:,:,1) = kron(ones(3,1),[1 0 0]);
step(:,:,2) = kron(ones(3,1),[0 1 0]);
step(:,:,3) = kron(ones(3,1),[0 0 1]);
step(:,:,4) = kron(ones(3,1),[1 0 1]);
step(:,:,5) = kron(ones(3,1),[0 1 1]);
h(1).CData = step(:,:,1);
h(2).CData = step(:,:,2);
h(3).CData = step(:,:,3);
h(4).CData = step(:,:,4);
h(5).CData = step(:,:,5);

ll = yline(significance_level,'k', 'Linewidth', 2); hold on

xlabel('$\psi_2$', 'interpreter','latex')

legend({['T = ', num2str(T_vec(1))], ['T = ', num2str(T_vec(2))], ['T = ', num2str(T_vec(3))], ['T = ', num2str(T_vec(4))], ['T = ', num2str(T_vec(5))]},'Orientation','horizontal','Location','NorthWest','NumColumns',2)

box off
grid

applyhatch(fH,'\-x./',[1 0 0;...
                       0 1 0; ...
                       0 0 1;...
                       1 0 1; ...
                       0 1 1])

fig_name = ['psi3_' num2str(psi3(i_z2))...
           '_p_' num2str(p_vec(i_p))]

curr_dir = pwd;
cd('ResultsMC')
cd(dir_case)
print([fig_name,'.pdf'],'-dpdf')
cd(curr_dir)  

close all
    end
end

end
