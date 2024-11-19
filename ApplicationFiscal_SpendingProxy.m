%% Setup
clc
clear all
close all
addpath(genpath('./functions'));
addpath(genpath('./data'));

%% Load and Prepare Data
load('data.mat')

date = MertensRavnData(:,1);
y = MertensRavnData(:,2:4);
z = MertensRavnData(:,7);
 
y = y(:,[2 1 3]);  % Reorder columns

%% Analyze Instrument Properties
fprintf('Skewness of z: %.4f\n', skewness(z))
fprintf('Kurtosis of z: %.4f\n', kurtosis(z))

%% VAR Estimation
lags = 4;
[T, n] = size(y);

% Prepare data for VAR
Y0 = y(1:lags,:);  % Initial conditions
shortY = y(lags+1:end,:);

dummy = zeros(T,1);
dummy(102) = 1;

tmpY = [Y0(end-lags+1:end,:); shortY];
X = zeros(T-lags, n*lags);
for i = 1:lags
    X(:, (lags-i+1)*n-n+1:(lags-i+1)*n) = tmpY(lags-i+1:end-i,:);
end

Xexo = [ones(T-lags,1), dummy(lags+1:end), (1:T-lags)', (1:T-lags).^2'];
X = [X, Xexo];

% Estimate VAR
Atilde = (X \ shortY)';
u = shortY - X * Atilde';
z = z(lags+1:end,:);

% Prepare data for further analysis
z = z';
U_hat = u';

%% Relevance Tests
significance_level = 0.1;

% Linear relevance test
lhs = U_hat(1,:)';
rhs = z';
linear_model = fitlm(rhs, lhs);
disp('Linear Relevance Test:')
disp(linear_model)

% Quadratic relevance test
lhs = U_hat(1,:)';
rhs = z.^2';
quadratic_model = fitlm(rhs, lhs);
disp('Quadratic Relevance Test:')
disp(quadratic_model)

% Detailed relevance tests for each variable
disp('Detailed Relevance Tests:')
for ik = 1:n
    inam = sprintf('r%d', ik);
    olsEst.(inam) = olsest(z', U_hat(ik,:)', true, true);
    
    fprintf('Variable %d:\n', ik)
    fprintf('  F-stat: %4.3f, p-value: %4.3f\n', olsEst.(inam).F, olsEst.(inam).Fpval)
    fprintf('  F-stat (robust): %4.3f, p-value: %4.3f\n', olsEst.(inam).Frobust, olsEst.(inam).Frobustpval)
    fprintf('  R^2: %4.3f, R^2 (adj): %4.3f\n', olsEst.(inam).R2, olsEst.(inam).R2adj)
end

% Detailed relevance tests (z^2) for each variable
disp('Detailed Relevance Tests (z^2):')
for ik = 1:n
    inam = sprintf('r%d', ik);
    olsEst.(inam) = olsest(z.^2', U_hat(ik,:)', true, true);
    
    fprintf('Variable %d:\n', ik)
    fprintf('  F-stat: %4.3f, p-value: %4.3f\n', olsEst.(inam).F, olsEst.(inam).Fpval)
    fprintf('  F-stat (robust): %4.3f, p-value: %4.3f\n', olsEst.(inam).Frobust, olsEst.(inam).Frobustpval)
    fprintf('  R^2: %4.3f, R^2 (adj): %4.3f\n', olsEst.(inam).R2, olsEst.(inam).R2adj)
end

 %% Initial Parameter Estimation
T = length(z);

% Estimate theta for z
step = U_hat * z' / T;
step2 = step / step(1);
theta_z = step2(2:end);

% Estimate theta for z^2
step = U_hat * z.^2' / T;
step2 = step / step(1);
theta_z2 = step2(2:end);

fprintf('Initial theta estimates:\n')
fprintf('theta_z:  '); fprintf('%8.4f ', theta_z); fprintf('\n')
fprintf('theta_z2: '); fprintf('%8.4f ', theta_z2); fprintf('\n\n')

%% Optimization Setup
options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10, ...
                       'StepTolerance', 1e-10, 'MaxIterations', 1000);

% Initial optimization without optimal weighting
W_initial = pinv(diag(diag(f_04(theta_z, U_hat, z) * f_04(theta_z, U_hat, z)' / size(U_hat, 2))));
loss_func_initial = @(theta) mean(f_04(theta, U_hat, z), 2)' * W_initial * mean(f_04(theta, U_hat, z), 2);
[theta_est_initial, ~] = fminunc(loss_func_initial, theta_z, options);

% Optimization with optimal weighting
W_optimal = pinv(f_04(theta_est_initial, U_hat, z) * f_04(theta_est_initial, U_hat, z)' / size(U_hat, 2));
loss_func_optimal = @(theta) mean(f_04(theta, U_hat, z), 2)' * W_optimal * mean(f_04(theta, U_hat, z), 2);
[theta_est_optimal, ~] = fminunc(loss_func_optimal, theta_est_initial, options);

thetastart = theta_est_optimal;

fprintf('Optimized theta estimate:\n')
fprintf('theta_est: '); fprintf('%8.4f ', thetastart); fprintf('\n\n')

%% Estimation and J-Test
[J_statistic, p_value, theta_est] = EstAndJTest_03(U_hat, z, thetastart);

fprintf('Standard GMM Estimation Results:\n')
fprintf('Estimated theta: '); fprintf('%8.4f ', theta_est); fprintf('\n')
fprintf('J-statistic: %8.4f\n', J_statistic)
fprintf('p-value: %8.4f\n', p_value)
fprintf('Reject H0 at %.1f%% level: %s\n\n', significance_level*100, ...
        logical2str(p_value < significance_level))

%% CUE Estimation (Non-conservative)
conservative = 'no';
[J_statistic_cue, p_value_cue, theta_est_cue] = EstAndJTest_03CUE(U_hat, z, thetastart, conservative);

fprintf('CUE Estimation Results (Non-conservative):\n')
fprintf('Estimated theta: '); fprintf('%8.4f ', theta_est_cue); fprintf('\n')
fprintf('J-statistic: %8.4f\n', J_statistic_cue)
fprintf('p-value: %8.4f\n', p_value_cue)
fprintf('Reject H0 at %.1f%% level: %s\n\n', significance_level*100, ...
        logical2str(p_value_cue < significance_level))

%% CUE Estimation (Conservative)
conservative = 'yes';
[J_statistic_cue_cons, p_value_cue_cons, theta_est_cue_cons] = EstAndJTest_03CUE(U_hat, z, thetastart, conservative);

fprintf('CUE Estimation Results (Conservative):\n')
fprintf('Estimated theta: '); fprintf('%8.4f ', theta_est_cue_cons); fprintf('\n')
fprintf('J-statistic: %8.4f\n', J_statistic_cue_cons)
fprintf('p-value: %8.4f\n', p_value_cue_cons)
fprintf('Reject H0 at %.1f%% level: %s\n\n', significance_level*100, ...
        logical2str(p_value_cue_cons < significance_level))

%% Estimation with Kurtosis
[J_statistic_kurt, p_value_kurt, theta_est_kurt] = EstAndJTest_03_wKurt(U_hat, z, thetastart);

fprintf('Estimation Results with Kurtosis:\n')
fprintf('Estimated theta: '); fprintf('%8.4f ', theta_est_kurt); fprintf('\n')
fprintf('J-statistic: %8.4f\n', J_statistic_kurt)
fprintf('p-value: %8.4f\n', p_value_kurt)
fprintf('Reject H0 at %.1f%% level: %s\n\n', significance_level*100, ...
        logical2str(p_value_kurt < significance_level))

%% CUE Estimation with Kurtosis (Non-conservative)
conservative = 'no';
[J_statistic_cue_kurt, p_value_cue_kurt, theta_est_cue_kurt] = EstAndJTest_03CUE_wKurt(U_hat, z, thetastart, conservative);

fprintf('CUE Estimation Results with Kurtosis (Non-conservative):\n')
fprintf('Estimated theta: '); fprintf('%8.4f ', theta_est_cue_kurt); fprintf('\n')
fprintf('J-statistic: %8.4f\n', J_statistic_cue_kurt)
fprintf('p-value: %8.4f\n', p_value_cue_kurt)
fprintf('Reject H0 at %.1f%% level: %s\n\n', significance_level*100, ...
        logical2str(p_value_cue_kurt < significance_level))

%% CUE Estimation with Kurtosis (Conservative)
conservative = 'yes';
[J_statistic_cue_kurt_cons, p_value_cue_kurt_cons, theta_est_cue_kurt_cons] = EstAndJTest_03CUE_wKurt(U_hat, z, thetastart, conservative);

fprintf('CUE Estimation Results with Kurtosis (Conservative):\n')
fprintf('Estimated theta: '); fprintf('%8.4f ', theta_est_cue_kurt_cons); fprintf('\n')
fprintf('J-statistic: %8.4f\n', J_statistic_cue_kurt_cons)
fprintf('p-value: %8.4f\n', p_value_cue_kurt_cons)
fprintf('Reject H0 at %.1f%% level: %s\n\n', significance_level*100, ...
        logical2str(p_value_cue_kurt_cons < significance_level))

%% Helper function
function str = logical2str(logic_val)
    if logic_val
        str = 'Yes';
    else
        str = 'No';
    end
end
