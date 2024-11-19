function [J_statistic,p_value,theta_est] = EstAndJTest_03CUE_wKurt(U_hat,z_used,thetastart,conservative)
 
    [n,T] = size(U_hat);
    
    options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);
 
 
      W = @(theta) pinv(f_04_wKurt(theta,U_hat, z_used)*f_04_wKurt(theta,U_hat, z_used)'/size(U_hat,2));
     %W = @(theta) pinv((f_04_wKurt(theta,U_hat, z_used)-mean(f_04_wKurt(theta,U_hat, z_used),2))*(f_04_wKurt(theta,U_hat, z_used)-mean(f_04_wKurt(theta,U_hat, z_used),2))'/size(U_hat,2));
    thisloss = @(theta) mean(f_04_wKurt(theta,U_hat, z_used),2)' * W(theta) * mean(f_04_wKurt(theta,U_hat, z_used),2);
    [theta_est,minloss,~,~,~,~]= fminunc( thisloss ,thetastart,options); 
    
 

 


    % Calculate J statistic
    J_statistic = T * minloss;   
    
    % Define the degrees of freedom 
    k = size(f_04_wKurt(theta_est,U_hat, z_used),1); % number of moments
    l = size(U_hat,1)-1;% number of unknown ( 
    if strcmp(conservative,'no')
        df = k - l;
    elseif strcmp(conservative,'yes')
        df = k;
    end

    % Calculate the p-value
    p_value = 1 - chi2cdf(J_statistic, df);
end

