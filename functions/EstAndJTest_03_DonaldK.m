function [J_statistic,p_value,theta_est] = EstAndJTest_03_DonaldK(U_hat,z_used,thetastart)
 
    [n,T] = size(U_hat);
    
    options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);
 
    thisloss_noW = @(theta) mean(f_04_Donald2K(theta,U_hat, z_used),2)' * mean(f_04_Donald2K(theta,U_hat, z_used),2);
    [theta_prelim,~,~,~,~,~]= fminunc( thisloss_noW ,thetastart,options); 
   
    W = pinv(f_04_Donald2K(theta_prelim,U_hat, z_used)*f_04_Donald2K(theta_prelim,U_hat, z_used)'/size(U_hat,2));

    thisloss = @(theta) mean(f_04_Donald2K(theta,U_hat, z_used),2)' * W * mean(f_04_Donald2K(theta,U_hat, z_used),2);
    [theta_est,minloss,~,~,~,~]= fminunc( thisloss ,thetastart,options); 
    
  
 %   W = @(theta) pinv(f_04_Donald2K(theta,U_hat, z_used)*f_04_Donald2K(theta,U_hat, z_used)'/size(U_hat,2));
 %   thisloss = @(theta) mean(f_04_Donald2K(theta,U_hat, z_used),2)' * W(theta) * mean(f_04_Donald2K(theta,U_hat, z_used),2);
  %  [theta_est,minloss,~,~,~,~]= fminunc( thisloss ,thetastart,options); 

    % Calculate J statistic
    J_statistic = T * minloss;   
    
    % Define the degrees of freedom 
     k = size(f_04_Donald2K(theta_est,U_hat, z_used),1); % number of moments
     l = size(U_hat,1)-1;% number of unknown ( 
     df = k - l;

     JK_statistic = (J_statistic - (df))/sqrt(2*df);


    % Calculate the p-value with chi
    p_value = 1 - chi2cdf(J_statistic, df);

    % Calculate the p-value with normal
   % p_value_one_side = 1 - normcdf(JK_statistic);
   % p_value_two_sided = 2 * min(p_value_one_side, 1 - p_value_one_side);
   % p_value = p_value_two_sided;
end

