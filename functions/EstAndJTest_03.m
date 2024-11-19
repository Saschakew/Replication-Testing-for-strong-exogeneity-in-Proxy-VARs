function [J_statistic,p_value,theta_est] = EstAndJTest_03(U_hat,z_used,thetastart)
 
    [n,T] = size(U_hat);
    
    options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);
 

    W = pinv( diag(diag( f_04(thetastart,U_hat, z_used)*f_04(thetastart,U_hat, z_used)'/size(U_hat,2) ) ));
    thisloss_noW = @(theta) mean(f_04(theta,U_hat, z_used),2)' *W* mean(f_04(theta,U_hat, z_used),2);
    [theta_est,~,~,~,~,~]= fminunc( thisloss_noW ,thetastart,options); 
     
    W = pinv(f_04(theta_est,U_hat, z_used)*f_04(theta_est,U_hat, z_used)'/size(U_hat,2));
    thisloss = @(theta) mean(f_04(theta,U_hat, z_used),2)' * W * mean(f_04(theta,U_hat, z_used),2);
    [theta_est,minloss,~,~,~,~]= fminunc( thisloss ,thetastart,options); 
    
 

 


    % Calculate J statistic
    J_statistic = T * minloss;   
    
    % Define the degrees of freedom 
    k = size(f_04(theta_est,U_hat, z_used),1); % number of moments
    l = size(U_hat,1)-1;% number of unknown ( 
    df = k - l;

    % Calculate the p-value
    p_value = 1 - chi2cdf(J_statistic, df);

 
end

