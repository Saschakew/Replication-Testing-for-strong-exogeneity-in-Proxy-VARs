function [J_statistic,p_value,all_param_est] = EstAndJTest_04(Y,z_used,p,const,trend,trend2,all_param_start)
 
    [K,T] = size(Y);
    
    options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);

    thisloss_noW = @(all_param) mean(f_05(all_param,Y, z_used,p,const,trend,trend2),2)' * mean(f_05(all_param,Y, z_used,p,const,trend,trend2),2);
    [all_param_prelim,~,~,~,~,~]= fminunc( thisloss_noW,all_param_start,options); 
    
    W = pinv(f_05(all_param_prelim,Y, z_used,p,const,trend,trend2)*f_05(all_param_prelim,Y, z_used,p,const,trend,trend2)'/size(z_used,2));

%     W = pinv(f_05(all_param_start,Y, z_used,p,const,trend,trend2)*f_05(all_param_start,Y, z_used,p,const,trend,trend2)'/size(z_used,2));

    thisloss = @(all_param) mean(f_05(all_param,Y, z_used,p,const,trend,trend2),2)' * W * mean(f_05(all_param,Y, z_used,p,const,trend,trend2),2);
    [all_param_est,minloss,~,~,~,~]= fminunc(thisloss,all_param_start,options); 

    % Calculate J statistic
    J_statistic = T * minloss;   
    
    % Define the degrees of freedom 
    k = size(f_05(all_param_est,Y, z_used,p,const,trend,trend2),1); % number of moments
    l = K*(K*p+const+trend+trend2)+K-1;% number of unknown ( 
    df = k - l;

    % Calculate the p-value
    p_value = 1 - chi2cdf(J_statistic, df);
end

