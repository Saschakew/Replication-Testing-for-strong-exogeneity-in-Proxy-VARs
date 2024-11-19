function [J_statistic,p_value,theta_est] = EstAndJTest_03_Donald(U_hat,z_used,thetastart)
 
    [n,T] = size(U_hat);
    
    options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);
    thisloss_noW = @(theta) mean(f_04_Donald(theta,U_hat, z_used),2)' * mean(f_04_Donald(theta,U_hat, z_used),2);
    [theta_prelim,~,~,~,~,~]= fminunc( thisloss_noW ,thetastart,options); 
     
%     step = [U_hat(2:end,:)*z_used'./(U_hat(1,:)*z_used') - theta_prelim;...
%             U_hat(2:end,:)*z_used.^2'./(U_hat(1,:)*z_used.^2') - theta_prelim];
%     W = pinv(step*step'/size(U_hat,2));

%     step = [(U_hat(2:end,:)-theta_prelim*U_hat(1,:)).*z_used;...
%             (U_hat(2:end,:)-theta_prelim*U_hat(1,:)).*z_used.^2];
%     W   = inv(step*step'/size(U_hat,2));

%     W = pinv(f_03(theta_prelim,U_hat,z_used)*f_03(theta_prelim,U_hat,z_used)');

% % %     step = zeros(2*(K-1));
% % %     step2 = f_04(theta_prelim,U_hat, z_used);
% % %     for tt = 1:size(U_hat,2)
% % %         step = step + step2(:,tt)*step2(:,tt)';
% % %     end
% % %     W = pinv(step/size(U_hat,2));   
    
    W = pinv(f_04_Donald(theta_prelim,U_hat, z_used)*f_04_Donald(theta_prelim,U_hat, z_used)'/size(U_hat,2));

    thisloss = @(theta) mean(f_04_Donald(theta,U_hat, z_used),2)' * W * mean(f_04_Donald(theta,U_hat, z_used),2);
    [theta_est,minloss,~,~,~,~]= fminunc( thisloss ,thetastart,options); 

    % Calculate J statistic
    J_statistic = T * minloss;   
    
    % Define the degrees of freedom 
    k = size(f_04_Donald(theta_est,U_hat, z_used),1); % number of moments
    l = size(U_hat,1)-1;% number of unknown ( 
    df = k - l;

    % df = 3*(n-1)^2;

    % Calculate the p-value
    p_value = 1 - chi2cdf(J_statistic, df);
end

