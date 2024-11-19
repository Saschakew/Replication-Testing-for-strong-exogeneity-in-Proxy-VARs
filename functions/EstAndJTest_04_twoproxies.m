function [J_statistic,p_value,theta_est] = EstAndJTest_04_twoproxies(U_hat,z1_used,z2_used,thetastart )
 
    [n,T] = size(U_hat);
    
    options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 1000);
    %% first step
    %thisloss_noW = @(theta) mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2)' *W* mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2);
    thisloss_noW = @(theta) mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2)' *  mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2);
    [theta_est,minloss,~,~,~,~]= fminunc( thisloss_noW ,thetastart,options); 
    
    %% two steps
     W = pinv(f_04_twoproxies(theta_est,U_hat, z1_used,z2_used)*f_04_twoproxies(theta_est,U_hat, z1_used,z2_used)'/size(U_hat,2));
     thisloss = @(theta) mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2)' *W* mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2);
     [theta_est,minloss,~,~,~,~]= fminunc( thisloss ,thetastart,options); 
    %% three steps
    if false
      W = pinv(f_04_twoproxies(theta_est,U_hat, z1_used,z2_used)*f_04_twoproxies(theta_est,U_hat, z1_used,z2_used)'/size(U_hat,2));
     thisloss = @(theta) mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2)' *W* mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2);
     [theta_est,minloss,~,~,~,~]= fminunc( thisloss ,thetastart,options);    
    end
    %% four steps
    if false
      W = pinv(f_04_twoproxies(theta_est,U_hat, z1_used,z2_used)*f_04_twoproxies(theta_est,U_hat, z1_used,z2_used)'/size(U_hat,2));
     thisloss = @(theta) mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2)' *W* mean(f_04_twoproxies(theta,U_hat, z1_used,z2_used),2);
     [theta_est,minloss,~,~,~,~]= fminunc( thisloss ,thetastart,options);    
    end
%%
    if false
        bootnumb=300;
        Best = getBrest(theta_est,restrictions);
        e = getInovations(U_hat',Best); 
        minloss_boot_save = NaN(1,bootnumb);
        for i_boot = 1:bootnumb
        
            eboot = NaN(size(e));
            
            indices = randi(T, 1, T);
            eboot(:,1:2) = e(indices,1:2);
            z1_boot = z1_used(indices);
            z2_boot = z2_used(indices);
            for i_shock = 3:n
                indices = randi(T, 1, T);
                eboot(:,i_shock) = e(indices,i_shock);
            end
            U_boot = Best * eboot';
    
            thisloss_Boot = @(theta) mean(fe(theta,U_boot, z1_boot,z2_boot,restrictions),2)' * mean(fe(theta,U_boot, z1_boot,z2_boot,restrictions),2);
            [~,minloss_boot,~,~,~,~]= fminunc( thisloss_Boot ,thetastart,options); 
    
            minloss_boot_save(i_boot)=minloss_boot;
    
        end
        p_value=sum(minloss<minloss_boot_save)/bootnumb;
    end
 
    % Calculate J statistic
    J_statistic = T * minloss;   
    
    % Define the degrees of freedom 
    k = size(f_04_twoproxies(theta_est,U_hat, z1_used,z2_used),1); % number of moments
    l = size(U_hat,1)-1;% number of unknown ( 
    df = k - l;

    % Calculate the p-value
    p_value = 1 - chi2cdf(J_statistic, df);
end

