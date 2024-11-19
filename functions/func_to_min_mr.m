function out = func_to_min_mr(param, Sigma_hat, theta_nonnorm)
    
    B = [theta_nonnorm, reshape(param,3,2)];
    ll = reshape(Sigma_hat,3^2,1) - reshape(B*B', 3^2,1);

end