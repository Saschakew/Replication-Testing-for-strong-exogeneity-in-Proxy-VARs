function out = zFC_get_param_from_corr_twoproxies(psi,corr_z_eps)
    out     = NaN(2,1);
    out(1)  = corr_z_eps(1) - psi(1)/sqrt(psi(1)^2 + 0.25*psi(1)^2 + psi(2)^2 + 1);
    out(2)  = corr_z_eps(2) - psi(2)/sqrt(psi(1)^2 + 0.25*psi(1)^2 + psi(2)^2 + 1);

    % out(1)  = corr_z_eps(1) - 0.5*psi(1)/sqrt(0.5^2*psi(1)^2 + 0.25^2*psi(1)^2 + psi(2)^2 + 1);
    % out(2)  = corr_z_eps(2) - psi(2)/sqrt(0.5^2*psi(1)^2 + 0.25^2*psi(1)^2 + psi(2)^2 + 1);
    
    % out = (sum(corr_z_eps) - sum(psi)/sqrt(1.25*psi(1)^2 + psi(2)^2 + 1))^2;

    % out     = NaN(3,1);
    % out(1)  = corr_z_eps(1)     - psi(1)/sqrt(psi(1)^2 + psi(2)^2 + psi(3)^2 + 1);
    % out(2)  = corr_z_eps(1)/2   - psi(2)/sqrt(psi(1)^2 + psi(2)^2 + psi(3)^2 + 1);
    % out(3)  = corr_z_eps(2)     - psi(3)/sqrt(psi(1)^2 + psi(2)^2 + psi(3)^2 + 1);
    
end