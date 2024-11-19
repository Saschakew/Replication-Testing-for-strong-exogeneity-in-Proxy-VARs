function out = zFC_get_param_from_corr_twoexogproxies(psi,corr_z_eps)
    out     = NaN(2,1);
    out(1)  = corr_z_eps(1) - psi(1)/sqrt(psi(1)^2 + psi(2)^2 + 1);
    out(2)  = corr_z_eps(2) - psi(2)/sqrt(psi(1)^2 + psi(2)^2 + 1);
end