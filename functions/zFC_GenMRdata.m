function [Y,z,U] = zFC_GenMRdata(Pi_true,p,rho_true,sigma_nu_sq_true,T,eps_skew, B_true)

K       = size(Pi_true,1);
p_true  = size(Pi_true,2)/K;


eps_kurt = max(3,eps_skew^2+2); 

eps     = [ ]; 
for i = 1:K
    eps1    = pearsrnd(0,1,eps_skew,eps_kurt,1,T+p);
    eps     = [eps; eps1]; 
end
U       = B_true*eps;

Y = zeros(K,T+p);
for tt = p+1:T+p
    step = [];
    for pp = 1:p_true
        step = [step;Y(:,tt-pp)];
    end
    Y(:,tt) = Pi_true*step + B_true*eps(:,tt);
end

z = eps(1,:) + rho_true*eps(2,:) + sqrt(sigma_nu_sq_true)*randn(1,T+p);
corr([z;eps(1:2,:)]')    ;
skewness(z) ;
end