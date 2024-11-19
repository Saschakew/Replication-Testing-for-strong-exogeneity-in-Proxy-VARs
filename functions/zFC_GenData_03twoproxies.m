function [Y,z1,z2,eps] = zFC_GenData_03twoproxies(K,T,B,A1,p_true,...
                                psi,...
                                p,...
                                eps_mean,eps_std,eps_skew,eps_kurt)       
    eps     = [ ]; 
    for i = 1:K
        eps1    = pearsrnd(eps_mean,eps_std,eps_skew,eps_kurt,1,T+p);
        eps     = [eps; eps1]; 
    end
    U       = B*eps;

    z1 = psi(1)*eps(1,:)   + psi(1)/2 *eps(2,:) + psi(2)*eps(3,:) + randn(1,T+p);   
    z2 = psi(1)/2*eps(1,:) + psi(1) *eps(2,:)   + psi(2)*eps(3,:) + randn(1,T+p);    
    
    if p_true>0
        Y               = NaN(K,T+p);
        Y(:,1:p_true)   = zeros(K,p_true);
        
        for tt = p_true+1:(T+p)
           step     = Y(:,(tt-p_true):(tt-1));
           Y(:,tt)   = A1*step(:) + U(:,tt);
    
        end
    else
        Y=U;
    end

end