function [Y,z,eps] = zFC_GenData_04(K,T,B,A1,p_true,...
                                psi,...
                                p,...
                                eps_mean,eps_std,eps_skew,eps_kurt, shock_NG)       
    if strcmp(shock_NG,'all')
        eps     = [ ]; 
        for i = 1:K
            eps_step    = pearsrnd(eps_mean,eps_std,eps_skew,eps_kurt,1,T+p);
            eps         = [eps; eps_step]; 
        end
    
    elseif strcmp(shock_NG,'one')
        eps         = randn(K,T+p);
        eps(1,:)    = pearsrnd(eps_mean,eps_std,eps_skew,eps_kurt,1,T+p);
    elseif strcmp(shock_NG,'two')
        eps         = randn(K,T+p);
        eps(2,:)    = pearsrnd(eps_mean,eps_std,eps_skew,eps_kurt,1,T+p);
    elseif strcmp(shock_NG,'three')
        eps         = randn(K,T+p);
        eps(3,:)    = pearsrnd(eps_mean,eps_std,eps_skew,eps_kurt,1,T+p);
    end
    U       = B*eps;

    z = psi(1)*eps(1,:) + psi(2)*eps(2,:) + randn(1,T+p);    
    
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