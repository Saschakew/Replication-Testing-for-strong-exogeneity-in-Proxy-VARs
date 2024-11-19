function [Y,z,eps] = zFC_GenData_06nonlin3(K,T,B,A1,p_true,...
                                psi,...
                                p,...
                                eps_mean,eps_std,eps_skew,eps_kurt)       
    eps     = [ ]; 
    for i = 1:K
        eps1    = pearsrnd(0,1,0,3,1,T+p);
        eps     = [eps; eps1]; 
    end
    U       = B*eps;
    

    z_ts = ones(T+1,1)*psi(3);
    z = zeros(T,1);
    a = 0.9;
    for t=1:T
        shock = psi(1) * eps(1,t);
        if psi(3) + a * z_ts(t) + shock >0
            z_ts(t+1) = psi(3) + a *  z_ts(t) + shock;
            z(t) = shock ;
        else 
            z(t) = max(shock, -( psi(3) + a * z_ts(t))* rand);
            z_ts(t+1) = a *  z_ts(t) + z(t);
        end   
        z(t) = z(t) + psi(2)*eps(2,t)  + randn;
    end
    z = z'; 
    

 
 

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