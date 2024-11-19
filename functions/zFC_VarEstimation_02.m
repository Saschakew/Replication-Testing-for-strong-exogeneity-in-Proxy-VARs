function [A_hat_var,U_hat_var,Sigma_hat_var] = zFC_VarEstimation_02(Y,p,const,trend,trend2)
    T       = size(Y,2);
    k       = size(Y,1);
    Y_used  = Y(:,p+1:end);
    T_used  = size(Y_used,2);

    if const == 1
        Z = ones(1,T_used);
    elseif const == 0
        Z = [];        
    end

    if trend == 1
        Z = [Z; 1:T_used];     
    end
    if trend2 == 1
        Z = [Z; (1:T_used).^2];     
    end
    
    for pp = 1:p
       Z = [Z; Y(:,p+1-pp:end-pp)]; 
    end
    
    if isempty(Z)
        A_hat_var = [];
        U_hat_var = Y_used;
        Sigma_hat_var   = U_hat_var*U_hat_var'/(T-k*p);
    else
        assert(size(Z,2) == T_used);
        assert(size(Z,1) == k*p+const+trend+trend2);
    
        A_hat_var       = Y_used*Z'/(Z*Z');
        U_hat_var       = Y_used-A_hat_var*Z;
        Sigma_hat_var   = U_hat_var*U_hat_var'/(T-k*p);
    end
end