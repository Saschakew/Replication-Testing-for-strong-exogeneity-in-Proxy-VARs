function [out] = f_05(all_param,Y,z_used,p,const,trend,trend2)
[K,T] = size(Y);

T_used = T-p;

assert(length(all_param) == K*(K*p+const+trend+trend2) + (K-1));

a_hat = all_param(1:K*(K*p+const+trend+trend2));
A_hat = reshape(a_hat,K,K*p+const+trend+trend2); 
theta = all_param(K*(K*p+const+trend+trend2)+1:end);

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

out_var = [];
for kk = 1:K
    out_var = [out_var; Z.*(Y(kk,p+1:end) - A_hat(kk,:)*Z)];
end

U_hat = Y(:,p+1:end) - A_hat*Z;

out_proxy = [(U_hat(2:end,:)-theta*U_hat(1,:)).*z_used;...
       (U_hat(2:end,:)-theta*U_hat(1,:)).*z_used.^2];

out = [out_var;...
       out_proxy];

end