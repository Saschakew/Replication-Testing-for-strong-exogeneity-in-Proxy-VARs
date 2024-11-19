function [out] = f_04_Donald(theta,U_hat,z_used)
[n,T] = size(U_hat);

% theta_hat_z = U_hat(2:end,:)*z_used'./(U_hat(1,:)*z_used');
% theta_hat_z_sq = U_hat(2:end,:)*z_used.^2'./(U_hat(1,:)*z_used.^2');
% 
% out = [theta - theta_hat_z;...
%        theta - theta_hat_z_sq];

% out = [(U_hat(2:end,:)-theta*U_hat(1,:)).*z_used;...
%        (U_hat(2:end,:)-theta*U_hat(1,:)).*z_used.^2];

t_vec = prctile(z_used, [100/4, 200/4, 300/4]);

q_K = [ones(1,T);...
        z_used;...
        z_used.^2;...
        ];

for ii = 1:3
    step = (z_used - t_vec(ii)).^ii;
    step(z_used - t_vec(ii)<0) = 0;
    q_K = [q_K; step]; 
end

out = [];
for tt = 1:T
    out = [out, kron((U_hat(2:end,tt)-theta*U_hat(1,tt)), q_K(:,tt))];
end


end