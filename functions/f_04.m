function [out] = f_04(theta,U_hat,z_used)
[n,T] = size(U_hat);

% theta_hat_z = U_hat(2:end,:)*z_used'./(U_hat(1,:)*z_used');
% theta_hat_z_sq = U_hat(2:end,:)*z_used.^2'./(U_hat(1,:)*z_used.^2');
% 
% out = [theta - theta_hat_z;...
%        theta - theta_hat_z_sq];

out = [(U_hat(2:end,:)-theta*U_hat(1,:)).*z_used;...
       (U_hat(2:end,:)-theta*U_hat(1,:)).*z_used.^2];

end