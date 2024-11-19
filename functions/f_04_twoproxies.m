function [out] = f_04_twoproxies(theta,U_hat,z1_used,z2_used)
[n,T] = size(U_hat);

theta1 = theta(1:n-2);
theta2 = theta(n-2+1:end);

out = [(U_hat(3:end,:)  -theta1*U_hat(1,:)    -theta2*U_hat(2,:)    )   .*z1_used;...
       (U_hat(3:end,:)  -theta1*U_hat(1,:)    -theta2*U_hat(2,:)    )   .*z1_used.^2;...
       (U_hat(3:end,:)  -theta1*U_hat(1,:)    -theta2*U_hat(2,:)    )   .*z2_used;...
       (U_hat(3:end,:)  -theta1*U_hat(1,:)    -theta2*U_hat(2,:)    )   .*z2_used.^2];

end