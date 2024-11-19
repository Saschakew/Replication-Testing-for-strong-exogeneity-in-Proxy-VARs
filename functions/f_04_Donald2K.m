function [out] = f_04_Donald2K(theta,U_hat,z_used)
[n,T] = size(U_hat);

  

% Compute K as the cube root of T
K = round(T^(1/3));

% Check if K is even
if mod(K, 2) == 0
    % If K is even, subtract 1 to make it odd
   K = K - 1;
end

percentiles = (1:K) * 100 / (K + 1);
t_vec = prctile(z_used, percentiles);

q_K = [ones(1,T);...
        z_used;...
        z_used.^2;...
        z_used.^3;  
        ];

step2 = z_used' - t_vec;
step2(step2<0) = 0;
step2=step2.^3 ; 
step2=step2';
q_K = [q_K; step2]; 

Udiff = (U_hat(2:end,:)-theta*U_hat(1,:));
num_rows = size(q_K, 1);  
num_rows2 = size(Udiff, 1);  
q_K_tiled = repmat(q_K, num_rows2, 1);
Udiff_tiled = repmat(Udiff, num_rows, 1);
out = Udiff_tiled .* q_K_tiled;

end