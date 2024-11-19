function K = zFC_CommutationMatrix(m, n)
%     [m, n] = size(A);
    I = reshape(1:m*n, [m, n]); % initialize a matrix of indices of size(A)
    I = I'; % Transpose it
    I = I(:); % vectorize the required indices
    K = eye(m*n); % Initialize an identity matrix
    K = K(I,:); % Re-arrange the rows of the identity matrix

%     assert(sum(vec(A') == K*vec(A) ) == m*n)
end