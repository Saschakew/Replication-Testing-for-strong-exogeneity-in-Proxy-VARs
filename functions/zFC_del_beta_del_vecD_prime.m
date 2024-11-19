function del_beta_del_vecD_prime = zFC_del_beta_del_vecD_prime(D_1,D_2,Q,K,K_1,K_2,N)

% D_1 = D_hat_m1_1;
% D_2 = D_hat_m1_2;
% Q   = Q_m1;
% K_1
% K_2
% N
% Compute del_beta_del_vecDprime DO write function for general m, CAN IT REALLY BE eye(K_1)???
del_vecD1QD1prime_del_vecD1prime = kron(D_1*Q,eye(K_1)) ...
                                + kron(eye(K_1),D_1*Q)*zFC_CommutationMatrix(K_1,N);
assert(size(del_vecD1QD1prime_del_vecD1prime,1) == K_1^2)                            
assert(size(del_vecD1QD1prime_del_vecD1prime,2) == K_1*N)                            

del_vecD2QD1prime_del_vecD1prime = kron(eye(K_1), D_2*Q)*zFC_CommutationMatrix(K_1,N);
assert(size(del_vecD2QD1prime_del_vecD1prime,1) == K_1*(K-K_1))                            
assert(size(del_vecD2QD1prime_del_vecD1prime,2) == K_1*N)                            

del_vecD2QD1primeD1QD1primeminus1_del_vecD2prime = kron(inv(D_1*Q*D_1')*D_1*Q,eye(K_2));
assert(size(del_vecD2QD1primeD1QD1primeminus1_del_vecD2prime,1) == K_1*(K-K_1))                            
assert(size(del_vecD2QD1primeD1QD1primeminus1_del_vecD2prime,2) == (K-K_1)*N)                            

del_vecD1vecD2_del_vecD1primevecD2prime = [zFC_CommutationMatrix(N,K_1) zeros(K_1*N,(K-K_1)*N);...
                                           zeros((K-K_1)*N,K_1*N) zFC_CommutationMatrix(N,K-K_1)];
assert(size(del_vecD1vecD2_del_vecD1primevecD2prime,1) == K*N)                            
assert(size(del_vecD1vecD2_del_vecD1primevecD2prime,2) == K*N)                            

del_vecD1primevecD2prime_del_vecDprime = zFC_CommutationMatrix(K,N);
assert(size(del_vecD1primevecD2prime_del_vecDprime,1) == K*N)                            
assert(size(del_vecD1primevecD2prime_del_vecDprime,2) == K*N)                            

% Now build the remaining derivatives out of these five
del_vecD2QD1primeD1QD1primeminus1_del_vecD1prime = kron(inv(D_1*Q*D_1'),eye(K-K_1))...
                                                    *del_vecD2QD1prime_del_vecD1prime...
                                                   -kron(inv(D_1*Q*D_1'),D_2*Q*D_1'*inv(D_1*Q*D_1'))...
                                                   *del_vecD1QD1prime_del_vecD1prime;   

del_vec_D2QD1primeD1QD1primeminus1_del_vecD1primevecD2prime = [del_vecD2QD1primeD1QD1primeminus1_del_vecD1prime del_vecD2QD1primeD1QD1primeminus1_del_vecD2prime];

del_beta_del_vecD_prime = del_vec_D2QD1primeD1QD1primeminus1_del_vecD1primevecD2prime*del_vecD1vecD2_del_vecD1primevecD2prime*del_vecD1primevecD2prime_del_vecDprime;



end