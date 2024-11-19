function [out] = getf_MR(b,u,z)
    e =   ( inv(getB_MR(b))  * u' )';
    
    ez1 = [   z.* e(:,2:end) ] ; % proxy condition
    
    e2 = e.^2-1; % variance one condition
    for i = 1: size(e,2)-1
        e2 = [e2 , e(:,i).*e(:,i+1:end)]; % add covariance conditions
    end
     
    
    out = [ ez1,  e2]';
    
    end