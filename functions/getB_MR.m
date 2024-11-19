function [B] = getB_MR(b)
n = sqrt(size(b,1));
B=reshape(b,n,n);
B(2,3)=0;
end

