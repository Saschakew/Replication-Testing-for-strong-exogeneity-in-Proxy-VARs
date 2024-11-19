clc
clear all
addpath(genpath('./functions'));
addpath(genpath('./data'));

%% Data
load('data.mat')

date = MertensRavnData(:,1);
y = MertensRavnData(:,2:4);
z = MertensRavnData(:,5);

select = z~=0;
%z = z(select);
%z = (z - mean(z))./sqrt(var(z));
%z = (z - mean(z))  ;

%%
skewness(z)
kurtosis(z)
 
%% VAR
lags = 4;

[T,n] = size(y);

Y0 = y(1:lags,:);  % save the first 4 obs as the initial conditions
shortY = y(lags+1:end,:); 

dummy=zeros(T,1);
dummy(102 )=1;
 
tmpY = [Y0(end-lags+1:end,:); shortY];
X = zeros(T-lags,n*lags); 
for i=1:lags
    X(:,(lags-i+1)*n-n+1:(lags-i+1)*n) = tmpY(lags-i+1:end-i,:);% X(:,(i-1)*n+1:i*n)
end 
Xexo = [ ones(T-lags,1) dummy(lags+1:end) (1:T-lags)'  (1:T-lags).^2'];
X = [ X  Xexo ];
Atilde = [X ] \ shortY;
Atilde = Atilde';
 

u = shortY-[X]*Atilde';
z = z(lags+1:end,:);

%%
[T,n] = size(u);

startB=[	0.012 0.001  0.023  
	        0.001  0.024  0 
	        -0.004  0.0020  0.007];
bstart = startB(:); 

options = optimoptions('fminunc', 'Display', 'off', 'OptimalityTolerance', 1e-10,'StepTolerance', 1e-10,'MaxIterations', 5000);
thisloss_noW = @(b) mean(getf_MR(b,u,z),2)' * mean(getf_MR(b,u,z),2);
[best,~,~,~,~,~]= fminunc( thisloss_noW ,bstart,options); 

Best = getB_MR(best);
 
round(Best,3)
