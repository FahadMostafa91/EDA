function [C, Large_crit, Perm_crit]=EqualCovtest(X11,X22,alpha,B1)

[n1,p]=size(X11);
[n2,p]=size(X22);

% Calculate Test Statistic
S1=cov(X11);
S2=cov(X22);

Sp=(n1-1)/(n1+n2-2)*S1+(n2-1)/(n1+n2-2)*S2;

Lam=(det(S1)/det(Sp))^((n1-1)/2)*(det(S2)/det(Sp))^((n2-1)/2);
M=-2*log(Lam);


% Calculate Large Sample Critical Value
u1=1/(n1-1)+1/(n2-1)- 1/(n1-1+n2-1);
u2=(2*p^2+3*p-1)/(6*(p+1)*(2-1));
u=u1*u2;

C=(1-u)*M;

nu=.5*p*(p+1)*(2-1);

Large_crit=chi2inv(0.95,nu);


% Obtain Permutation Critical Value
D1=X11-mean(X11);
D2=X22-mean(X22);

D=[D1;D2];
Cp=zeros(B1,1);

for ind=1:B1
    permind=randperm(n1+n2);

    Dp1=D(permind(1:n1),:);
    Dp2=D(permind((n1+1):(n1+n2)),:);
    
    Sp1=cov(Dp1);
    Sp2=cov(Dp2);
    
    Spp=(n1-1)/(n1+n2-2)*Sp1+(n2-1)/(n1+n2-2)*Sp2;
    
    Lamp=(det(Sp1)/det(Spp))^((n1-1)/2)*(det(Sp2)/det(Spp))^((n2-1)/2);
    Cp(ind)=-2*log(Lamp)*(1-u);
end

Perm_crit=quantile(Cp,1-alpha);



