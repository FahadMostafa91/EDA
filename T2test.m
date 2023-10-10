function [T2, crit_chi, crit_T2, crit_boot]=T2test(X,mu,alpha,B)

[n, p]=size(X);

xbar=mean(X)';
S=cov(X);

T2=n*(xbar-mu)'*S^-1*(xbar-mu);

crit_chi=chi2inv(1-alpha,p);

crit_T2=(n-1)*p/(n-p)*finv(1-alpha,p,n-p);

crit_boot=99999999;

if B>0
    
    bsamp=randi(n,n,B);
    T2b=zeros(B,1);
    for ind=1:B
        Xb=X(bsamp(:,ind),:);
        xbarb=mean(Xb)';
        Sb=cov(Xb);
        
        T2b(ind)=n*(xbarb-xbar)'*Sb^-1*(xbarb-xbar);
    end
    
    crit_boot=quantile(T2b,1-alpha);
    
end