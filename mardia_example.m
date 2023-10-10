% Simulate Data from Null distribution
n=10;
p=5;
mu=ones(p,1);
Sigma=eye(p);
alpha=.05;

x=mvnrnd(mu,Sigma,n);

[H stats] = mardiatest(x,alpha)


%% Simulate Data from non-Null distribution
n=1000;
p=5;
mu=ones(p,1);
Sigma=eye(p);
alpha=.05;

y=trnd(1,n,p);

[H stats] = mardiatest(y,alpha)