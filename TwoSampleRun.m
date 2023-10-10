mu1=[1 3 5]';
mu2=[1 3 8]';

Sig1=[3 1.5 .2;
     1.5 2  .1;
     .2 .1 4];
 
Sig2=[4 -1.5 .2;
     -1.5 3  .1;
     .2 .1 3];
 
 d=[0 0 0]';
 
 alpha=.05;

n1=20;
n2=29;
p=3;
 
X1=mvnrnd(mu1,Sig1,n1);
X2=mvnrnd(mu2,Sig2,n2);

B=20000;

%% Test for Equal Means

[T2pool, T2, Pool_crit, Unpool_crit, Large_crit, Perm_crit, Boot_crit]=TwoSampleT2test(X1,X2,d,alpha,B)

%% Test for Equal Covariances

[M, Large_crit, Perm_crit]=EqualCovtest(X1,X2,alpha,B)

