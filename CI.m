%hypothesis test for mean of cancer B and cancer M

p=9; % number of covariates
%conB = table2array(conB(:, 3:32));
%conM = table2array(conM(:, 3:32));
n1=length(X11); % sample size for expHealthy
n2=length(X22); % sample size for expPatient
mu_he=mean(X11);
mu_pt=mean(X22);
mu_diff=mu_he-mu_pt;
S1=cov(X11);
S2=cov(X22);
S=(1/n1)*S1+(1/n2)*S2;
invS=inv(S);
conT2=(mu_diff)*invS*(mu_diff)'; % T^2 test statistics
% % calcualting critical value
nu1=(1/n1)*((trace((1/n1)*S1*invS)^2)+(trace((1/n2)*S1*invS)^2));
nu2=(1/n2)*((trace((1/n1)*S2*invS)^2)+(trace((1/n2)*S2*invS)^2));
nu=(p+p^2)/(nu1+nu2);
crtcon=(nu*p)/(nu-p+1)*finv(1-alpha,p,nu-p+1);
%
% % confidence intervals
concrt=chi2inv(0.95,p);
for i=1:p
 conlow(i)=mu_diff(i)-sqrt(concrt*S(i,i));
 conupper(i)=mu_diff(i)+sqrt(concrt*S(i,i));
end
%