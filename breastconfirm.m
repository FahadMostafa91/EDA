% Taking only the quantitative variables
B=breast(:,1:10);
confind=setdiff(1:116,expind);
C=breast(confind,:);
% Extract data and save it as a new name:X
X1=B(confind,1:9);
Y1=B(confind,:);
% Plot the sample data; https://www.mathworks.com/help/stats/gplotmatrix.html
figure()
gplotmatrix(X1)
S1=cov(X1);
imagesc(S1)
% Correlation matrix for new data 
R1 = corrcoef(X1);
imagesc(R1)
imagesc(abs(R1))
%Plot the correlation matrix; https://www.mathworks.com/help/econ/corrplot.html
hold on;
corrplot(R1)
% Find the eigenvalue and eigenvector of X
 [eigenvector1, eigenvalue1]=eig(R1);
 % Make the diagonal elements into decending order
 s2=sort(diag(eigenvalue1),'descend');
 %Total Variation
 TVP1=trace(eigenvalue1);
 %Proportion of variation
 var1=s2/TVP1;
 % to see variation we can use scree plot
  plot(var1)
  title('scree-plot')
   %Cumulative sum property
 cp1=cumsum(s2/TVP1);
% Plotting variation
figure()
pareto(var1)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

% To show the Principal direction ## https://www.mathworks.com/help/stats/biplot.html
figure()
Z1 = zscore(X1); % Standardized data
[coefs1,score1] = pca(Z1);
biplot(coefs1(:,1:3))
vbls = {'MCP.1','Resistin','Adiponectin','Leptin','HOMA','Insulin','Glucose','BMI','Age'}; % Labels for the variables
biplot(coefs1(:,1:3),'Scores',score1(:,1:3),'Color','b','Marker','o','VarLabels',vbls);


Y1=B(confind,:);
response1=Y1(:,10);
Hind1=find(response1==1);
HC1=X1(Hind1,:);
Pind1=find(response1==2);
PT1=X1(Pind1,:);
%%%Mardia test for checking Multivariate Normality
alpha=0.05; %define alpha
[H1 stats1] = mardiatest(HC1, alpha) 
[H1 stats1] = mardiatest(PT1, alpha) 

% Equal Covarience Test
X11=HC1;X22=PT1;alpha=0.05;B1=20000;p=9; 

