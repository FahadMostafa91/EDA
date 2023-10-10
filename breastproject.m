% Taking only the quantitative variables
B=breast(:,1:9);
%B= table2array(breast);
%Random sampling of data set
expind = randsample(116,39)';
% Extract data and save it as a new name:X
X=B(expind,:);
% Plot the sample data; https://www.mathworks.com/help/stats/gplotmatrix.html
gplotmatrix(X)
boxplot(X,'PlotStyle','compact')
%To check normality of data set
probplot('normal',X)
% Correlation matrix for new data 
R = corrcoef(X);
%Plot the correlation matrix; https://www.mathworks.com/help/econ/corrplot.html
corrplot(R)
% Find the eigenvalue and eigenvector of X
 [eigenvector eigenvalue]=eig(R);
 % Make the diagonal elements into decending order
 s1=sort(diag(eigenvalue),'descend');
 %Total Variation
 TVP=trace(eigenvalue);
 %Proportion of variation
 var=s1/TVP;
 % to see variation we can use scree plot
  plot(var)
  title('scree-plot')
 %Cumulative sum property
 cp=cumsum(s1/TVP)
% Plotting eigenvalues and eigenvector

