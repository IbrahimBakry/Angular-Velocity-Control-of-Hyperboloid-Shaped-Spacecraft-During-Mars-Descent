% Here is an example of using fitnlm(). For simplicity, none of
% of the fitted parameters are actually nonlinear!
% Define the data to be fit
x = [0 5 10 15 20]'; % Explanatory variable
y = [0 -0.02 -0.04 -0.05 -0.06]';
% Tabulate the data
tbl = table(x,y);
% %  Fit the model
% Define function that will be used to fit data
% (F is a vector of fitting parameters)
f = @(F,x) F(1).*sind(x);
beta0 = [1];
mdl = fitnlm(tbl,f,beta0);
% Calculate the model values at the empirical x
y_predicted = predict(mdl,x);
% Plot the data and fit
figure
plot(x,y,'*',x,y_predicted,'g');
legend('data','fit')