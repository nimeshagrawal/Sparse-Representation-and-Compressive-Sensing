n = input('Please enter value of signal length(n):');   %n = size of input signal
p = input('Please enter value of probability(p):');     %p = probability of success
e = input('Please enter value of threshold(e):');        %e = thresold for residual error
%input signal generator
A = seq(n,p);                             
X = double.empty(n,0);                   

for i = 1:length(A)
    if A(i) == 1
        X(i,1) = randn;
    else
        X(i,1) = 0;
    end
end

%Sparsity indicator
h = 0:n;
q = binopdf(h,n,p);
figure
bar(h,q,1)
xlabel("Sparsity level")
ylabel("Probability")
title("N = " + n + " p = " + p)


error = [];                     %recovery error initialization
K_hat = [];                     %Sparsity of recovered signals
for m = [200,400, 600]          %m = number of measurements
% Inputs
K = sum(A);                     %sparsity of input signal
phi = randn(m,n);               %random matrix
y = phi*X;                      %measurement vector


%OMP recovery 
%Initialization
k = 0;                          %iteration count
Xhat = zeros(n,1);              %estimate of x
yhat = zeros(m,1);              %Approximation of y
rhat = y;                       %residual
lambda = [];                    %Suppport of xhat
%loop
 while norm(rhat) > e
    k = k + 1;
    
    %sweep
    [value, location] = max(abs(phi'*rhat));
    
    %update support
    
    lambda = [lambda, location];
    
    %update provisional solution
    %Selecting columns corresponding to updated support and finding
    %least square solution
    phi_k = phi(:,lambda);
    x_k = lsqminnorm(phi_k, y);
    for i = 1:length(lambda)
        Xhat(lambda(i)) = x_k(i);               %Updating xhat with least square solution with 
    end                                         %correspoding support
  
    %update residual
    yhat = phi*Xhat;
    rhat = y - yhat;
 end
%Comparison plot for original and recovered signal
figure;hold on;grid on;
stem(1:length(A),X,'r');
stem(1:length(A),Xhat,'b');
legend('X','Xhat')
title("N = " + n + " M = " + m + " p = " + p)
error = [error, norm(X - Xhat)];           %Reccovery error
K_hat = [K_hat, nnz(Xhat)];                %Sparsity for recovered signals
end

% Display values
error
K_hat
K
%Function for bernoulli sequence generator
function A = seq(n,p)
    A = (rand(n,1) < p);
end

%Note: After running the code, it will prompt the user to input certain
%values. After running through the given number of measurements set i.e [200 400 600], it
%will display 4 plots(3 plots of comparison and one for binomial distribution 
%to check typical range of sparsity) which are overlapping over each other so the user
%need to drag the images. Also the display values(error, K_hat and K) can be seen in the
%command window. For higher sparsity level computation time is quite large.
