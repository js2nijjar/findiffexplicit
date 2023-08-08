function oPrice = finDiffExplicit(X,S0,r,sig,Svec,tvec,oType)

% Implementation of the code follows closely from the original author Phil Goddard

% Inputs
% X - strike
% S0 - stock price
% r - risk free interest rate
% sig - volatility
% Svec - vector of stock prices
% tvec - vector of times
% oType - 'CALL' or 'PUT'

% Output
% oPrice - option price

% number of intervals
J = length(Svec)-1;
N = length(tvec)-1;

% length of time interval
dt = tvec(2)-tvec(1);

% coefficients aj, bj, cj
j = 1:J-1;
sig2 = sig*sig;
j2 = j.*j;
aj = 0.5*dt*(sig2*j2-r*j);
bj = 1-dt*(sig2*j2+r);
cj = 0.5*dt*(sig2*j2+r*j);

% pre-allocate output
price(1:J+1,1:N+1) = nan;

% boundary conditions
switch oType
    case 'CALL'
        price(:,end) = max(Svec-X,0);
        price(1,:) = 0;
        price(end,:) = (Svec(end)-X)*exp(-r*tvec(end:-1:1));
    case 'PUT'
        price(:,end) = max(X-Svec,0);
        price(1,:) = (X-Svec(end))*exp(-r*tvec(end:-1:1));
        price(end,:) = 0;
end

% tridiagonal matrix
A = diag(bj);
A(2:J:end) = aj(2:end); % terms below the diagonal
A(J:J:end) = cj(1:end-1); % terms above the diagonal

% constants not included in A
K = [aj(1); cj(end)];

% price at interior nodes
for i = N:-1:1
    price(2:end-1,i) = A*price(2:end-1,i+1);
    price([2 end-1],i) = price([2 end-1],i) + K.*price([1 end],i+1);
end

% option price
oPrice = interp1(Svec,price(:,1),S0);
