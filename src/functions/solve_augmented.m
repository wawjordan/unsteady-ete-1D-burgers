function p = solve_augmented(A,y,weight,ind)
n = length(y);
% m = size(A,2);

w = eye(n);
w(ind,ind) = weight;

p = (A'*w*A)\(A'*w*y);

% w = eye(n);
% w(n,n) = 1/weight;
% 
% augA = [w,A;A',zeros(m)];
% augb = [y;zeros(m,1)];
% 
% rx = augA\augb;
% 
% p = rx(n+1:n+m);

% [Q,R] = qr(A);
% 
% RHS = Q'*y;
% 
% RHS2 = RHS(1:n);
% 
% p = R\RHS2;


end