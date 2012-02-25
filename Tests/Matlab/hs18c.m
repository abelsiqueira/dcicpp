f = @(x) 0.01*x(1)^2 + x(2)^2;
c = @(x) [x(1)*x(2) - 25];
bl = [2; 0];
bu = [50, inf];
cl = [0];
cu = [inf];

jacob = @(x) [x(2), x(1)];
jacobz = @(x,s) [jacob(x), -1];
h = @(x,s) c(x) - s;
Diag = @(x,s) diag([(bu(1) - x(1))*(x(1) - bl(1)); x(2) - bl(2); s]);
A = @(x,s) jacobz(x,s)*Diag(x,s);
% hess = @(x,y,mu) [0.02 + mu/(x(1) - 2)^2, y; y, 2 + mu/x(2)^2];
hess = @(x,s,y,mu) [0.02, y, 0; y, 2, 0; 0, 0, 0] + mu*diag([1/(x(1) - bl(1))^2; 1/(bu(1) - x(1))^2; 1/(x(2) - bl(2))^2; 1/(s - cl(1))^2]);
% hess = @(x,y,mu) [0.02, y; y, 2];
B = @(x,s,y,mu) Diag(x,s)*hess(x,y,mu)*Diag(x,s);
% B = @(x,y,mu) Diag(x)*hess(x,y,mu)*Diag(x) + mu*eye(2);
% g = @(x,mu) Diag(x) * [0.02 * x(1) - mu/(x(1) - 2); 2*x(2) - mu/x(2)];
g = @(x,s,mu) Diag(x,s) * ([0.02 * x(1); 2*x(2); 0] - mu*[1/(x(1) - bl(1)) - 1/(bu(1) - x(1)); 1/(x(2) - bl(2)); 1/(s - cl(1))]);
gp = @(x,s,y,mu) g(x,s,mu) + A(x,s)'*y;
lls = @(x,s,mu) -inv(A(x,s)*A(x,s)')*A(x,s)*g(x,s,mu);

% q = @(x, y, mu, d) 0.5 * dot(d, B(x,y,mu)*d) + dot(d, g(x,mu));

x = [15.8133; 1.58087];
s = 6.12225e-7;
y = -0.200057;
mu = 1.87908e-5;
solx = [sqrt(250); sqrt(2.5)];
sols = 0;
soly = -0.2;
% x = [2.48; 2];

% d = -A(x)'*h(x);
% dn = -A(x)'*inv(A(x)*A(x)')*h(x);
% d2 = -jacob(x)'*h(x);
% d2n = -jacob(x)'*inv(jacob(x)*jacob(x)')*h(x);

% hold off
% subplot(2,2,1)
% fplot( @(t) norm(h(x + t*d)), [0,1],'k');
% subplot(2,2,2)
% fplot( @(t) norm(h(x + t*dn)), [0,1],'r');
% subplot(2,2,3)
% fplot( @(t) norm(h(x + t*d2)), [0,1],'g');
% subplot(2,2,4)
% fplot( @(t) norm(h(x + t*d2n)), [0,1],'b');
