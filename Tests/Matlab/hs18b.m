f = @(x) 0.01*x(1)^2 + x(2)^2;
c = @(x) [x(1)*x(2) - 25];
bl = [2; 0];
cl = [0];
cu = [inf];

jacob = @(x) [x(2), x(1)];
h = @(x) c(x);
Diag = @(x) diag([x(1) - bl(1); x(2) - bl(2)]);
A = @(x) jacob(x)*Diag(x);
hess = @(x,y,mu) [0.02 + mu/(x(1) - 2)^2, y; y, 2 + mu/x(2)^2];
% hess = @(x,y,mu) [0.02, y; y, 2];
B = @(x,y,mu) Diag(x)*hess(x,y,mu)*Diag(x);
% B = @(x,y,mu) Diag(x)*hess(x,y,mu)*Diag(x) + mu*eye(2);
g = @(x,mu) Diag(x) * [0.02 * x(1) - mu/(x(1) - 2); 2*x(2) - mu/x(2)];
gp = @(x,y,mu) g(x,mu) + A(x)'*y;
lls = @(x,mu) -inv(A(x)*A(x)')*A(x)*g(x,mu);

q = @(x, y, mu, d) 0.5 * dot(d, B(x,y,mu)*d) + dot(d, g(x,mu));

x = [15.81558703864949; 1.5807190684058807];
y = -0.19983015401231807;
mu = 0.0036606906644387096;
solx = [sqrt(250); sqrt(2.5)];
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
