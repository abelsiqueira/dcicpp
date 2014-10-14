f = @(x) 0.01*x(1)^2 + x(2)^2;
c = @(x) [x(1)*x(2) - 25; x(1)^2 + x(2)^2 - 25];
bl = [2; 0];
bu = [50; 50];
cl = [0;0];
cu = [inf; inf];

jacob = @(x) [x(2), x(1); 2*x(1), 2*x(2)];

jacobz = @(x,s) [jacob(x), -eye(2)];
h = @(x,s) c(x) - s;

Diag = @(x,s) diag([
    (bu(1) - x(1))*(x(1) - bl(1));
    (bu(2) - x(2))*(x(2) - bl(2));
    s(1);
    s(2)]);
A = @(x,s) jacobz(x,s)*Diag(x,s);

q = @(x, y, mu, d) 0.5 * dot(d, B(x,y,mu)*d) + dot(d, g(x,mu));

% x = [2; 2.9669];
% s = [0.991267; 0.997745];
% x = [3.07266; 3.94338];
% s = [0.991093; 0.997836];
% x = [3.306544; 4.809837];
% s = [0.991271; 0.997744];
x = [2.48; 2];
s = [1; 1];
xc = [3.52801; 4.66875];
sc = [0.999862; 0.999942];
mu = 1;
solx = [sqrt(250); sqrt(2.5)];
soly = -0.2;


d = -A(x,s)'*h(x,s);
dn = -A(x,s)'*inv(A(x,s)*A(x,s)')*h(x,s);
d2 = -jacobz(x,s)'*h(x,s);
d2n = -jacobz(x,s)'*inv(jacobz(x,s)*jacobz(x,s)')*h(x,s);

% hold off
% subplot(2,2,1)
% fplot( @(t) norm(h(x + t*d(1:2), s + t*d(3:4))), [0,1],'k');
% subplot(2,2,2)
% fplot( @(t) norm(h(x + t*dn(1:2), s + t*dn(3:4))), [0,1],'r');
% subplot(2,2,3)
% fplot( @(t) norm(h(x + t*d2(1:2), s + t*d2(3:4))), [0,1],'g');
% subplot(2,2,4)
% fplot( @(t) norm(h(x + t*d2n(1:2), s + t*d2n(3:4))), [0,1],'b');
