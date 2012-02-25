f = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
c = @(x) [x(2)^2 - x(1); x(1)^2 - x(2)];
bl = [-0.5; -inf];
bu = [0.5; 1];
cl = [0;0];
cu = [inf; inf];

jacob = @(x) [-1, 2*x(2); 2*x(1), -1];

phi = @(x,s,mu) f(x) - mu*( log(s(1)) + log(s(2)) + log(x(1) - blx(1)) + log(bu(1) - x(1)) + log(bu(2) - x(2)) );
jacobz = @(x,s) [jacob(x), -eye(2)];
h = @(x,s) c(x) - s;

Diag = @(x,s) diag([(bu(1) - x(1))*(x(1) - bl(1));
    bu(2) - x(2); s(1); s(2)]);
A = @(x,s) jacobz(x,s)*Diag(x,s);

x = [-0.499999; 0.122719];
s = [3.87e-6; 0.315753];

d = -A(x,s)'*h(x,s);
dn = -A(x,s)'*inv(A(x,s)*A(x,s)')*h(x,s);
d2 = -jacobz(x,s)'*h(x,s);
d2n = -jacobz(x,s)'*inv(jacobz(x,s)*jacobz(x,s)')*h(x,s);

hold off
fplot( @(t) norm(h(x + t*d(1:2), s + t*d(3:4))), [0,1],'k');
hold on
fplot( @(t) norm(h(x + t*dn(1:2), s + t*dn(3:4))), [0,1],'r');
fplot( @(t) norm(h(x + t*d2(1:2), s + t*d2(3:4))), [0,1],'g');
fplot( @(t) norm(h(x + t*d2n(1:2), s + t*d2n(3:4))), [0,1],'b');
