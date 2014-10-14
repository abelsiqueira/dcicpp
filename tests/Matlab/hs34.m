f = @(x) -x(1);
c = @(x) [x(2) - exp(x(1)); x(3) - exp(x(2))];
bl = [0; 0; 0];
bu = [100; 100; 10];
cl = [0; 0];
cu = [inf; inf];

jacob = @(x) [-exp(x(1)), 1, 0; 0, -exp(x(2)), 1];
jacobz = @(x,s) [jacob(x), -eye(2)];
h = @(x,s) c(x) - s;
Diag = @(x,s) diag([(x(1) - bl(1))*(bu(1) - x(1));
                    (x(2) - bl(2))*(bu(2) - x(2));
                    (x(3) - bl(3))*(bu(3) - x(3));
                    s(1) - cl(1);
                    s(2) - cl(2) ]);
A = @(x,s) jacobz(x,s)*Diag(x,s);
hess = @(x,s,y,mu) y(1)*diag([-exp(x(1));0;0;0;0]) + y(2)*diag([0;-exp(x(2));0;0;0]) + mu*diag([1/(x(1)-bl(1)) + 1/(bu(1) - x(1)); 1/(x(2)-bl(2)) + 1/(bu(2) - x(2)); 1/(x(3)-bl(3)) + 1/(bu(3) - x(3)); 1/(s(1) - cl(1)); 1/(s(2) - cl(2))]);
B = @(x,s,y,mu) Diag(x,s)*hess(x,s,y,mu)*Diag(x,s);

g = @(x,s,mu) Diag(x,s) * ([-1;0;0;0;0] - mu* [1/(x(1) - bl(1)) - 1/(bu(1) - x(1)); 1/(x(2) - bl(2)) - 1/(bu(2) - x(2)); 1/(x(3) - bl(3)) - 1/(bu(3) - x(3)); 1/(s(1) - cl(1)); 1/(s(2) - cl(2))]);

gp = @(x,s,y,mu) g(x,s,mu) + A(x,s)'*y;
lls = @(x,s,mu) -inv(A(x,s)*A(x,s)')*A(x,s)*g(x,s,mu);

x = [3.92e-8; 1.44401; 7.78085];
s = [1.5879; 3.85295];
y = [-0.664936; -0.317939];

solx = [log(log(10)); log(10); 10];
sols = [0; 0];
soly = 1;
