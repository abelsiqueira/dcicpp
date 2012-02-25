function plotcontour (fun, x, y)

M = numel(x);
N = numel(y);

Z = zeros (M,N);

for i = 1:M
  for j = 1:N
    Z(i,j) = feval (fun, [x(i); y(j)]);
  end
end

[X,Y] = meshgrid(x, y);
hold off
contour (X, Y, Z)
hold on
