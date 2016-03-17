source('map_feature.m');

function out = regular_matrix(x, lambda);
  out = eye(length(x(1,:)));
  out(1,1) = 0;
  out = lambda * out;
end

function out = normal_solution(x, y, lambda)
  out = inv(x' * x + regular_matrix(x, lambda)) * x' * y;
end

function plot_normal_soln(x, y, n, lambda, fmt='-')
  theta = normal_solution(x, y, lambda);
  printf('lambda: %i\n', lambda)
  printf('theta: [\n');
  for i=1:(n+1)
    printf('\t%f\n', theta(i));
  end
  printf(']\n');
  fn_plot = @(v) theta(1) + theta(2)*v + theta(3)*(v.^2) + ...
                 theta(4)*(v.^3) + theta(5)*(v.^4) + theta(6)*(v.^5);
  fplot(fn_plot, [0,1], fmt);
end

function main_run()
  y = load('ex5Liny.dat');
  m = length(y);
  x = [load('ex5Linx.dat')];
  plot(x, y, 'o');
  hold on;
  x = [ones(m, 1), x, x .^ 2, x .^ 3, x .^ 4, x .^ 5];
  n = length(x(1, :)) - 1;
  printf('%i\n', n);

  plot_normal_soln(x, y, n, 0, '-r');
  plot_normal_soln(x, y, n, 1, '-g');
  plot_normal_soln(x, y, n, 10, '-b');
end
