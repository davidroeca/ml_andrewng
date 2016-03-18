source('map_feature.m');

function v = g(z)
  v = 1.0 ./ (1.0 + exp(-z));
end

function out = cost_function(theta, x, y, lambda, m)
  reg_term = (0.5 * lambda / m) .* theta' * theta;
  vec_1 = ones(m, 1);
  h = g(x * theta);
  out = -1 / m * (y' * log(h) + (vec_1 - y)' * log(vec_1 - h)) + reg_term;
end

function out = gradj(theta, x, y, lambda, m)
  h = g(x * theta);
  reg_term = (lambda / m) .* theta;
  reg_term(1,1) = 0;
  out = (1 / m) .* x' * (h - y) + reg_term;
end

function out = hessian(theta, x, lambda, m)
  h = g(x * theta);
  reg_term = (lambda / m) .* eye(length(x(1,:)));
  reg_term(1,1) = 0;
  out = (1 / m) .* x' * diag(h) * diag(1 - h) * x + reg_term;
end

function out = newton_step(theta, x, y, lambda, m)
  out = theta - inv(hessian(theta, x, lambda, m)) * gradj(theta, x, y, lambda, m);
end

function out = newtons_method(theta, x, y, lambda, m, iter_max=15)
  out = theta;
  for i=1:iter_max
    out = newton_step(out, x, y, lambda, m);
  end
end

function newtons_method_plot_print(theta, x, y, lambda, m, plot_fmt)
  printf('lambda: %f ------------------------------------\n', lambda);
  theta = newtons_method(theta, x, y, lambda, m);
  printf('norm theta: %f\n', norm(theta));
  u = linspace(-1, 1.5, 200);
  v = linspace(-1, 1.5, 200);
  z = zeros(length(u), length(v));
  for i = 1:length(u)
    for j = 1:length(v)
      z(j, i) = map_feature(u(i), v(j)) * theta;
    end
  end
  z = z';
  contour(u, v, z, [0, 0], 'LineWidth', 2, plot_fmt);
end

function part_2(script=false)
  y = load('ex5Logy.dat');
  m = length(y);
  x = load('ex5Logx.dat');

  pos = find(y);
  neg = find(y == 0);

  if (!script)
    figure;
    plot(x(pos, 1), x(pos, 2), '+');
    hold on;
    plot(x(neg, 1), x(neg, 2), 'o');

    x = map_feature(x(:,1), x(:,2));
    n = length(x(1, :));
    theta = zeros(n,1);
    newtons_method_plot_print(theta, x, y, 0, m, '-r')
    newtons_method_plot_print(theta, x, y, 1, m, '-g')
    newtons_method_plot_print(theta, x, y, 10, m, '-b')
  end
end
