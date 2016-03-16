function out = cost_function(theta, x, y, m)
  s = 0; % Start the sum
  for i=1:m
    t_dot_x = theta' * x(i, :)';
    s += -y(i, :) * log(t_dot_x) - (1 - y(i, :)) * log(1 - g(t_dot_x));
  end
  vall = size(s)
  printf('%i, %i\n', vall(1), vall(2));
  out = 1 ./ m .* s;
end

function out = gradj(theta, x, y, m)
  n_1 = length(x(1,:));
  s = zeros(n_1, 1);
  for i=1:m
    t_dot_x = theta' * x(i, :)';
    s += (g(t_dot_x) - y(i, :)) * x(i, :)';
  end
  out = s ./ m;
end

function out = hessian(theta, x, m)
  n_1 = length(x(1,:));
  s = zeros(n_1, n_1);
  for i=1:m
    t_dot_x = theta' * x(i, :)';
    s += g(t_dot_x) * (1 - g(t_dot_x)) * x(i, :)' * x(i, :);
  end
  out = s ./ m;
end

function newtheta = newton_step(theta, x, y, m)
  newtheta = theta - inv(hessian(theta, x, m)) * gradj(theta, x, y, m);
end

function v = g(z)
  v = 1.0 ./ (1.0 + exp(-z));
end

function main_run(script = false)
  y = load('ex4y.dat');
  m = length(y);
  x = [ones(m, 1), load('ex4x.dat')];
  if (!script)
    pos = find(y == 1);
    neg = find(y == 0);
    plot(x(pos, 2), x(pos, 3), '+');
    hold on;
    plot(x(neg, 2), x(neg, 3), 'o');
  end
end
