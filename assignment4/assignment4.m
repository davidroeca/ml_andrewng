function out = cost_function(theta, x, y, m)
  s = 0; % Start the sum
  for i=1:m
    t_dot_x = theta' * x(i, :)';
    s += -y(i, :) .* log(g(t_dot_x)) - (1 - y(i, :)) .* log(1 - g(t_dot_x));
  end
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
    s += g(t_dot_x) .* (1 - g(t_dot_x)) * x(i, :)' * x(i, :);
  end
  out = s ./ m;
end

function newtheta = newton_step(theta, x, y, m)
  newtheta = theta - inv(hessian(theta, x, m)) * gradj(theta, x, y, m);
end

function finaltheta = newtons_method(theta, x, y, m, iters)
  next_theta = theta;
  for i=1:iters
    next_theta = newton_step(next_theta, x, y, m);
  end
  finaltheta = next_theta;
end

function v = g(z)
  v = 1.0 ./ (1.0 + exp(-z));
end

function main_run(script = false)
  y = load('ex4y.dat');
  m = length(y);
  x = [ones(m, 1), load('ex4x.dat')];

  n = length(x(1,:));
  theta_guess = zeros(n, 1);
  iters = 7; % Plotted previously--converges at around 4 iterations

  theta = newtons_method(theta_guess, x, y, m, iters);
  p_reject = 1 - g(theta' * vec([1., 20., 80.]));
  printf('That student would have a %f%% chance of rejection\n', p_reject * 100);

  if (!script)
    pos = find(y == 1);
    neg = find(y == 0);
    plot(x(pos, 2), x(pos, 3), '+');
    hold on;
    plot(x(neg, 2), x(neg, 3), 'o');
    x_rng = min(x(:,2)):max(x(:,2));
    y_rng = (-theta(1) - x_rng * theta(2))./theta(3)
    plot(x_rng, y_rng);
  end
end
