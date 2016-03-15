function newtheta = descend(theta, alpha, x, y, m)
  newtheta = theta - alpha ./ m * x' * (x * theta - y);
end

function out = cost_func(x, theta, y, m)
  out = 0.5 / m .* (x * theta - y)'*(x * theta - y);
end

function j_vec = descend_n_times(n, theta, alpha, x, y, m)
  for i=1:n
    theta = descend(theta, alpha, x, y, m);
    j_vec(i, :) = cost_func(x, theta, y, m)';
  end
end

function theta = steepest_descent(theta, alpha, x, y, m, threshold = 10e-3)
  j_now = cost_func(x, theta, y, m);
  do
    j_prior = j_now;
    theta = descend(theta, alpha, x, y, m);
    j_now = cost_func(x, theta, y, m);
  until (abs(j_now - j_prior) < threshold)
end

function part_1_plots(n, theta, alpha, x, y, m, script)
  J = descend_n_times(n, theta, alpha, x, y, m);
  if (!script)
    xlabel('Number of iterations');
    ylabel('Cost J');
    hold on;
    alpha = 0.03;
    p = 'b-';
    J = descend_n_times(n, theta, alpha, x, y, m);
    plot(1:n, J(1:n), p);
    alpha = 0.1;
    p = 'r-';
    J = descend_n_times(n, theta, alpha, x, y, m);
    plot(1:n, J(1:n), p);
    alpha = 0.3;
    p = 'k-';
    J = descend_n_times(n, theta, alpha, x, y, m);
    plot(1:n, J(1:n), p);
  end
end

function theta_norm = normal_eq(x, y)
  theta_norm = inv(x' * x) * x' * y;
end
  
function main_run(script=false, n=100)
  y = load('ex3y.dat'); 
  m = length(y);
  x = [ones(m, 1), load('ex3x.dat')]; 

  sigma = std(x);
  mu = mean(x);
  x(:, 2) = (x(:,2) - mu(2))./ sigma(2);
  x(:, 3) = (x(:,3) - mu(3))./ sigma(3);
  %-------------------------------------------------------------
  % Part 1
  %-------------------------------------------------------------
  theta = zeros(size(x(1,:)))';
  alpha = 0.01;
  part_1_plots(n, theta, alpha, x, y, m, script);
  theta = steepest_descent(theta, alpha, x, y, m)
  price = theta(1);
  val = [1, 1650, 3]
  for i=2:3
    price += (val(i) - mu(i)) / sigma(i) * theta(i);
    printf('%f * %f\n', val(i), theta(i));
  end
  printf('price = %f\n', price);
  %-------------------------------------------------------------
  % Normal Equations
  %-------------------------------------------------------------

  x_norm = [ones(m, 1), load('ex3x.dat')]; 
  theta_norm = normal_eq(x_norm, y);
  price_norm = 0.0;
  for i=1:3
    price_norm += val(i) * theta_norm(i);
    printf('%f * %f\n', val(i), theta_norm(i));
  end
  printf('price_norm = %f\n', price_norm);
end
