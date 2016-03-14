function out = cost_func(x_mat, theta, y, m)
  out = 1 / (2 * m) * (x_mat * theta - y)'*(x_mat * theta - y)
end

function main_run
  y = load('ex3y.dat'); 
  m = length(y);
  x = [ones(m, 1), load('ex3x.dat')]; 

  sigma = std(x);
  mu = mean(x);
  x(:, 2) = (x(:,2) - mu(2))./ sigma(2);
  x(:, 3) = (x(:,3) - mu(3))./ sigma(3);
end 
