function newtheta = descend(theta, alpha, x, y, m)
  s = vec([0, 0]);
  for i = 1:m
    % Computing x(i, :)' because x(i) is represented horizontally
    s += (theta' * x(i, :)' - y(i)) * x(i, :)';
  end
  newtheta = theta - alpha / m * s;
end

function theta = gradient_descent(theta, alpha, x, y, m)
  for i = 1:1500
    theta = descend(theta, alpha, x, y, m);
  end
end

function height = estimate_height(theta, age)
  height = theta(1) + age * theta(2);
end

function val = J(theta, x, y, m)
  s = 0;
  for i = 1:m
    % Computing x(i, :)' because x(i) is represented horizontally
    s += (theta' * x(i, :)' - y(i)) ** 2
  end
  val = 1 / (2 * m) * s
end
  

function main_run(script = false)
  x = load('ex2x.dat');
  y = load('ex2y.dat');

  if (!script)
    figure;
    plot(x, y, 'o');
    ylabel('Height in meters');
    xlabel('Age in years');
  end
  %------------------------------------------------------------------
  % Problem 1 -------------------------------------------------------
  %------------------------------------------------------------------
  printf('Problem 1:\n');
  m = length(y);
  x = [ones(m, 1), x];
  alpha = 0.07;
  theta = vec([0, 0]);
  newtheta = descend(theta, alpha, x, y, m);
  printf('\tTheta = [%f, %f]\n', newtheta(1), newtheta(2));
  %------------------------------------------------------------------
  % Problem 2 -------------------------------------------------------
  %------------------------------------------------------------------
  printf('Problem 2:\n');
  theta = gradient_descent(theta, alpha, x, y, m);
  if (!script)
    hold on;
    plot(x(:,2), x * theta, '-');
    legend('Training data', 'Linear regression');
  end
  printf('\tTheta = [%f, %f]\n', theta(1), theta(2));
  %------------------------------------------------------------------
  % Problem 3 -------------------------------------------------------
  %------------------------------------------------------------------
  printf('Problem 3:\n');
  printf('\tEstimated height of 3.5 yo %f\n', estimate_height(theta, 3.5))
  printf('\tEstimated height of 7 yo %f\n', estimate_height(theta, 7))
  %------------------------------------------------------------------
  % Understanding J(theta)-------------------------------------------
  %------------------------------------------------------------------
  J_vals = zeros(100, 100);
  theta0_vals = linspace(-3, 3, 100);
  theta1_vals = linspace(-1, 1, 100);
  for i = 1:length(theta0_vals)
    for j = 1:length(theta1_vals)
      t = [theta0_vals(i), theta1_vals(j)];
      J_vals(i, j) = J(t', x, y, m);
    end
  end
  if (!script)
    J_vals = J_vals'
    figure;
    surf(theta0_vals, theta1_vals, J_vals);
    xlabel('\theta_0');
    ylabel('\theta_1');
  end
end
