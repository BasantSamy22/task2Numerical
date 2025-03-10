% Numerical Solution of ODE using Euler Methods

% Parameters
t0 = 0;      % Initial time
tf = 2;      % Final time
y0 = 1;      % Initial condition
h = 0.04;    % Step size
t = t0:h:tf; % Time vector
N = length(t);

% Initialize solutions
y_forward = zeros(1, N);
y_modified = zeros(1, N);
y_backward = zeros(1, N);
y_rk4 = zeros(1, N);
y_ab2 = zeros(1, N);

y_forward(1) = y0;
y_modified(1) = y0;
y_backward(1) = y0;
y_rk4(1) = y0;
y_ab2(1) = y0;

% Forward Euler Method
for i = 1:N-1
    y_forward(i+1) = y_forward(i) + h * (-50 * y_forward(i) + sin(t(i)));
end

% Modified Euler Method (Heun's Method)
for i = 1:N-1
    y_pred = y_modified(i) + h * (-50 * y_modified(i) + sin(t(i))); % Predictor Step
    y_modified(i+1) = y_modified(i) + (h/2) * ((-50 * y_modified(i) + sin(t(i))) + (-50 * y_pred + sin(t(i+1)))); % Corrector Step
end

% Backward Euler Method
for i = 1:N-1
    y_backward(i+1) = (y_backward(i) + h * sin(t(i+1))) / (1 + 50 * h);
end

% Runge-Kutta 4th Order (RK4)
for i = 1:N-1
    y_rk4(i+1) = y_rk4(i) + h * (-50 * y_rk4(i) + sin(t(i))) + (h/2) * (-50 * (y_rk4(i) + h * (-50 * y_rk4(i) + sin(t(i)))) + sin(t(i) + h));
end

% Adams-Bashforth 2-step (AB2)
% First step uses Euler method
y_ab2(2) = y_ab2(1) + h * (-50 * y_ab2(1) + sin(t(1)));
for i = 2:N-1
    y_ab2(i+1) = y_ab2(i) + (h/2) * (3 * (-50 * y_ab2(i) + sin(t(i))) - (-50 * y_ab2(i-1) + sin(t(i-1))));
end

% Plot Results
figure;
plot(t, y_forward, 'r', 'LineWidth', 1.5); hold on;
plot(t, y_modified, 'b', 'LineWidth', 1.5);
plot(t, y_backward, 'g', 'LineWidth', 1.5);
plot(t, y_rk4, 'm', 'LineWidth', 1.5);
plot(t, y_ab2, 'c', 'LineWidth', 1.5);
legend('Forward Euler', 'Modified Euler', 'Backward Euler' ,'RK4', 'Adams-Bashforth 2');
xlabel('Time t'); ylabel('y(t)');
title('Numerical Solutions using Euler Methods');
grid on;
