clear;
addpath("Functions");
% точное решение
u_exact = @(hi,t)(t.^3+3.*t.^2+1).*sin(hi); 

% Начальное условие
initialCondition = @(hi)sin(hi);

% Граничные условия
leftBound = @(t)0;
rightBound = @(t)0;

% alpha(t)
% alpha = @(t)(2+sin(t))./4;
% alpha_der = @(t)cos(t)./4;
% alpha = @(t)1-t.^2;
% alpha_der = @(t)-2.*t;
alpha = @(t)exp(-t);
alpha_der = @(t)-exp(-t);

% f(hi,t)
f = @(hi,t) ( 6.*t.^(3-alpha(t))./gamma(4-alpha(t)) + ...
   6.*t.^(2-alpha(t))./gamma(3-alpha(t)) + t.^3+3.*t^2+1 ).*sin(hi);

% Границы для t
a_T = 0;
b_T = 1;
M = 6400; % Колво разбиений отрезка [a_T; b_T]

% Границы для hi
a_x = 0;
b_x = pi;
N = round(80*pi); % количества разбиений каждого отрезка


% Расчеты
tau = (b_T-a_T)/M;
h = (b_x-a_x)/N;

t = a_T:tau:b_T;
hi = a_x:h:b_x;

u = zeros(N+1, M+1); % первая для xi, вторая для t
tic
u(:, 1) = initialCondition(hi);
for n = 1:M
    % вычислим sigma для данной итерации
    sigma = get_sigma_Newton(alpha, alpha_der, t, n, tau, 0.001);
%     sigma = get_sigma(alpha, t, n, tau, 0.001);
    % зададим матрицу и столбец свободных членов для метода прогонки
    matrix = zeros(N+1, N+1);
    matrix_f = zeros(N+1, 1);

    % необходимые значения
    step_coef = sigma./h.^2;
    sigmaStep = t(n)+sigma.*tau;
    gamma_value = tau.^(-alpha(sigmaStep))./gamma(2-alpha(sigmaStep));
    
    % вычислим "c"
    c_values = calculate_c(n, alpha, sigma, sigmaStep);

    unknown_coef = kaputo_der_unknown_coef(c_values, gamma_value);

    % загоним в матрицу краевые условия
    matrix(1,1) = 1;
    matrix(N+1,N+1) = 1;
    matrix_f(1) = leftBound(t(n+1));
    matrix_f(N+1) = rightBound(t(n+1));

    % загоним в матрицу остальные точки
    for i = 2:N
        known_part = kaputo_der_known_part(c_values, i, n, u, gamma_value);
        % для начала заполним саму матрицу
        matrix(i, i-1) = step_coef;
        matrix(i, i) = -1.*(unknown_coef + 2.*step_coef);
        matrix(i, i+1) = step_coef;
        
        % а затем столбец свободных членов
        matrix_f(i) = known_part - f(hi(i), sigmaStep) ...
            - (1-sigma).*(u(i+1, n)-2.*u(i, n)+u(i-1, n))./h.^2;
    end
    % находим корни и записываем их на новый слой
    roots = SweepMethod(matrix, matrix_f);
%     roots = tridiagonal(matrix, matrix_f);
%     roots = matrix\matrix_f;
    u(:, n+1) = roots;
end
toc

delta = abs(u_exact(hi', t)-u);
plot(hi, u(:, length(t)), hi, u_exact(hi', t(length(t))));
legend("Численное решение", "Точное решение");

function known_part = kaputo_der_known_part(c_values, i, j, u, gamma_value)
        known_part = (sum(c_values(2:j) .* (u(i, j:-1:2)-u(i, j-1:-1:1)))-c_values(1).*u(i, j)) ...
        .* gamma_value;    
end

function unknown_coef = kaputo_der_unknown_coef(c_values, gamma_value)
    unknown_coef = c_values(1).* gamma_value;
end