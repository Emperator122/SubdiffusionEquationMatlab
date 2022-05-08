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
alpha = @(t)(2+sin(t))./4;

% f(hi,t)
f = @(hi,t) ( 6.*t.^(3-alpha(t))./gamma(4-alpha(t)) + ...
   6.*t.^(2-alpha(t))./gamma(3-alpha(t)) + t.^3+3.*t^2+1 ).*sin(hi);

% Границы для t
a_T = 0;
b_T = 1;
M = 640; % Колво разбиений отрезка [a_T; b_T]

% Границы для hi
a_x = 0;
b_x = pi;
N = round(640*pi); % количества разбиений каждого отрезка


% Расчеты
tau = (b_T-a_T)/M;
h = (b_x-a_x)/N;

t = a_T:tau:b_T;
hi = a_x:h:b_x;

t_layer = length(t);

u = zeros(N+1, M+1); % первая для xi, вторая для t
tic
u(:, 1) = initialCondition(hi);
for n = 1:M
    step_coef = 1./h.^2;
    a_values = calculate_a(n, alpha, t);
    gamma_value = tau.^(-alpha(t(n+1))) ./ gamma(2-alpha(t(n+1)));
    unknown_coef = kaputo_der_unknown_coef(a_values, gamma_value);
    % зададим матрицу и столбец свободных членов для метода прогонки
    matrix = zeros(N+1, N+1);
    matrix_f = zeros(N+1, 1);
    
    % загоним в матрицу краевые условия
    matrix(1,1) = 1;
    matrix(N+1,N+1) = 1;
    matrix_f(1) = leftBound(t(n+1));
    matrix_f(N+1) = rightBound(t(n+1));
    
    % загоним в матрицу остальные точки
    for i = 2:N
        known_part = kaputo_der_known_part(a_values, i, n, u, gamma_value);
        % для начала заполним саму матрицу
        matrix(i, i-1) = step_coef;
        matrix(i, i) = -1.*(unknown_coef + 2./h.^2);
        matrix(i, i+1) = step_coef;
        
        matrix_f(i) = known_part - f(hi(i), t(n+1));
    end
    % находим корни и записываем их на новый слой
    %     roots = SweepMethod(matrix, matrix_f);
    roots = tridiagonal(matrix, matrix_f);
%     roots = matrix\matrix_f;
    u(:, n+1) = roots;
end
toc
delta = abs(u_exact(hi', t)-u);
plot(hi, u(:, t_layer), hi, u_exact(hi', t(t_layer)));
legend("Численное решение", "Точное решение");


function known_part = kaputo_der_known_part(a_values, i, j, u, gamma_value)
        known_part = (sum(a_values(2:j) .* (u(i, j:-1:2)-u(i, j-1:-1:1)))-a_values(1).*u(i, j)) ...
        .* gamma_value;    
end

function unknown_coef = kaputo_der_unknown_coef(a_values, gamma_value)
    unknown_coef = a_values(1).* gamma_value;
end