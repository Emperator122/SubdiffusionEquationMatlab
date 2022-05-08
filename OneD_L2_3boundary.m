clear;
addpath("Functions");
% точное решение
u_exact = @(hi,t)(t.^3+3.*t.^2+1).*(hi.^4+10.*hi+3); 

% Начальное условие
initialCondition = @(hi)hi.^4+10.*hi+3;

% для третьего краевого
beta_1 = @(t)3+t.^2;
beta_2 = @(t)t-1;
mu_1 = @(t)3.*t.^2.*(t.^3+3.*t.^2+1)-(t.^3+3.*t.^2+1);
mu_2 = @(t)14.*t.*(t.^3+3.*t.^2+1);


% alpha(t)
% alpha = @(t)(2+sin(t))./4;
alpha = @(t)1-t.^2;
% alpha = @(t)exp(-t);

% f(hi,t)
f = @(hi,t) ( 6.*t.^(3-alpha(t))./gamma(4-alpha(t)) + ...
   6.*t.^(2-alpha(t))./gamma(3-alpha(t))).*(hi.^4+10.*hi+3) - ...
   12.*hi.^2.*(t.^3+3.*t.^2+1);

% Границы для t
a_T = 0;
b_T = 1;
M = 320; % Колво разбиений отрезка [a_T; b_T]

% Границы для hi
a_x = 0;
b_x = 1;
N = 320; % количества разбиений каждого отрезка


% Расчеты
tau = (b_T-a_T)/M;
h = (b_x-a_x)/N;

t = a_T:tau:b_T;
hi = a_x:h:b_x;

t_layer = length(t);

u = zeros(N+1, M+1); % первая для xi, вторая для t

u(:, 1) = initialCondition(hi);
for n = 1:M
    step_coef = 1./h.^2;
    % зададим матрицу и столбец свободных членов для метода прогонки
    matrix = zeros(N+1, N+1);
    matrix_f = zeros(N+1, 1);
    
    % загоним в матрицу краевые условия
    mu_1_wave = 2.*mu_1(t(n+1))./h + f(hi(1), t(n+1));
    [known_part, unknown_coef] = kaputo_der(1, n, alpha, u, tau, t);
    matrix(1,1) = unknown_coef+2.*step_coef+2.*1./h.*beta_1(t(n+1));
    matrix(1,2) = -1.*2.*step_coef;
    matrix_f(1) = -1.*known_part+mu_1_wave;

    mu_2_wave = 2.*mu_2(t(n+1))./h + f(hi(N+1), t(n+1));
    [known_part, unknown_coef] = kaputo_der(N+1, n, alpha, u, tau, t);
    matrix(N+1,N) = -1.*2.*step_coef;
    matrix(N+1,N+1) = unknown_coef+2.*step_coef+2.*1./h.*beta_2(t(n+1));
    matrix_f(N+1) = -known_part+mu_2_wave;
    
    % загоним в матрицу остальные точки
    for i = 2:N
        [known_part, unknown_coef] = kaputo_der(i, n, alpha, u, tau, t);
        % для начала заполним саму матрицу
        matrix(i, i-1) = step_coef;
        matrix(i, i) = -1.*(unknown_coef + 2./h.^2);
        matrix(i, i+1) = step_coef;
        
        % а затем столбец свободных членов
        matrix_f(i) = known_part - f(hi(i), t(n+1));
    end
    % находим корни и записываем их на новый слой
    roots = matrix\matrix_f;
    u(:, n+1) = roots;
end

delta = abs(u_exact(hi', t)-u);
plot(hi, u(:, t_layer), hi, u_exact(hi', t(t_layer)));
legend("Численное решение", "Точное решение");
% err_E = E(u, u_exact(hi', t))

function [known_part, unknown_coef] = kaputo_der(i, j, alpha, u, tau, t)
    gamma_coef = tau.^(-alpha(t(j+1))) ./ gamma(2-alpha(t(j+1)));
    known_part = -1.*a(1,j, alpha, t).*u(i, j);
    for k = 2:j
        known_part = known_part + a(k, j, alpha, t).*(u(i, j-k+2) - u(i, j-k+1));
    end
    known_part = known_part .* gamma_coef;
    
    unknown_coef = a(1, j, alpha, t).*gamma_coef;
end

% "a" переписана под индексацию с единицы
function result = a(k, n, alpha, t)
    if k == 1
       result = 1; 
    else
        k = k-1;
        result = (k+1).^(1-alpha(t(n+1)))-k.^(1-alpha(t(n+1)));
    end
end