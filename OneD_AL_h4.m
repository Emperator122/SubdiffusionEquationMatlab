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
M = 400; % Колво разбиений отрезка [a_T; b_T]

% Границы для hi
a_x = 0;
b_x = pi;
N = round(20*pi); % количества разбиений каждого отрезка


% Расчеты
tau = (b_T-a_T)/M;
h = (b_x-a_x)/N;

t = a_T:tau:b_T;
hi = a_x:h:b_x;

u = zeros(N+1, M+1); % первая для xi, вторая для t

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

    % вычислим "c"
    c_values = calculate_c(n, alpha, sigma, sigmaStep);

    % загоним в матрицу краевые условия
    matrix(1,1) = 1;
    matrix(N+1,N+1) = 1;
    matrix_f(1) = leftBound(t(n+1));
    matrix_f(N+1) = rightBound(t(n+1));
    
    % загоним в матрицу остальные точки
    for i = 2:N
        [known_part, unknown_coef] = kaputo_der(c_values, i, n, alpha, u, tau, sigmaStep);
        % для начала заполним саму матрицу
        matrix(i, i-1) = step_coef-unknown_coef./12;
        matrix(i, i) = -1.*(10.*unknown_coef./12 + 2.*step_coef);
        matrix(i, i+1) = step_coef-unknown_coef./12;
        
        % а затем столбец свободных членов
        sum = known_part;
        sum = sum - A_func(i, hi, f, sigmaStep);
        sum = sum - (1-sigma).*(u(i+1, n)-2.*u(i, n)+u(i-1, n))./h.^2;
        
        % результат вычислений закидываем в столбец с.ч.
        matrix_f(i) = sum;
    end
    % находим корни и записываем их на новый слой
    roots = matrix\matrix_f;
    u(:, n+1) = roots;
end

delta = abs(u_exact(hi', t)-u);
%err_E = E(u, u_exact(hi', t))
plot(hi, u(:, length(t)), hi, u_exact(hi', t(length(t))));
legend("Численное решение", "Точное решение");


% Оператор A на известном слое сетки
function result = A_mesh(i, n, u)
    if i == 1 || i == size(u, 1)
        result = u(i);
        return;
    end
    result = (u(i+1,n)+10.*u(i,n)+u(i-1,n))./12;
end

% Оператор A для произвольной функции
function result = A_func(i, x, f, t_point)
    
    if i == 1 || i == length(x)
        if(t_point == Inf)
           result = f(x(i)); 
        else
           result= f(x(i), t_point);
        end
        return;
    end

    if(t_point == Inf)
        result = (f(x(i+1))+10.*f(x(i))+f(x(i-1)))./12;
    else
        result = (f(x(i+1), t_point)+10.*f(x(i), t_point)+f(x(i-1), t_point))./12;
    end
end

function [known_part, unknown_coef] = kaputo_der(c_values, i, j, alpha, u, tau, sigmaStep)
    known_part = -c_values(1).* A_mesh(i, j, u);
        for k = 2:j
            known_part = known_part + c_values(k).*...
                (A_mesh(i, j-k+2, u) - A_mesh(i, j-k+1, u));
        end
    known_part = known_part .* tau.^(-alpha(sigmaStep)) ./ gamma(2-alpha(sigmaStep));
    
    unknown_coef = c_values(1).* tau.^(-alpha(sigmaStep)) ./ gamma(2-alpha(sigmaStep));
end