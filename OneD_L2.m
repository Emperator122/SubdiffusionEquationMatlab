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
M = 100; % Колво разбиений отрезка [a_T; b_T]

% Границы для hi
a_x = 0;
b_x = pi;
N = round(100*pi); % количества разбиений каждого отрезка


% Расчеты
tau = (b_T-a_T)/M;
h = (b_x-a_x)/N;

t = a_T:tau:b_T;
hi = a_x:h:b_x;

u = zeros(N+1, M+1); % первая для xi, вторая для t

u(:, 1) = initialCondition(hi);
for n = 1:M
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
        step_coef = 1./h.^2;
        gamma_coef = tau.^(-alpha(t(n+1))) ./ gamma(2-alpha(t(n+1)));
        % для начала заполним саму матрицу
        matrix(i, i-1) = step_coef;
        matrix(i, i) = -1.*(a(1, n, alpha, t).*gamma_coef + 2./h.^2);
        matrix(i, i+1) = step_coef;
        
        % а затем столбец свободных членов
        % для заполнения столбца с.ч. необходимо вычислять сумму
        % вычислим ее
        sum = -1.*a(1,n, alpha, t).*u(i, n);
        for k = 2:n
            sum = sum + a(k, n, alpha, t).*(u(i, n-k+2) - u(i, n-k+1));
        end
        sum = sum .* gamma_coef;
        sum = sum - f(hi(i), t(n+1));
        
        % результат вычислений закидываем в столбец с.ч.
        matrix_f(i) = sum;
    end
    % находим корни и записываем их на новый слой
    roots = matrix\matrix_f;
    u(:, n+1) = roots;
end

delta = abs(u_exact(hi', t)-u);
% err_E = E(u, u_exact(hi', t))



% "a" переписана под индексацию с единицы
function result = a(k, n, alpha, t)
    if k == 1
       result = 1; 
    else
        k = k-1;
        result = (k+1).^(1-alpha(t(n+1)))-k.^(1-alpha(t(n+1)));
    end
end