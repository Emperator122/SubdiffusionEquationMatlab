function result = SweepMethod(matrix,f)
%Метод прогонки
    %matrix - матрица системы
    %f - столбец свободных членов
    n = length(f);
    alpha = zeros(n,1);
    beta = zeros(n,1);
    %Прямой ход
    alpha(2) = (-1 .* matrix(1, 2)) ./ matrix(1, 1);
    beta(2) = f(1) ./ matrix(1, 1);
    for i = 2:n-1
        alpha(i + 1) = (-1 .* matrix(i, i + 1)) ./ (matrix(i, i) + ...
                       matrix(i, i - 1) .* alpha(i));
        beta(i + 1) = (-1 .* matrix(i, i - 1) .* beta(i) + f(i)) ./ ...
                      (matrix(i, i) + matrix(i, i - 1) .* alpha(i));
    end
    %Обратный хол
    result = zeros(n,1);
    result(n) = (-1 .* matrix(n, n - 1) .* beta(n) + f(n)) ./ ...
                (matrix(n, n) + matrix(n, n - 1) .* alpha(n));
    for i = n-1:-1:1
       result(i) = alpha(i + 1) .* result(i + 1) + beta(i + 1); 
    end
end

