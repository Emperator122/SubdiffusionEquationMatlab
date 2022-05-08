function result = calculate_a(n, alpha, t)
    alpha_value = alpha(t(n+1));
    result = zeros(1, n);
    for k = 1:n
        result(k) = a(k, alpha_value);
    end
end

% "a" переписана под индексацию с единицы
function result = a(k, alpha_value)
    if k == 1
       result = 1; 
    else
        k = k-1;
        result = (k+1).^(1-alpha_value)-k.^(1-alpha_value);
    end
end