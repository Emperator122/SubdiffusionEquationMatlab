function result = calculate_c(n, alpha, sigma, sigmaStep)
    alpha_value = alpha(sigmaStep);
    result = zeros(1, n);
    for k = 1:n
        result(k) = c(k, n, sigma, alpha_value);
    end
end

% "c" переписана под индексацию с единицы
function result = c(k, n, sigma, alpha_value)
    if n == 1
        result = sigma.^(1-alpha_value);
    else
        if k == 1
           result = ((1+sigma).^(2-alpha_value) - ...
               sigma.^(2-alpha_value))./(2-alpha_value) - ...
               ((1+sigma).^(1-alpha_value) - ...
               sigma.^(1-alpha_value))./2;
        elseif k == n
            k = k-1;
            result = -1.*((k+sigma).^(2-alpha_value) - ...
                (k+sigma-1).^(2-alpha_value))./(2-alpha_value) + ...
                (3.*(k+sigma).^(1-alpha_value) - ...
                (k+sigma-1).^(1-alpha_value))./2;
        else
            k = k-1;
            result = ((k+sigma+1).^(2-alpha_value) - ...
                2.*(k+sigma).^(2-alpha_value) + ...
                (k+sigma-1).^(2-alpha_value))./(2-alpha_value) - ...
                ((k+sigma+1).^(1-alpha_value) - 2.* ...
                (k+sigma).^(1-alpha_value) + ...
                (k+sigma-1).^(1-alpha_value))./2;
        end
    end
end