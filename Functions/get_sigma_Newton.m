function result = get_sigma_Newton(alpha, alpha_der, t, n, tau, eps)
    g = @(t_n, sigma,alpha)sigma - (1 - 0.5 .* alpha(t_n+sigma.*tau));
    g_der = @(t_n, sigma, alpha_der)1 + 0.5 .* tau .* alpha_der(t_n+sigma.*tau);
    sigma = 3./4;
    temp = sigma;
    while true
        sigma = temp-g(t(n), temp, alpha)./g_der(t(n), temp, alpha_der); 
        if abs(temp-sigma) <= eps
           break; 
        end
        temp = sigma;
    end
    result = sigma;
end

