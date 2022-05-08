clear;

matrix = zeros(20,4);

tau = 1/20;
t = 0:tau:1;
alpha = @(t)(2+sin(t))./4;
alpha_der = @(t)cos(t)./4;
g = @(t_n, sigma,alpha)sigma - (1 - 0.5 .* alpha(t_n+sigma.*tau));
g_der = @(t_n, sigma, alpha_der)1 + 0.5 .* tau .* alpha_der(t_n+sigma.*tau);

for i=1:length(t)
    t_n = t(i);
    matrix(i,:) = get_sigma_Newton(t_n, g, g_der, alpha, alpha_der);
end




function result = get_sigma_Newton(t_n, g, g_der, alpha, alpha_der)   
    result = zeros(4,1);
    sigma = 3./4;
    temp = sigma;
    k = 0;
    eps = 0.001;
    while k<4
        k = k + 1;
        sigma = temp-g(t_n, temp, alpha)./g_der(t_n, temp, alpha_der); 
        result(k) = sigma;
        if abs(temp-sigma) <= eps
           %break; 
        end
        temp = sigma;
    end
    %result = sigma;
end