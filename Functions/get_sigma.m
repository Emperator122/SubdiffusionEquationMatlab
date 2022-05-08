function result = get_sigma(alpha, t, n, tau, eps)
    sigma = 3/4;
    temp = sigma;
    while true
        sigma = 1-alpha(t(n)+temp.*tau)./2; 
        if abs(temp-sigma) <= eps
           break; 
        end
        temp = sigma;
    end
    result = sigma;
end

