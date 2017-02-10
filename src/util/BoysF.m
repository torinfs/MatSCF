function y = BoysF(m, T)
    %gamma(x) = integral from 0 to inf of t^(x-1) exp(-t) dt
    
    %gammainc(x,a) = 1 ./ gamma(a) .* integral from 0 to x of t^(a-1) exp(-t) dt
    %Therefore, lower gamma should be gammainc(x,a).*gamma(a)
    
    %gammainc(x, a, 'upper') = 1 - gammainc(x, a)
    %Therefore, upper gamma should be gammainc(x, a, 'upper').*gamma(a)
    
    %We also need to set some threshold to switch from the Boys function to
    %the limiting behavior
    for i = 1:numel(m)
        if T(i) < 1e-8
            y = 1/(2*m(i) + 1);
        else
            y = (gammainc(T(i), m(i)+0.5)*gamma(m(i)+0.5))/(2*T(i)^(m(i)+0.5));
        end
    end
end

