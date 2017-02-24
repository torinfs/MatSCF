function y = boysF(T, m)
    %gamma(x) = integral from 0 to inf of t^(x-1) exp(-t) dt
    
    %gammainc(x,a) = 1 ./ gamma(a) .* integral from 0 to x of t^(a-1) exp(-t) dt
    %Therefore, lower gamma should be gammainc(x,a).*gamma(a)
    
    %gammainc(x, a, 'upper') = 1 - gammainc(x, a)
    %Therefore, upper gamma should be gammainc(x, a, 'upper').*gamma(a)
    
    %We also need to set some threshold to switch from the Boys function to
    %the limiting behavior
    if T < 1e-10
        y = 1./(2.*m+1);
    else
        mp = m + 0.5;
        y = gammainc(T, mp, 'lower').*gamma(mp)./(2*T.^mp);
    end
end

