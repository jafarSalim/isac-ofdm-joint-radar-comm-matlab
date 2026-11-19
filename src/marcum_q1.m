function q = marcum_q1(a, b)
%MARCUM_Q1  First-order Marcum Q function by its Bessel series (no toolbox).
q = zeros(size(a));
for i = 1:numel(a)
    if a(i) == 0, q(i) = exp(-b^2/2); continue; end
    x = a(i) * b;  s = 0;
    if a(i) < b
        for k = 0:400
            t = (a(i)/b)^k * besseli(k, x, 1) * exp(x - (a(i)^2 + b^2)/2);
            s = s + t;
            if k > 5 && t < 1e-16, break; end
        end
        q(i) = s;
    else   % Q1(a,b) = 1 - sum_{k>=1} (b/a)^k I_k(ab) e^{-(a^2+b^2)/2}
        for k = 1:400
            t = (b/a(i))^k * besseli(k, x, 1) * exp(x - (a(i)^2 + b^2)/2);
            s = s + t;
            if k > 5 && t < 1e-16, break; end
        end
        q(i) = 1 - s;
    end
end
end
