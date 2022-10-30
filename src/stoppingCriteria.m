function stop = stoppingCriteria(i, s, eps, f, nu, x_k, x_kplus1)
    %--------------
    % Documentation
    %--------------
    % This function determines if the stopping criteria is
    % reached or not where, for a given value of i, one has:
    %
    % 1 - max_value(abs(grad_f(x_k+1))) < eps
    %
    % 2 - norm(grad_f)^2 < eps
    %
    % 3 - |f(x_k+1) - f(x_k)| < nu
    %
    stop = false;

    switch i
        case 1
            if max(abs(s)) < eps
                stop = true;
            end

        case 2
            if norm(s, 2) < eps
                stop = true;
            end

        case 3
            res = abs(f(x_kplus1(1), x_kplus1(2)) - f(x_k(1), x_k(2)) );
            if res < nu
                stop = true;
            end
    end
end
