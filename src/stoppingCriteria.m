function stop = stoppingCriteria(i, x_k, x_k_1, functionID, eps, nu)
    %--------------
    % Documentation
    %--------------
    % This function determines if the stopping criteria is
    % reached or not where, for a given value of i, one has:
    %
    % 1 - max_value(grad_f(x_k+1)) < eps
    % 
    % 2 - sum(grad_f^2) < eps
    %
    % 3 - |f(x_k+1) - f(x_k)| < nu
    %
    % Initialization of the SC value
    stop = false;

    switch i
        case 1
            % Step 1 - Computing grad
            grad = getObjFGradVal(x_k_1, functionID);

            % Step 2 - Check if max value is < eps
            if max(grad) < eps
                stop = true;
            end

        case 2
            % Step 1 - Computing grad
            grad = getObjFGradVal(x_k_1, functionID);

            % Step 2 - Check if sum of all grad element to ^2 is < eps
            if sum(grad.^2) < eps
                stop = true;
            end
        
        case 3
            % Step 1 - Computing residual
            res = abs(getObjFVal(x_k_1, functionID) - getObjFVal(x_k, functionID));

            % Step 2 - Check if residual is small enough
            if res < nu
                stop = true;
            end
    end
end


