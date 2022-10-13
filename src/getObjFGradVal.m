function grad_f = getObjFGradVal(X, functionID)
    %--------------
    % Documentation
    %--------------
    % This function returns the values of nabla(f_i), evaluated
    % at the point X = (x, y)
    %
    % Retreiving coordinates
    x = X(1);
    y = X(2);

    % Stores the value of the gradient
    grad_f = [0, 0];

    % Evaluating f_i at (x, y)
    switch functionID
        case 1
            grad_f(1) =   4 * x - 3 * y - 3;  % dx
            grad_f(2) = - 3 * x + 4 * y + 10; % dy
        case 2
            grad_f(1) = 4 * x^3 + 3 * x^2 - 4 * x - 2; % dx
            grad_f(2) = 2 * y;                         % dy
    end
end