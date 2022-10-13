function fval = getObjFVal(X, functionID)
    %--------------
    % Documentation
    %--------------
    % This function returns the value of f_i, evaluated
    % at the point X = (x, y)
    %
    % Retreiving coordinates
    x = X(1);
    y = X(2);

    % Evaluating f_i at (x, y)
    switch functionID
        case 1
            fval = 2 * x^2 - 3 * x * y + 2 * y^2 - 2 * x + 10 * y - 1;
        case 2
            fval = x^4 + x^3 - 2 * x^2 - 2 * x + y^2;
    end
end