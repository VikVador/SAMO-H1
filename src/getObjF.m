function f = getObjF(X, functionID)
    %--------------
    % Documentation
    %--------------
    % This function is used to define the symbolic objective function
    %
    x = X(1); y = X(2);

    switch functionID
        case 1
            f = 2 * x^2 - 3 * x * y + 2 * y^2 - 2 * x + 10 * y - 1;
        case 2
            f = x^4 + x^3 - 2 * x^2 - 2 * x + y^2;
        case 3
            f = y * sin(x) - x * cos(y);
    end
end