function fval = getObjFVal(x,functionID)
    switch functionID
        case 1
            fval = 2*x(1)^2 - 3*x(1)*x(2) + 2*x(2)^2 - 2*x(1) + 10*x(2) - 1;
        case 2
            fval = x(1)^4 + x(1)^3 - 2*x(1)^2 - 2*x(1) + x(2)^2;
    end
end