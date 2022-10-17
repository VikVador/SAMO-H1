function plotObjF(functionID)
    %--------------
    % Documentation
    %--------------
    % This function has for purpose to plot a 3D surface to 
    % visualize the objective function of our statement
    %
    % Creation of the domain of vizualization
    [x, y] = meshgrid(-10 : 0.5 : 10, -10 : 0.5 : 10);

    % Computing the value of the objective function
    switch functionID
        case 1
            fvalue = 2 * x.^2 - 3 * x .* y + 2 * y.^2 - 2 * x + 10 * y - 1;
        case 2
            fvalue = x.^4 + x.^3 - 2 * x.^2 - 2 * x + y.^2;
    end

    % Plotting the 3D surface (Reducing color intensity and removing edges)
    plt = surf(x, y, fvalue, 'FaceAlpha', 0.8);
    plt.EdgeColor = 'none';
    xlabel('x [-]');
    ylabel('y [-]');
    zlabel("f_" + int2str(functionID) + "(x, y)");
    
    % Saving the result
    saveas(plt, "../graphs/surfaces/obj_fun_" + int2str(functionID) + ".pdf");
end