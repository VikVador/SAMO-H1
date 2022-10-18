function plotObjF(f)
    %--------------
    % Documentation
    %--------------
    % This function has for purpose to plot a 3D surface to 
    % visualize the objective function of our statement
    %
    % Creation of the domain of vizualization
    [x, y] = meshgrid(-10 : 0.5 : 10, -10 : 0.5 : 10);

    % Computing the value of the objective function
    
    fvalue = f(x, y);

    % Plotting the 3D surface (Reducing color intensity and removing edges)
    plt = surf(x, y, fvalue, 'FaceAlpha', 0.8);
    plt.EdgeColor = 'none';
    xlabel('x [-]');
    ylabel('y [-]');
    zlabel("f_" + int2str(functionID) + "(x, y)");
    
    % Saving the result
    saveas(plt, "../graphs/surfaces/obj_fun_" + int2str(functionID) + ".pdf");
end