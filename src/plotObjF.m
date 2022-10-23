function plotObjF(f, functionID)
    %--------------
    % Documentation
    %--------------
    % Plot a 3D surface to visualize the objective function
    %
    % Plotting the 3D surface (Reducing color intensity and removing edges)
    plt = fsurf(f, [-10 10 -10 10], 'FaceAlpha', 0.5);
    plt.EdgeColor = 'none';
    xlabel('x [-]');
    ylabel('y [-]');
    zlabel("f_" + int2str(functionID) + "(x, y)");
    
    % Saving the result
    saveas(plt, "../graphs/surfaces/obj_fun_" + int2str(functionID) + ".pdf");
end