function [] = plotEvolutionObjF(x, f, functionID, parameters, save_plot)
    %--------------
    % Documentation
    %--------------
    % Plot the evolution of the objective function
    % throughout the optimization process
    %
    % Used to plot evolution
    iterations = linspace(1, size(x, 2), size(x, 2));

    % Plotting
    plt = plot(iterations, f(x(1, :), x(2, :)));
    xlabel('Iterations [-]');
    ylabel("f_" + int2str(functionID) + "(x, y)");
    
    % Name of the different method for the plot saving
    methodName = ["SDM"; "CGM"; "BFGS"];
    
    % Saving the plot
    if save_plot
        saveas(plt,"../graphs/results/evolution/evol_"   + methodName(parameters(1)) + ...
                   "_(" + int2str(parameters(2)) + "," + int2str(parameters(3))      + ")_"+ ...
                   "f_value_" + sprintf('%.6f', f(x(1, end), x(2, end)))             + "_" + ...
                   "nb_iter_"       + int2str(parameters(8))                         + "_" + ...
                   "nb_iter_total_" + int2str(parameters(4))                         + "_" + ...
                   "eps_"           + sprintf('%.6f', parameters(5))                 + "_" + ...
                   "nu_"            + sprintf('%.6f', parameters(6))                 + "_" + ...
                   "sc_"            +  int2str(parameters(7))                        + ".pdf");
    end
end