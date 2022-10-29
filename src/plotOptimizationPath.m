function [] = plotOptimizationPath(x, fun, functionID, parameters, save_plot)
    %--------------
    % Documentation
    %--------------
    % This function has for purpose to plot iso-curves to 
    % visualize the objective function of our statement as well
    % as the optimization path of our optimization method
    %
    % Definition of the domain size
    lb = -10;
    up = 10;

    xi = lb : 0.25 : up;
    f = zeros(length(xi), length(xi));

    % Computing objective function's v
    for i = 1:length(xi)
        f(:, i) = fun(xi(i), xi(:));
    end

    % Initialization of figure
    figure('name','Optimization path')
    hold on
    xlabel('x')
    ylabel('y')
    axis([lb up lb up])
    title('Optimization path')

    % Plotting iso-curve
    switch functionID
        case 1
            [C,h]=contour(xi,xi,f,[0:1:2 2:6:20 20:20:80 80:40:200 200:100:2000]);
            clabel(C,h);
        case 2
            [C,h]=contour(xi,xi,f,([-10:2:2 2:4:100]));
            clabel(C,h);
    end
    hold on
    
    % Adding iso-curve value
    ind=1;
    for i=1:size(x,2)-1
        plot(x(1,i),x(2,i),'.c','markersize',30)
        plot([x(1,i) x(1,i+1)],[x(2,i) x(2,i+1)],'c','linewidth',2)
        text(x(1,i),x(2,i),num2str(ind-1),'horizontalalignment','center','verticalalignment','middle')
        ind=ind+1;
    end
    plt = plot(x(1,end),x(2,end),'.c','markersize',30);
    text(x(1,end),x(2,end),num2str(ind-1),'horizontalalignment','center','verticalalignment','middle')
    
    % Contains all the methods name to save the plot
    methodName = ["SDM"; "CGM"; "BFGS"];
    
    % Saving the plot
    if save_plot
        saveas(plt,"../graphs/results/2D/opti_path_"   + methodName(parameters(1))   + ...
                   "_(" + int2str(parameters(2)) + "," + int2str(parameters(3))      + ")_"+ ...
                   "f_value_" + sprintf('%.6f', fun(x(1, end), x(2, end)))    + "_" + ...
                   "nb_iter_"       + int2str(parameters(8))                         + "_" + ...
                   "nb_iter_total_" + int2str(parameters(4))                         + "_" + ...
                   "eps_"           + sprintf('%.6f', parameters(5))                 + "_" + ...
                   "nu_"            + sprintf('%.6f', parameters(6))                 + "_" + ...
                   "sc_"            +  int2str(parameters(7))                        + ".pdf");
    end
end