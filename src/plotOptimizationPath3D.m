function [] = plotOptimizationPath3D(X, functionID)
    %--------------
    % Documentation
    %--------------
    % This function has for purpose to plot a 3D surface to 
    % visualize the objective function of our statement as well
    % as the optimization path of our OM
    %
    %--------
    % Pathway
    %--------
    % Contains the number of iterations of our algorithm
    nb_iter = size(X, 2);

    % Stores the solution
    solution = zeros(1, nb_iter);

    % Computing solution value at each step
    for i  = 1 : nb_iter
        solution(i) = getObjFVal(X(:, i), functionID);
    end

    %--------
    % Surface
    %--------
    % Creation of the domain of vizualization for the surface
    [xs, ys] = meshgrid(-100:0.5:100, -100:0.5:100);

    % Computing the value of the objective function
    switch functionID
        case 1
            fvalue = 2 * xs.^2 - 3 * xs .* ys + 2 * ys.^2 - 2 * xs + 10 * ys - 1;
        case 2
            fvalue = xs.^4 + xs.^3 - 2 * xs.^2 - 2 * xs + ys.^2;
    end

    %-------------------------------------------
    % Plotting the surface and optimization path
    %-------------------------------------------
    % Setting the plot to dim = 3 for a 3-D plot
    view(3);
    hold on;

    % Pathway
    plot3(X(1, :), X(2, :), solution, 'Color','r','LineWidth', 3)

    % Storing initial and ending value
    init_point  = [X(1, 1), X(2, 1), solution(1)];
    final_point = [X(1, end), X(2, end), solution(end)];

    % Creating point as text to plot it
    init_text = "(" + sprintf('%.2f', init_point(1)) + "," + ...
                      sprintf('%.2f', init_point(2)) + "," + ...
                      sprintf('%.2f', init_point(3)) + ")";

    final_text = "(" + sprintf('%.2f', final_point(1)) + "," + ...
                       sprintf('%.2f', final_point(2)) + "," + ...
                       sprintf('%.2f', final_point(3)) + ")";

    % Plotting points
    scatter3(init_point(1) , init_point(2)  , init_point(3), 'filled', "SizeData", 100)
    scatter3(final_point(1), final_point(2), final_point(3),'*',"SizeData", 100)

    % Adding labels
    text(init_point(1), init_point(2), init_point(3) + 0.5, init_text)
    text(final_point(1), final_point(2), final_point(3) + 0.5,final_text)
    
    % Surface
    plt = surf(xs, ys, fvalue, 'FaceAlpha', 0.5);
    plt.EdgeColor = 'none';
    
end