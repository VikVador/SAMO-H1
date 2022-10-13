function method = terminal(i, parameters)
    %--------------
    % Documentation
    %--------------
    % A simple function to display information over the terminal
    %
    switch i

        % Initialization
        case 1
            disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
            disp("%                                                                         %");
            disp("%               Structural and Multidisciplinary Optimization             %");
            disp("%                                                                         %");
            disp("%                      H1 - Unconstrained Optimization                    %");
            disp("%                                                                         %");
            disp("% @ Arnaud Rémi                                        @ Victor Mangeleer %");
            disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
            disp(" ");
            disp('Which optimization method do you want to use ?')
            disp(" ");
            disp('     1 - Steepest descent method')
            disp(" ");
            disp('     2 - Conjugate gradients method with Fletcher-Reeves update rule')
            disp(" ");
            disp('     3 - BFGS Quasi-Newton method')
            disp(" ");
            prompt = 'Choice:';
            method = input(prompt);
            
        % Parameters of the simulation
        case 2
            clc;
            functionsName = ["Steepest descent method";
                             "Conjugate gradients method with Fletcher-Reeves update rule";
                             "BFGS Quasi-Newton method"];
            disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
            disp("%                                                                         %");
            disp("%               Structural and Multidisciplinary Optimization             %");
            disp("%                                                                         %");
            disp("%                      H1 - Unconstrained Optimization                    %");
            disp("%                                                                         %");
            disp("% @ Arnaud Rémi                                        @ Victor Mangeleer %");
            disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
            disp(" ");
            disp("-------------------");
            disp("     Parameters    ");
            disp("-------------------");
            disp(" ");
            disp("Method  : " + functionsName(parameters(1)));
            disp(" ");
            disp("Start   : (" + int2str(parameters(2)) + "," + int2str(parameters(3)) + ")");
            disp(" ");
            disp("Nb iter : " + int2str(parameters(4)));
            disp(" ");
            disp("Epsilon : " + sprintf('%.6f', parameters(5)));
            disp(" ");
            disp("Nu      : " + sprintf('%.6f', parameters(6)));
            disp(" ");
            disp("SC Crit : " + int2str(parameters(7)));
            disp(" ");

         % Plotting the results
        case 3
            disp("-------------------");
            disp("      Plotting     ");
            disp("-------------------");
            disp(" ");
            disp("Plotting objective functions : ...");
            disp(" ");

        case 4
            clc;
            terminal(2, parameters);
            disp("-------------------");
            disp("      Plotting     ");
            disp("-------------------");
            disp(" ");
            disp("Plotting objective functions : done");
            disp(" ");
            close all;
            disp("-------------------");
            disp("    Optimizing     ");
            disp("-------------------");
            disp(" ");
            disp("Solution                : ...");
        
        case 5
            clc;
            terminal(2, parameters);
            disp("-------------------");
            disp("      Plotting     ");
            disp("-------------------");
            disp(" ");
            disp("Plotting objective functions : done");
            disp(" ");
            close all;
            disp("-------------------");
            disp("    Optimizing     ");
            disp("-------------------");
            disp(" ");
            disp("Solution                 : found");
            disp(" ");
    end
end










