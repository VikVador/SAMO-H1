function method = terminal(index, parameters, ls_method)
    %--------------
    % Documentation
    %--------------
    % A simple function to display information over the 
    % terminal during the whole optimization process
    %
    switch index

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
            sc_name = ["max(grad_f) < eps";
                       "norm(grad_f, 2) < eps";
                       "f(x_(k+1)) - f(x_k) < nu"];
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
            disp("Method         : " + functionsName(parameters(1)));
            disp(" ");
            switch ls_method
                case 'NR'
                    disp("Alpha comput.  : Newton Raphson");
                case 'S'
                    disp("Alpha comput.  : Secant");
                case 'D'
                    disp("Alpha comput.  : Dichotomy");
                case 'BB'
                    disp("Alpha comput.  : Black Box");
                case 'DIV'
                    disp("Alpha comput.  : Divergent Serie");
                case 'CQ'
                    disp("Alpha comput.  : Convex quadratic function");
            end
            disp(" ");
            disp("Stop. Crit.    : " + sc_name(parameters(7)));
            disp(" ");
            disp("Max iterations : " + int2str(parameters(4)));
            disp(" ");
            disp("X_0            : (" + int2str(parameters(2)) + ", " + int2str(parameters(3)) + ")");
            disp(" ");
            disp("Epsilon        : " + sprintf('%.10f', parameters(5)));
            disp(" ");
            disp("Nu             : " + sprintf('%.10f', parameters(6)));
            disp(" ");
            
         % Plotting the results (1)
        case 3
            disp("-------------------");
            disp("      Plotting     ");
            disp("-------------------");
            disp(" ");
            disp("Plotting objective functions : ...");
            disp(" ");

        % Plotting the results & optimizing (2)
        case 4
            clc;
            terminal(2, parameters, ls_method);
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
            disp("Solution                        : ...");
        
        % Showing the results
        case 5
            clc;
            terminal(2, parameters, ls_method);
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
            if parameters(4) == parameters(8)
                disp("Solution                         : Maximum number of iterations reached");
            else
                disp("Solution                         : Done");
            end
            disp(" ");
    end
end










