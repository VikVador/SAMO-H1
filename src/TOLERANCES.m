%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               Structural and Multidisciplinary Optimization             %
%                                                                         %
%                      H1 - Unconstrained Optimization                    %
%                                                                         %
% @ Arnaud Rémi                                        @ Victor Mangeleer %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc
%% ----------
%  Parameters
%  ----------
% Defines the objective function
functionID = 2;

% Defines the line search method
ls_method = 'S';

% Defines the stopping criteria
SC_index = 2;

% Defines if the surface plot of f is shown or not
s_plot = false;

% Defines if the plots are saved or not (ARNAUD)
save_plot = false;

% Initial point 
xinit_values   = [10, 10; -10, -10; 10, -10; -10, 10];

% Maximum number of iterations
MaxIter = 300;

% Maximum number of iterations to compute alpha
MaxIter_alpha = 100;    

% Tolerances for the stoping criteria 
Epsilon = linspace(2, 20, 19);
Nu      = linspace(2, 20, 19);    

%% --------------
%  Initialization
%  --------------
% Stores all the parameters to be displayed on command window
parameters = [0, 10, 10, MaxIter, Epsilon, Nu, SC_index, 0];

% Retreiving the method trough CW
method = 1;

% Updating the method used
parameters(1) = method;



%% --------
%  Symbolic
%  --------
syms x1 x2 alpha;

% Definition of a symbolic objective f, gradient and Hessian matrix.
X         = [x1 x2];
f(x1, x2) = getObjF(X, functionID);
grad_f    = gradient(f);
H_f       = hessian(f);

% Stores the total number of iterations made by the compu. of alpha
alpha_iters = 0;
   
figure();

for xi = 1 : 4

   xinit = [xinit_values(xi, 1), xinit_values(xi, 2)];

    % Further initialization
    n       = 2;                    % Dimension of the problem
    xinit   = reshape(xinit, 2, 1); % To be sure that it's a column vector
    x       = zeros(n, MaxIter);    % Initialization of vector x
    x(:, 1) = xinit;                % Put xinit in vector x.

    % Contains the number of iterations for each epsilon
    nb_iterations = zeros(size(Epsilon));

    for t = 1 : size(Epsilon, 2)

        cur_eps = Epsilon(t)

        for i = 1 : MaxIter
        
            % 1 - Finding steepest descent direction                    
            s = -[grad_f(x(1, i), x(2, i))];
        
            % 2 - Convergence check (/!\ SC evaluated at grad(x_(k + 1)))
            if i ~= 1 && stoppingCriteria(SC_index, s, 10^-cur_eps, f, Nu, x(:, i - 1), x(:, i))
                parameters(8) = i;
                break;
            end
        
            % 3 - Computing alpha
            phi(alpha)             = f(x(1, i) + alpha * s(1), x(2, i) + alpha * s(2));
            [alpha_opt, alpha_it]  = find_alpha(phi, ls_method, method, MaxIter_alpha, 0.1, i, H_f, s, s);
            
            % Updating the number of iterations to compute alpha
            alpha_iters = alpha_iters + alpha_it;
        
            % 4 - Updating x
            x(1, i + 1) = x(1, i) + alpha_opt * s(1);
            x(2, i + 1) = x(2, i) + alpha_opt * s(2);    
        end
        
        x = x(:, 1 : i);   

        if parameters(8) == 0
            parameters(8) = MaxIter;
        end
            
        % Updating iters
        nb_iterations(t) = parameters(8);
        double(f(x(1, end), x(2,end)))

    end

    if xi == 1
        semilogy(Epsilon, nb_iterations, '-', 'Color',  [0 0.4470 0.7410], 'LineWidth', 2)
    elseif xi == 2
        semilogy(Epsilon, nb_iterations, '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
    elseif xi == 3
        semilogy(Epsilon, nb_iterations, ':', 'Color',  [0 0.4470 0.7410],'LineWidth', 2)
    else
        semilogy(Epsilon, nb_iterations, '-.', 'Color', [0 0.4470 0.7410],'LineWidth', 2)
    end
    hold on;
    grid on;
end

xlabel('-log($\varepsilon$) [-]', 'FontSize', 18, 'interpreter', 'Latex');
ylabel("Number of iterations [-]", 'Fontsize', 18);
















