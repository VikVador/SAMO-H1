%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               Structural and Multidisciplinary Optimization             %
%                                                                         %
%                      H1 - Unconstrained Optimization                    %
%                                                                         %
% @ Arnaud Rémi                                        @ Victor Mangeleer %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Statement:
% ----------
%   Solve the minimization problem
%
%             min f(x,y)
%
%   using the following optimization methods
%
%   1) Steepest descent
%   2) Conjugate gradients with Fletcher-Reeves update rule
%   3) BFGS Quasi-Newton
%
% Code wiki:
%-----------
% - For the linesearch method:
%      
%   Newton raphson            ---> 'NR'
%   Secant                    ---> 'S'
%   Dichotomy                 ---> 'D'
%   Black box                 ---> 'BB'
%   Divergent serie           ---> 'DIV'
%   Convex quadratic function ---> 'CQ'
%
% - For the stopping criteria
%
%   max(grad_f)         < eps
%   norm(grad_f, 2)     < eps
%   f(x_(k+1)) - f(x_k) < nu
%

clearvars; close all; clc
%% ----------
%  Parameters
%  ----------
% Defines the objective function
functionID = 2;

% Defines the line search method
ls_method = 'D';

% Defines the stopping criteria
SC_index = 2;

% Defines if the surface plot of f is shown or not
s_plot = true;

% Initial point 
xinit   = randi([-10 10], 1, 2);

% Maximum number of iterations
MaxIter = 100;

% Tolerances for the stoping criteria 
Epsilon = 1e-8;         
Nu      = 1e-4;       

%% --------------
%  Initialization
%  --------------
% Stores all the parameters to be displayed on command window
parameters = [0, xinit(1), xinit(2), MaxIter, Epsilon, Nu, SC_index, 0];

% Retreiving the method trough CW
method = terminal(1, parameters, ls_method);

% Updating the method used
parameters(1) = method;

% Displaying all the optimization parameters
terminal(2, parameters, ls_method)

% Further initialization
n       = 2;                    % Dimension of the problem
xinit   = reshape(xinit, 2, 1); % To be sure that it's a column vector
x       = zeros(n, MaxIter);    % Initialization of vector x
x(:, 1) = xinit;                % Put xinit in vector x.

%% --------
%  Symbolic
%  --------
syms x1 x2 alpha;

% Definition of a symbolic objective f, gradient and Hessian matrix.
X         = [x1 x2];
f(x1, x2) = getObjF(X, functionID);
grad_f    = gradient(f);
H_f       = hessian(f);

%% --------
%  Plotting
%  --------
% Plot the 3D surface of the objective function being optimized
%
terminal(3, parameters, ls_method);

if s_plot
    plotObjF(f, functionID);
    waitforbuttonpress();
end

%% ------------
%  Optimization
%  ------------
terminal(4, parameters, ls_method);

% Stores the total number of iterations made by the compu. of alpha
alpha_iters = 0;

switch method
    case 1        
        for i = 1 : MaxIter

            % 1 - Finding steepest descent direction                    
            s = -[grad_f(x(1, i), x(2, i))];
          
            % 2 - Computing alpha (/!\ divergent serie is not a LS method)
            phi(alpha) = f(x(1, i) + alpha * s(1), x(2, i) + alpha * s(2));
            [alpha_opt, alpha_it]  = find_alpha(phi, ls_method, 0.1, i, H_f, s);
            alpha_iters = alpha_iters + alpha_it;

            % 3 - Updating x
            x(1,i+1) = x(1,i) + alpha_opt * s(1);
            x(2,i+1) = x(2,i) + alpha_opt * s(2);
            

                
            % 4 - Convergence check
            if stoppingCriteria(SC_index, s, Epsilon, f, Nu, x(:,i), x(:,i+1))
                parameters(8) = i; % Count iterations
                break;
            end

        end
        
        x = x(:, 1 : i); 
        
    case 2        
        for i = 2 : MaxIter

            % Conjugate gradients with Fletcher-Reeves update rule

        end
        
        x = x(:, 1 : i); 
        
    case 3        
        for i = 2 : MaxIter
            
            % BFGS Quasi-Newton

        end
        
        x = x(:, 1 : i);      
end

%% ----------------------------
%  Plotting and showing results
%  ----------------------------
terminal(5, parameters, ls_method);

disp("Iterations (Method Optimization) : " + int2str(parameters(8)));
disp(" ");
disp("Iterations (Alpha Computation)   : " + int2str(alpha_iters));
disp(" ");
disp("Iterations (Total)               : " + int2str(alpha_iters + parameters(8)));
disp(" ");
fprintf('Optimal solution                 : (x = %4.3f, y = %4.3f)\n', x(1, end), x(2, end));
disp(" ");
fprintf('Objective function value         : %4.3f.\n' , f(x(1, end), x(2, end)));
disp(" ");
disp(" ");
disp(" ");

% Plotting optimization path (2D & 3D) and evolution of f
plotEvolutionObjF(x, f, functionID, parameters)
plotOptimizationPath(x, f, functionID, parameters);
plotOptimizationPath3D(x, f, functionID, parameters);
















