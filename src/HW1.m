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
clearvars; close all; clc

%% ----------
%  Parameters
%  ----------
% Defines the objective function that will be used for the optimization
functionID = 1;

% Defines the used line search method
ls_method = 'NR';

% Defines the stopping criteria that will be used during optimization
SC_index = 2;

% Defines if the surface plot are shown or not
s_plot = false;

% Initial point 
xinit   = randi([-10 10], 1, 2);

% Maximum number of iterations
MaxIter = 100;

% Tolerance for the stoping criteria (1 & 2)
Epsilon = 1e-5;         

% Tolerance for the stoping criteria (3)
Nu = 1e-4;       

%% --------------
%  Initialization
%  --------------
% Stores all the parameters to be displayed on command window
parameters = [0, xinit(1), xinit(2), MaxIter, Epsilon, Nu, SC_index, 0];

% Retreiving the method
method = terminal(1, parameters);

% Updating the method used
parameters(1) = method;

% Displaying all the optimization parameters
terminal(2, parameters)

% Further initialization
n       = 2;                    % Dimension of the problem
xinit   = reshape(xinit, 2, 1); % To be sure that it's a column vector
x       = zeros(n, MaxIter);    % Initialization of vector x
x(:, 1) = xinit;                % Put xinit in vector x.

%% --------
%  Symbolic
%  --------
syms x1 x2 alpha;

X         = [x1 x2];
f(x1, x2) = getObjFVal(X, functionID);
grad_f    = gradient(f);
H_f       = hessian(f);

%% --------
%  Plotting
%  --------
% This section has for purpose to plot the 3D surface for f1 and f2
%
terminal(3, parameters);

if s_plot
    plotObjF(f);
    waitforbuttonpress();
    plotObjF(f);
    waitforbuttonpress();
end

%% ------
%  Method
%  ------
terminal(4, parameters);

switch method
    case 1        
        for i=1:MaxIter
            %---
            % Finding steepest descent direction
            %---           
            s = -[grad_f(x(1,i),x(2,i))];
          
            %---
            % Line search methods
            %---
            phi(alpha) = f(x(1,i)+alpha*s(1) , x(2,i)+alpha*s(2));
            alpha_opt = find_alpha(phi, ls_method, 0.1, i, H_f, s);
                        
            %---
            % Update x
            %---
            x(1,i+1) = x(1,i) + alpha_opt*s(1);
            x(2,i+1) = x(2,i) + alpha_opt*s(2);

            %---
            % Check convergence
            %---
            if stoppingCriteria(SC_index, s, Epsilon, f, Nu, x(:,i), x(:,i+1))
                % Count iterations
                parameters(8) = i;
                break;
            end

        end
        
        x = x(:, 1 : i); 
        
    case 2        
        for i = 2 : MaxIter
            %%%% ---------------------------------------------------------------------------
            %%%% ADD YOUR CODE
            %%%%
        end
        
        x = x(:, 1 : i); 
        
    case 3        
        for i = 2 : MaxIter
            %%%% ---------------------------------------------------------------------------
            %%%% ADD YOUR CODE
            %%%%
        end
        
        x = x(:, 1 : i);      
end

%% ----------------------------
%  Plotting and showing results
%  ----------------------------
% Adding the actual number of iterations
terminal(5, parameters);

disp("Number of iterations     : " + int2str(parameters(8)));
disp(" ");
fprintf('Optimal solution         : (x = %4.3f, y = %4.3f)\n', x(1, end), x(2, end));
disp(" ");
fprintf('Objective function value : %4.3f.\n' ,getObjFVal(x(:, end),functionID));
disp(" ");
disp(" ");
disp(" ");

% Plotting the optimization path (2D & 3D) and saving it
plotOptimizationPath(x, functionID, parameters);
plotOptimizationPath3D(x, functionID, parameters);
















