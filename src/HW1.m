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
functionID = 1;

% Defines the line search method
ls_method = 'NR';

% Defines the stopping criteria
SC_index = 2;

% Defines if the surface plot of f is shown or not
s_plot = true;

% Defines if the plots are saved or not (ARNAUD)
save_plot = true;

% Initial point 
xinit   = [10, 10];

% Maximum number of iterations
MaxIter = 100;

% Maximum number of iterations to compute alpha
MaxIter_alpha = 100;    

% Tolerances for the stoping criteria 
Epsilon = 1e-5;         
Nu      = 1e-5;       

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
    plotObjF(f, functionID, save_plot);
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

            % 2 - Convergence check (/!\ SC evaluated at grad(x_(k + 1)))
            if i ~= 1 && stoppingCriteria(SC_index, s, Epsilon, f, Nu, x(:, i - 1), x(:, i))
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
        

    case 2 

        % 1 - Initialization of gradient and direction
        g1 = [grad_f(x(1, 1), x(2, 1))];
        d  = -g1;

        for i = 1 : MaxIter
            
            % 2 - Computing alpha
            phi(alpha)             = f(x(1, i) + alpha * d(1), x(2, i) + alpha * d(2));
            [alpha_opt, alpha_it]  = find_alpha(phi, ls_method, method, MaxIter_alpha, 0.1, i, H_f, g1, d);
            alpha_iters            = alpha_iters + alpha_it;
            
            % 3 - Updating x
            x(1, i + 1) = x(1, i) + alpha_opt * d(1);
            x(2, i + 1) = x(2, i) + alpha_opt * d(2);   
            
            % 4 - Computing new gradient
            g0 = g1;
            g1 = [grad_f(x(1, i + 1), x(2, i + 1))];

            % 5 - Convergence check
            if i ~= 1 && stoppingCriteria(SC_index, g1, Epsilon, f, Nu, x(:, i), x(:, i + 1))
                parameters(8) = i;
                break;
            end

            % 6 - Computing beta using Method of Fletcher and Reeves
            beta = vpa(norm(g1, 2)/norm(g0, 2))^2;

            % 7 - Computing new direction
            d(1) = -g1(1) + beta * d(1);
            d(2) = -g1(2) + beta * d(2);

        end
        
        x = x(:, 1 : i+1); 
        
    case 3        
        H = eye(2);
        for i=1:MaxIter
            if i==1
                    g(:,i) = [grad_f(x(1,i),x(2,i))];
                    phi(alpha) = f(x(1,i)-alpha*g(1) , x(2,i)-alpha*g(2));
                    alpha_opt = find_alpha(phi, ls_method, method, MaxIter_alpha, 0.1, i, H_f, -g(:,i), -g(:,i));
                    x(1,i+1) = x(1,i) - alpha_opt*g(1);
                    x(2,i+1) = x(2,i) - alpha_opt*g(2);
                continue; 
            end
            
            g(:,i) = [grad_f(x(1,i),x(2,i))];
            gamma_k = g(:,i) - g(:,i-1);
            delta_k = x(:,i) - x(:,i-1);
            
            if i ~= 1 && stoppingCriteria(SC_index, g(:,i), Epsilon, f, Nu, x(:, i - 1), x(:, i))
               parameters(8) = i;
               break;
            end
            
            % Terms in H's update
            a = 1 + ( transpose(gamma_k) * H * gamma_k )/...
                    (transpose(delta_k)*gamma_k);
            A = a*(delta_k * transpose(delta_k))/(transpose(delta_k)*gamma_k);
            B = ( delta_k * transpose(gamma_k) * H...
                     +H * gamma_k * transpose(delta_k) )/...
                                    (transpose(delta_k)*gamma_k);
            H = H + A - B;
            d = H*g(:,i);
            if(g(:,i) ~= [0;0])
                phi(alpha) = f(x(1,i)-alpha*d(1) , x(2,i)-alpha*d(2));
                alpha_opt = find_alpha(phi, ls_method, method, MaxIter_alpha, 0.1, i, H_f, g(:,i), d); 
            else
                alpha_opt = 0;
            end
            
            x(:,i+1) = x(:,i) - alpha_opt*H*g(:,i);

        end
        
        x = x(:, 1 : i);      
end

%% ----------------------------
%  Plotting and showing results
%  ----------------------------
terminal(5, parameters, ls_method);

% Fixing iteration value due to new position of SC in the loop
if parameters(8) == 0
    parameters(8) = MaxIter;
end

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
plotEvolutionObjF(x, f, functionID, parameters, save_plot)
plotOptimizationPath(x, f, functionID, parameters, save_plot);
plotOptimizationPath3D(x, f, functionID, parameters, save_plot);













