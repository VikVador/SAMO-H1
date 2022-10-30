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
functionID = 1;

% Maximum number of iterations
MaxIter = 100;

% Maximum number of iterations to compute alpha
MaxIter_alpha = 100;

% Step size for mesh
n = 0.5;

% Size of the square domain ([-bound, bound])
bond = 10;

% Precision of the classification in %
class_precision = 5;

% NOTE: this mean if the optimize value of f found is inside the interval
%       [f_real_opti * (1 - precision/100), f_real_opti * (1 + precision/100)]
%       the value will be considered as an optimum.

%% --------------
%  Initialization
%  --------------
method    = 1;
ls_method = 'NR';
SC_index  = 2;
Epsilon   = 1e-5;
Nu        = 1e-4;

% Initialization of the mesh
[x_dom, y_dom]       = meshgrid(-bond : n : bond);
[f_value, categorie] = meshgrid(-bond : n : bond);

%% --------
%  Symbolic
%  --------
syms x1 x2 alpha;

% Definition of a symbolic objective f, gradient and Hessian matrix.
X         = [x1 x2];
f(x1, x2) = getObjF(X, functionID);
grad_f    = gradient(f);
H_f       = hessian(f);

% Opti values
f1_opt       = -22.143;
f2_opt_local = -0.138;
f2_opt_abs   = -2.038;

% For the classification boundaries
classificator_u = 1 - class_precision/100;
classificator_d = 1 + class_precision/100;

%% ------------
%  Optimization
%  ------------
% Displaying information over terminal (1)
disp("-------------------");
disp("Computing solutions");
disp("-------------------");

% Counts the total number of iteration
nb_iterations = 0;

% Measure timing
time = 0;

% Looping over the domain
for xi = 1 : size(x_dom, 1)

    % Displaying information over terminal (2)
    disp(int2str(xi) + "/" + int2str(size(x_dom, 1)) + " - Time left = " + ...
         num2str(time * (size(x_dom, 1) - xi)) + " [s] = " + ...
         num2str(time * (size(x_dom, 1) - xi)/60) + " [min]");
    tic;

    for yi = 1 : size(y_dom, 1)

        % New initial point
        xinit   = [x_dom(xi, yi) y_dom(xi, yi)];

        %------------------------------------------------------------------
        %                                SDM
        %------------------------------------------------------------------
        n       = 2;
        xinit   = reshape(xinit, 2, 1);
        x       = zeros(n, MaxIter);
        x(:, 1) = xinit;
        for i = 1 : MaxIter
            s = -[grad_f(x(1, i), x(2, i))];
            if i ~= 1 && stoppingCriteria(SC_index, s, Epsilon, f, Nu, x(:, i - 1), x(:, i))
                nb_iterations = i;
                break;
            end
            phi(alpha) = f(x(1, i) + alpha * s(1), x(2, i) + alpha * s(2));
            [alpha_opt, alpha_it]  = find_alpha(phi, ls_method, method, MaxIter_alpha, 0.1, i, H_f, s, s);
            x(1, i + 1) = x(1, i) + alpha_opt * s(1);
            x(2, i + 1) = x(2, i) + alpha_opt * s(2);
        end
        x = x(:, 1 : i);
        %------------------------------------------------------------------
        %                                SDM
        %------------------------------------------------------------------
        % Updating init value of f (to create starting point regions)
        f_value(xi, yi) = f(xinit(1), xinit(2));

        % Computing final value of f for classification
        x_end = x(1, end); y_end = x(2, end);
        f_end_value = double(f(x_end, y_end));

        % Classification of f1 (abs or impossible to converge)
        if functionID == 1
            if (f1_opt * classificator_d) < f_end_value && f_end_value < (f1_opt * classificator_u)
                categorie(xi, yi) = 1;
            else
                categorie(xi, yi) = 0;
            end
        end

        % Classification of f (abs, local or impossible to converge)
        if functionID == 2
            if (f2_opt_abs * classificator_d) < f_end_value && f_end_value < (f2_opt_abs * classificator_u)
                categorie(xi, yi) = 1;
            elseif (f2_opt_local * classificator_d) < f_end_value && f_end_value < (f2_opt_local * classificator_u)
                categorie(xi, yi) = 2;
            else
                categorie(xi, yi) = 0;
            end
        end
    end

    % Update time
    time = toc;
end

%% --------------------------------------
%  Plotting the domain of starting points
%  --------------------------------------
% Defining possible dot colors
colors = [1 0 0;                 % impossible to converge = 0
          0.4660 0.6740 0.1880;  % abs
          0.9290 0.6940 0.1250]; % local

figure();

% Looping over all possible results (abs = 1, local = 2 or impossible to converge = 0)
for xi = 1 : size(x_dom, 1)
    for yi = 1 : size(x_dom, 2)

        % Skip (0, 0)
        if x_dom(xi, yi) == 0 && y_dom(xi, yi) == 0
            continue
        end

        % Current categorie
        categ = categorie(xi, yi) + 1;

        % Current color
        color = colors(categ, :);

        % Plotting point
        scatter3(x_dom(xi, yi), y_dom(xi, yi), f_value(xi, yi), ...
        'MarkerEdgeColor','k', 'MarkerFaceColor', color, "SizeData", 200);
        hold on;
    end
end

% Adding axes and legend
xlabel('x_{init} [-]');
ylabel('y_{init} [-]');
zlabel("f_" + int2str(functionID) + "(x, y)");
set(gca,'fontsize', 20);
set(gcf,'position',[50, 50, 750, 600]);
