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

% Maximum number of iterations
MaxIter = 20;

% Step size for mesh
n = 0.5;

% Boundary
bond = 10;

% Tolerances for the stoping criteria 
Epsilon = 1e-4;         
Nu      = 1e-4;  

% Precision for classification in %
class_precision = 5;

%% --------------
%  Initialization
%  --------------
method    = 1;
ls_method = 'NR';
SC_index  = 2;

% Initialization of the mesh
[x_dom, y_dom]       = meshgrid(-bond:n:bond);
[f_value, categorie] = meshgrid(-bond:n:bond);

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
% Counts the total number of iteration
nb_iterations = 0;

% Looping over the domain
for xi = 1 : size(x_dom, 1)
    for yi = 1 : size(y_dom, 1)
            
        % New initial point
        xinit   = [x_dom(xi, yi) y_dom(xi, yi)]        

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
            [alpha_opt, alpha_it]  = find_alpha(phi, ls_method, 0.1, i, H_f, s);
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
        
        % Classification of f (abs, local or impossible to converge)
        if functionID == 1
            % Opti abs
            if (f1_opt * classificator_d) < f_end_value && f_end_value < (f1_opt * classificator_u)
                categorie(xi, yi) = 1;
            else
                categorie(xi, yi) = 0;
            end
        end


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
end

%%
close all;
figure();

% Plotting all the subsurfaces
for i = 0 : 2
    
    % Creation of sub-matrices
    x_sub = zeros(size(x_dom));
    y_sub = zeros(size(y_dom));
    f_sub = zeros(size(f_value));

    % Applying mask
    x_sub(categorie == i) = x_dom(categorie == i);
    y_sub(categorie == i) = y_dom(categorie == i);
    f_sub(categorie == i) = f_value(categorie == i);

    % Colors
    colors = [1 0 0; 0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250];

    surf(x_sub, y_sub, f_sub, 'EdgeColor','none', 'FaceColor', colors(i + 1, :))
    hold on;

end

% Adding axes
xlabel('x [-]');
ylabel('y [-]');
zlabel("f_" + int2str(functionID) + "(x, y)");

% Adding position of optimums
if functionID == 1
    scatter3(-3.144, 4.859, -22.143 + 0.1,'o',"SizeData", 200)
else
    scatter3(0.922,  0, -0.138 + 0.1,'o',"SizeData", 200)
    scatter3(-1.232, 0, -2.038 + 0.1,'o',"SizeData", 200)
end
















