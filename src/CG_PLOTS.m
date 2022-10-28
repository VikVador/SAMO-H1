%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               Structural and Multidisciplinary Optimization             %
%                                                                         %
%                      H1 - Unconstrained Optimization                    %
%                                                                         %
% @ Arnaud Rémi                                        @ Victor Mangeleer %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Documentation
% -------------
% This script allows us to generate all the graphs we need for the report

clearvars; close all; clc
%% ----------
%  Parameters
%  ----------
% Defines the objective function
functionID = 1;

% Initial points
xinit_values = [10 10; -10 -10]; 

% Maximum number of iterations in the main loop
MaxIter = 14;    

% Maximum number of iterations to compute alpha
MaxIter_alpha = 50;    

% Defining stopping criteria and alpha used for our comparison
SC_values    = [1, 2, 3];
ALPHA_values = ["NR", "S", "D", "BB", "DIV", "CQ"];

% Tolerances for the stoping criteria 
Epsilon = 1e-5;         
Nu      = 1e-5; 

%  --------
%  Symbolic
%  --------
syms x1 x2 alpha;

% Definition of a symbolic objective f, gradient and Hessian matrix.
X         = [x1 x2];
f(x1, x2) = getObjF(X, functionID);
grad_f    = gradient(f);
H_f       = hessian(f);

%  ------------------------------------------------------------------------
%                         Others (don't need to look)
%  ------------------------------------------------------------------------
% Stores the number of iterations/calls during the optimization process
iter_call   = zeros(6, 3);

% For the plots and terminal
alpha_name = ["Newton raphson"; "Secant"; "Dichotomy"; "Black Box"; "Divergent serie"; "Convex quadratic function"];
sc_name    = ["max(abs(grad_f)) < eps"; "norm(grad_f, 2) < eps"; "f(x_(k+1)) - f(x_k) < nu"];

% Information over terminal (1)
disp("Epsilon : " + sprintf('%.10f', Epsilon));
disp(" ");
disp("Nu : " + sprintf('%.10f', Nu));
disp(" ");

%% ------------------------------------------------------------------------
%                                Optimization
%  ------------------------------------------------------------------------
% Looping over stopping criteria
for s_val = 1 : size(SC_values, 2)

    % Current stopping criteria observed
    SC_index = SC_values(s_val);

    % Information over terminal (1)
    disp(" "); disp("SC = " + sc_name(s_val)); plt = figure();

    % Looping over init points
    for xi = 1 : 2
    
        % Current initial point value for the optimization
        xinit = [xinit_values(xi, 1), xinit_values(xi, 2)];

        % Information over terminal (2)
        disp("X_0 : (" + int2str(xinit(1)) + ", " + int2str(xinit(2)) + ")"); disp(" ");

        % Looping over alpha computation criteria
        for a = 1 : size(ALPHA_values, 2)

            % Current method used to compute alpha
            ls_method = ALPHA_values(a);

            % Stores the number of iterations/calls during the optimization process
            alpha_iters = 0;
            opti_iters  = 0;
            f_calls     = 0;

            %--------------------------------------------------------------
            %                           SECURITY
            %--------------------------------------------------------------
            % Security (1) - Do not apply CQ on second f
            if ls_method == "CQ" && functionID == 2
                continue;
            end

            %--------------------------------------------------------------
            %                           SECURITY
            %--------------------------------------------------------------
            % Information over terminal (2)
            disp("|--> Method = " + alpha_name(a));

            % Initialization
            n       = 2;                    
            xinit   = reshape(xinit, 2, 1); 
            x       = zeros(n, MaxIter);   
            
            % 1 - Initialization of gradient and direction
            g1 = [grad_f(x(1, 1), x(2, 1))];
            d  = -g1;
    
            for i = 1 : MaxIter
                
                % 2 - Computing alpha
                phi(alpha)                        = f(x(1, i) + alpha * d(1), x(2, i) + alpha * d(2));
                [alpha_opt, alpha_it, f_call_it]  = find_alpha(phi, ls_method, 2, MaxIter_alpha, 0.1, i, H_f, g1, d);
                
                % Updating the number of iterations
                alpha_iters = alpha_iters + alpha_it;
                f_calls     = f_calls + f_call_it;
    
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
                beta = norm(g1, 2)/norm(g0, 2);
    
                % 7 - Computing new direction
                d(1) = -g1(1) + beta * d(1);
                d(2) = -g1(2) + beta * d(2);

                % Updating the number of iterations
                if SC_index == 3
                    f_calls = f_calls + 1 + 2 + 1;
                else
                    f_calls = f_calls + 1 + 1;
                end

                % We don't know how many calls for the BB
                if ls_method == "BB"
                    f_calls = 0;
                end
            end

            x = x(:, 1 : i); 
                
            %  ------------------------------------------------------------
            %                   Others (don't need to look)
            %  ------------------------------------------------------------
            % If MaxIter reached, the value of opti_iters is not updated !
            % Thus, we do it manually here
            if opti_iters == 0
                opti_iters = MaxIter;
            end
            
            % Updating number of iterations
            iter_call(a, 1) = opti_iters;
            iter_call(a, 2) = alpha_iters;
            iter_call(a, 3) = f_calls;

            % Information over terminal (3)
            fprintf('   |--> X_opti = (x = %4.3f, y = %4.3f)\n', x(1, end), x(2, end));
            fprintf('   |--> f_opti = %4.3f.\n' , f(x(1, end), x(2, end)));

            % -------------------------------------------------------------
            %       Plotting (1) - Evolution curve of f w.r.t nb_iter
            % -------------------------------------------------------------
            % Number of optimization iterations to get there
            iterations = linspace(1, size(x, 2), size(x, 2));
        
            % Contains all the colors that will be used to plot ! Indeed,
            % since we are plotting the evolution curve for all methods for
            % both starting points on A SAME GRAPH, we needed to make sur
            % that the same colors were used for the 2 curves ('-' and
            % '--') representing the same method but at different starting
            % points ! This is the way that I found working !
            colormap = ["ebdc78", "ffb55a", "4421af", "0d88e6", "77dd77", "b30000"];

            % Converts hexadecimal value to rgb triplet ! Matlab does not
            % know how to read hexadecimal :'(
            cur_color = sscanf(colormap(a),'%2x%2x%2x',[1 3])/255;
            
            % Plotting evolution curve (shape differs depends on starting point)
            if xi == 1
                plot(iterations, f(x(1, :), x(2, :)), "-", 'Color', cur_color, 'LineWidth', 2);
            else
                plot(iterations, f(x(1, :), x(2, :)), "--", 'Color', cur_color, 'LineWidth', 2);
            end

            % Make sure to add all curves on same plot
            hold on;
        end
    end
    % -------------------------------------------------------------
    %       Plotting (1) - Making graph look pretty
    % -------------------------------------------------------------
    xlabel('Number of iterations [-]', 'FontSize', 18);
    ylabel("f_" + int2str(functionID) + "(x, y)", 'Fontsize', 18);
    if functionID == 1
        h = legend('Newton-Raphson','Secant', 'Dichotomy', 'Black Box', 'Divergent serie', 'Convex quadratic function');
    else
        h = legend('Newton-Raphson','Secant', 'Dichotomy', 'Black Box', 'Divergent serie');
    end
    set(h,'FontSize', 18);
    set(gca,'fontsize', 18);
    grid on;
    xlim([1, MaxIter]);
    saveas(plt, "../graphs/report/CG/" + int2str(functionID) + "/f" + int2str(functionID) + "_" + ...
                sc_name(s_val) + "_(" + int2str(xinit_values(1,1)) + "," + int2str(xinit_values(1,2)) + ")_" + ...
                "_(" + int2str(xinit_values(2,1)) + "," + int2str(xinit_values(2,2)) + ")_" + ".pdf");

    % -------------------------------------------------------------
    %     Plotting (2) - Bar plot of iterations and calls to f
    % -------------------------------------------------------------
    plt = figure();
    b = bar(iter_call);
    ylim([0, max(max(iter_call)) * 1.1]) 
    % Note : Increase by 10 % the vertical height for no 
    % overlap between text and the top of the plot figure

    % Adding value on top of the bars ! This is for bar 1
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(b(1).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', 15)

    % This is for bar 2
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YEndPoints;
    labels2 = string(b(2).YData);
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', 15)
   
    % This is for bar 3
    xtips3 = b(3).XEndPoints;
    ytips3 = b(3).YEndPoints;
    labels3 = string(b(3).YData);
    text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom', 'FontSize', 15)

    % Making the plot good looking
    set(h,'FontSize', 18);
    set(gca,'fontsize', 18);

    % Adding method name as labels
    set(gca, 'XTickLabel',{'Newton-Raphson','Secant', 'Dichotomy', 'Black Box', 'Divergent serie', 'Conv. quad. fun.'}, 'FontSize', 12)
    xtickangle(20);
    ylabel('Number of iterations/calls to f [-]', 'FontSize', 20);
    grid on;
    saveas(plt, "../graphs/report/CG/" + int2str(functionID) + "/bar_f" + int2str(functionID) + "_" + ...
                sc_name(s_val) + "_(" + int2str(xinit_values(1,1)) + "," + int2str(xinit_values(1,2)) + ")_" + ...
                "_(" + int2str(xinit_values(2,1)) + "," + int2str(xinit_values(2,2)) + ")_" + ".pdf");
end















