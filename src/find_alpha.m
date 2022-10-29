function [alpha, iter, calls_to_f] = find_alpha(phi, method, opti_method, max_it, h, k, H, s, d)
    %--------------
    % Documentation
    %--------------
    % Computes the value of alpha where:
    % - method      = 'NR'  --> Newton-Raphson
    %                 'S'   --> Secant
    %                 'D'   --> Dichotomy
    %                 'BB'  --> Black-box 'solve' function @Mathworks
    %                 'DIV' --> Use divergent serie ~ 1/(k+1) for alpha_k
    %                 'CQ'  --> Formula for convex quadratic functions
    %
    % - opti_method = index of optimization method used (needed for "CQ")
    %
    % - h  : the step used for finding bounds for Secant and Dichotomy methods
    % - k  : the iteration at which find_alpha is called
    % - H  : the Hessian matrix of f
    % - s  : the gradient of f at the considered point
    % - d  : the direction of f at the considered point (used in conjugate
    %        gradient method) Note: d = s in the steepest descent method,
    %        it's important in the case the method used to compute alpha is
    %        the "CQ".
    %
    % Stores the number of iterations to compute alpha as well as the
    % number of calls to the objective function f
    iter       = 0;
    calls_to_f = 0;

    switch method

        case 'NR'
            % First and second derivatives of phi
            dphi = diff(phi);
            d2phi = diff(dphi);      

            % Initialization of alpha's
            alpha0 = 0;
            previous_alpha = alpha0;
            alpha = 1;
            
            % Iterative update (adding maximum number of it in loop condition to reduce computation time)
            while(abs(alpha - previous_alpha) > 1e-5 && iter < max_it)
                previous_alpha = alpha0;
                % prevent from division by zero
                if(vpa(d2phi(alpha0)) == 0)
                    alpha0 = alpha0 - vpa(dphi(alpha0))/(vpa(d2phi(alpha0))+0.01);
                else
                    alpha0 = alpha0 - vpa(dphi(alpha0))/vpa(d2phi(alpha0));
                end
                alpha = alpha0;
                iter = iter + 1;
                calls_to_f = calls_to_f + 2;                               % 2 calls : dphi and d2phi
            end

        case 'D'
            % Derivative of phi
            dphi = diff(phi);

            % Initialization of alpha's
            alpha_k = h;
            while(dphi(alpha_k) < 0)
                alpha_k = alpha_k + h;
                calls_to_f = calls_to_f + 1;                               % 1 call : dphi
            end

            % Iterative update
            alpha_kminus1 = 0;
            alpha = 0.5*(alpha_k + alpha_kminus1);
            while(max( abs(alpha - alpha_k) ,...
                       abs(alpha - alpha_kminus1) ) > 1e-5)
                   
                if(vpa(dphi(alpha)) < 0)
                    alpha_kminus1 = alpha;
                else
                    alpha_k = alpha;
                end
                calls_to_f = calls_to_f + 1;                               % 1 call : dphi in if condition

                alpha = 0.5*(alpha_k + alpha_kminus1);

                iter = iter + 1;
            end

        case 'S'
            % Derivative of phi
            dphi = diff(phi);

            % Initialization of alpha's
            alpha_k = h;
            while(dphi(alpha_k) < 0)
                alpha_k = alpha_k + h;
                calls_to_f = calls_to_f + 1;                               % 1 call : dphi 
            end
            alpha = alpha_k/2;

            % Iterative update
            alpha_kminus1 = 0;
            while(max(abs(alpha - alpha_k) , abs(alpha - alpha_kminus1)) > 1e-5)
                dphi_alpha_k = dphi(alpha_k); % Save 1 call

                % Division by 0 for f2 using the CG in (10, 10)
                denominator = vpa(dphi_alpha_k) - vpa(dphi(alpha_kminus1));
                if denominator == 0
                    denominator = 0.1;
                end

                alpha = alpha_k - vpa(dphi_alpha_k) * (alpha_k - alpha_kminus1)/( denominator );
                alpha_kminus1 = alpha_k;
                alpha_k = alpha;
                iter = iter + 1;
                calls_to_f = calls_to_f + 3;                               % 2 calls : dphi_alpha_k and dphi(alpha_kminus1)
            end

        case 'BB'
            alpha = solve(diff(phi) == 0);
            alpha = min(alpha([isAlways(alpha>0)] )); % Take the minimum positive alpha
            iter = iter + 1;
            calls_to_f = 0;                                                % 1 call : dphi even if we don't know exactly since BB

        case 'DIV'
            alpha = 1/(k + 1);
            iter = iter + 1;

        case 'CQ'
            % Note: (-1)^(opti_method + 1), there is a + if we are in the
            % SDM and -1 in the case of the CG method.
            alpha = (-1)^(opti_method + 1) * (transpose(d) * s)./(transpose(d) * H * d) ;
            iter = iter + 1;
    end
end












