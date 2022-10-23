function [alpha, iter] = find_alpha(phi, method, h, k, H, s)
    %--------------
    % Documentation
    %--------------
    % Computes the value of alpha where:
    % - method = 'NR'  --> Newton-Raphson
    %            'S'   --> Secant
    %            'D'   --> Dichotomy
    %            'BB'  --> Black-box 'solve' function @Mathworks
    %            'DIV' --> Use divergent serie ~ 1/(k+1) for alpha_k
    %            'CQ'  --> Formula for convex quadratic functions
    %
    % - h : the step used for finding bounds for Secant and Dichotomy methods
    % - k : the iteration at which find_alpha is called
    % - H : the Hessian matrix of f
    % - s : the gradient of f at the considered point
    % 
    % Stores the number of iterations
    iter = 0;

    switch method

        case 'NR'
            % First and second derivatives of phi
            dphi = diff(phi);
            d2phi = diff(dphi);      

            % Initialization of alpha's
            alpha0 = 0;
            previous_alpha = alpha0;
            alpha = 1;

            % Iterative update
            while(abs(alpha - previous_alpha) > 1e-5)
                previous_alpha = alpha0;
                alpha0 = alpha0 - vpa(dphi(alpha0))/vpa(d2phi(alpha0));
                alpha = alpha0;
                iter = iter + 1;
            end

        case 'D'
            % Derivative of phi
            dphi = diff(phi);

            % Initialization of alpha's
            alpha_k = h;
            while(dphi(alpha_k) < 0)
                alpha_k = alpha_k + h;
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
            end
            alpha = alpha_k/2;

            % Iterative update
            alpha_kminus1 = 0;
            while(max(abs(alpha - alpha_k) , abs(alpha - alpha_kminus1)) > 1e-5)
                alpha = alpha_k - vpa(dphi(alpha_k)) * (alpha_k - alpha_kminus1)/( vpa(dphi(alpha_k)) - vpa(dphi(alpha_kminus1)) );
                alpha_kminus1 = alpha_k;
                alpha_k = alpha;
                iter = iter + 1;
            end

        case 'BB'
            alpha = solve(diff(phi) == 0);
            alpha = min(alpha( [isAlways(alpha>0)] )); % Take the minimum positive alpha
            iter = iter + 1;

        case 'DIV'
            alpha = 1/(k + 1);
            iter = iter + 1;

        case 'CQ'
            alpha = norm(s, 2)./(transpose(s) * H * s) ;
            iter = iter + 1;
    end
end












