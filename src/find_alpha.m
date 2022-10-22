function [alpha] = find_alpha(phi, method, h, k, H, s)
    %--------------
    % Documentation
    %--------------
    % find_alpha returns the minimizer alpha of phi(alpha).
    % method == 'NR'  --> Newton-Raphson
    % method == 'S'   --> Secant
    % method == 'D'   --> Dichotomy
    % method == 'BB'  --> Black-box 'solve' function @Mathworks
    % method == 'DIV' --> Use divergent serie ~ 1/(k+1) for alpha_k
    % method == 'CQ'  --> Formula for convex quadratic functions
    % h is the step used for finding bounds for Secant and Dichotomy methods
    % k is the iteration at which find_alpha is called
    % H is the Hessian matrix of f
    % s is the gradient of f at the considered point

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
            end

        case 'BB'
            alpha = solve(diff(phi) == 0);
            % Take the minimum positive alpha
            alpha = min(alpha( [isAlways(alpha>0)] ));

        case 'DIV'
            alpha = 1/(k + 1);

        case 'CQ'
            alpha = norm(s, 2)./(transpose(s) * H * s) ;
    end
end











