function alpha_value = alpha(i, index, alpha_const, functionID, grad)
    %--------------
    % Documentation
    %--------------
    % This function computes the value of the step length:
    %
    % 1 - alpha = alpha_const   (Constant value)
    % 
    % 2 - alpha = 1/(index + 1) (Divergence Serie)
    %
    % 3 - alpha = alpha* found by lineSearchMethod where:
    %
    %       - alpha* (strictly convex quadratic function) if functionID = 1
    %
    %       - alpha* (line search method) if functionID diff. from 1
                  
    %
    switch i
        case 1
           alpha_value = alpha_const;

        case 2
           alpha_value = 1/(index + 1);
            
        %--------------------
        % Line search method
        %--------------------
        case 3
            
           % CASE 1 - Strictly convex quadratic function
           if functionID == 1
                
               % Hardcoded A matrix (Not pretty, I know)
               A = [4, -3 ; -3, 4];
                
               % Computing alpha
               alpha_value = norm(grad)/(grad * A * transpose(grad));
              
           % CASE 2 - Line search method
           % https://www.mathworks.com/matlabcentral/answers/506524-line-search-algorithm-help
           else
              disp("A IMPLEMENTER") 
           end
    end
end













