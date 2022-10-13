%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          MECA0027 Structural and Multidisciplinary Optimization         %
%                     Unconstrained Optimization                          %
%                    University Of Liège, Belgium                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solve the minimization problem
%
%             min f(x,y)
%
% using the following optimization methods
%
% 1) Steepest descent
% 2) Conjugate gradients with Fletcher-Reeves update rule
% 3) BFGS Quasi-Newton
%
clearvars; close all; clc

%% Parameters --------------------------------------------------------------------------
functionID = 1;

xinit   = [0 0];         % initial point
MaxIter = 100;           % Maximum number of iterations
Epsilon = 1e-5;          % Tolerance for the stop criteria

%% Initialization

disp('Which optimization method do you want to use? Press:')
disp('     1 for Steepest descent method')
disp('     2 for Conjugate gradients method with Fletcher-Reeves update rule')
disp('     3 for BFGS Quasi-Newton method')
prompt = '';
method = input(prompt);

n=2;                        % Dimension of the problem
xinit=reshape(xinit,2,1);   % To be sure that it's a column vector
x=zeros(n,MaxIter);         % Initialization of vector x
x(:,1)=xinit;               % Put xinit in vector x.


% Symbolic function and variables
syms x1 x2;
X = [x1 x2];
f = getObjFVal(X, functionID);

% Gradient computation
df_dx1 = diff(f,'x1')
df_dx2 = diff(f,'x2')
grad_f = [df_dx1 ; df_dx2];

% Hessian computation
d2f_dx1dx1 = diff(df_dx1, 'x1')
d2f_dx1dx2 = diff(df_dx1, 'x2')
d2f_dx2dx1 = diff(df_dx2, 'x1')
d2f_dx2dx2 = diff(df_dx2, 'x2')
H = [d2f_dx1dx1 d2f_dx1dx2;
     d2f_dx2dx1 d2f_dx2dx2];

%% Methods -----------------------------------------------------------------------------

switch method
    case 1
        disp('You chose the steepest descent method.')
        
        for i=1:MaxIter
            %%%% ---------------------------------------------------------------------------
            %%%% ADD YOUR CODE
            %%%%
        end
        
        x=x(:,1:i); %Remove the zero elements due to the initialization step
        
    case 2
        disp('You chose the conjugate gradients method with Fletcher-Reeves update rule.')
        
        for i=1:MaxIter
            %%%% ---------------------------------------------------------------------------
            %%%% ADD YOUR CODE
            %%%%
        end
        
        x=x(:,1:i); %Remove the zero elements due to the initialization step
        
    case 3
        disp('You chose the BFGS Quasi-Newton method.')
        
        for i=1:MaxIter
            %%%% ---------------------------------------------------------------------------
            %%%% ADD YOUR CODE
            %%%%
        end
        
        x=x(:,1:i); %Remove the zero elements due to the initialization step
        
end

%% Plot the function with the optimization path and the results

fprintf('The optimal point is: x = %4.3f, y = %4.3f.\n', x(1,end), x(2,end))
fprintf('The objective function value is: %4.3f.\n',getObjFVal(x(:,end),functionID))
plotOptimizationPath(x,functionID)

