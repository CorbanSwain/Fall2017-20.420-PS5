function swain_corban_huang_mod
%% Introduction
% This function implements a toy example for the curve fitting procedure
% used in Huang et al. for an implementation in 20.420. The example will
% not be identical to what Huang et al. did, but is meant to demonstrate
% how nlinfit (and ODE solvers) fundamentally work and interact with the
% functions you/we write to make the coding workflow make more sense.

%% "Experimental Data" Generation
% We're going to generate "experimental data" according to a recursive
% formula. Each new data point will be some formula applied to a previous
% one. In this way, the data generation process is ITERATIVE (like an ODE solver!).

n = 5;  % Number of cycles to complete
x = zeros(6, 1);
x(1) = 5; % x will store our six data points (five iterations and x0)

% parameters
a = 2; % scaling factor
b = 0.5; % exponent
c = 10; % constant offset
    % These will be the parameters we'll later fit. See how I used them to
    % generate the "data" in the function "Generate_Data" below.
params = [a, b, c];
    % package them all up
    
% Iteratively generate data:
for step = 1:n
    % Call the helper function "Generate_Next_Point" which applies the
    % recursive formula. Generate_Next_Point is like ODE_fun.
    x(step+1) = Generate_Next_Point(params, x(step));
end

% Visualize data
figure(1)
subplot 121
plot(0:n, x, '-x', 'LineWidth', 3)
xlabel('Iteration Number')
ylabel('Data Value')
title('Raw Data')

%% Curve fitting an iterative process
% Now assume we have this data, x. We know the general formula describing
% the data (x_n = (a*x_{n-1})^b+c). What we don't know are the values for
% a, b, and c. We have to use nlinfit to fit these parameters.

global iterated_parameters
iterated_parameters = ones(1, 3);
    % This is an example of using a global variable. It may be useful with
    % the Huang experimental data, though it is not necessary to use global
    % variables to pass data (such as the Cf18 data. You can pass that
    % using addition inputs to the functions if you prefer). The purpose of
    % this variable is to track how the parameters change during curve
    % fitting to see how nlinfit interacts with an iterative data
    % generation function.

% Use nlinfit to fit parameters. This is where the first trick comes in (in
% my opinion). How do we set up the model_fun to generate iterative data?
% We can't give nlinfit "Generate_Next_Point" as a model function, because
% it needs to compare the real data to more than one fitted point. The
% model function MUST take the parameters and independent variable (as well
% as anything else we want to pass it), implement the relationship between
% input and output that we code, and return predicted values in a vector
% the same size as our input data.
[p_pred] = nlinfit(0:5', x, @(A, B)model_fun(A, B, x(1)), ones(1, 3));
    % The "rules" of nlinfit still apply. A detailed (hopefully)
    % description:
    
    % The first input is our independent variable. In this case, the
    % independent variable is iteration number.
    
    % The second input is the "measured" data. This is that data nlinfit
    % will try to match by changing parameters. For us, this is the vector
    % x.
    
    % The third input is the model function. The model function also has
    % rules that it must follow. The first input to the model function MUST
    % be the parameters, and the second input to the model function MUST be
    % the independent variable. You can see below in the model_fun function
    % that this is the case. Now, I need to pass in an extra value to
    % model_fun. Our data is iterative, and needs (1) the parameters, (2)
    % the number of iterations to complete, and (3) the initial x value to
    % start from. nlinfit knows what x(1) is, but we also need to give that
    % value to model_fun. I did that using the syntax above. The @ before
    % model_fun tells nlinfit that A and B, the first two inputs to
    % model_fun, are the two inputs it's usually expecting (parameters and
    % independent variable. nlinfit will take this information and stick
    % the parameter guesses in the first input and the iteration vector in
    % the second input. The parenthetical stuff after model_fun tells
    % nlinfit that it has to give more info to model_fun. That extra info
    % is x(1), the starting point for the iterations.
    
    % The fourth input MUST be the initial guess for the parameters. My
    % favorite initial guess is a vector of ones. If that doesn't work, I
    % usually shed a tear for the failure of laziness, then go and pick new
    % guesses.
    
% Given estimated parameter values, generate and visualize predicted x
% values. Note that we can use model_fun for this data generation because
% that's exactly what it's designed to do!
% Generate predicted data
[x_pred] = model_fun(p_pred, 0:5', x(1));

figure(1)
subplot 122
plot(0:5, x_pred, '-x', 'LineWidth', 3)
xlabel('Iteration Number')
ylabel('Data Value')
title('Predicted or Fit Data')

% Plot the iterated_parameters to see how nlinfit changes the parameters
% over time.
figure(2)
plot(iterated_parameters, '-', 'LineWidth', 3)
xlabel('Fitting Iteration Number')
ylabel('Parameter Value')
title('Parameter values during curvefitting')
legend('a = 2', 'b = 0.5', 'c = 10', 'location', 'best')

%% Trickier example: when you don't want all of the data points
% This is relevant when you have to iterate over more values than you have
% as input data. For example, if the comparison/experimental data is for t
% = 1, 5, 10, but you have to evaluate for t = 1, 2, 3, 4, 5, 6, 7, 8, 9,
% 10.
iterated_parameters = [];
iterated_parameters = ones(1, 3);
tricky_iteration_vector = [2, 3, 4, 5, 10, 12, 13, 14];
tricky_x = model_fun(params, tricky_iteration_vector, 2);
figure(3)
plot(tricky_iteration_vector, tricky_x, '-x', 'LineWidth', 3)

[p_pred2] = nlinfit(tricky_iteration_vector, tricky_x,...
    @(A, B)model_fun(A, B, tricky_x(1)), ones(1, 3));
[All_Predicted] = model_fun(p_pred2, 2:14, 2);
figure(3)
hold on
plot(2:14, All_Predicted, '-x', 'LineWidth', 3)
xlabel('Iteration Number')
ylabel('Data Value')
title('Tricky Example')
legend('Raw Data', 'Fit Data')

% Plot the iterated_parameters to see how nlinfit changes the parameters
% over time.
figure(4)
plot(iterated_parameters, '-', 'LineWidth', 3)
xlabel('Fitting Iteration Number')
ylabel('Parameter Value')
title('Parameter values during curvefitting for tricky example')
legend('a = 2', 'b = 0.5', 'c = 10', 'location', 'best')

% fprintf('Real Values: %g\n', params)
% fprintf('Predicted Values problem 1: %g\n', p_pred)
% fprintf('Predicted Values tricky problem: %g\n', p_pred2)
end

%% Generate_Next_Point helper function
function [x_n2] = Generate_Next_Point(p, x_n1)
% This function calculates the next data point, x_n2, given the parameters,
% p, and the current data point, x_n1.

a = p(1); %a, scaling factor
b = p(2); %b, exponent
c = p(3); %c, constant offset

x_n2 = (a*x_n1)^b+c;

end

%% model_fun helper function
function [x_final] = model_fun(parameters, iteration_vector, initial_x)

global iterated_parameters
iterated_parameters(end+1, :) = parameters;

% disp(parameters)
    % Uncomment the above line if you want to see how the parameters change
    % during the curve fitting procedure!

iteration_start = iteration_vector(1);
iteration_end = iteration_vector(end);
    % We're implementing a recursive formula, so we must know what number
    % iteration we start on (this MUST correspond to initial_x), and what
    % number iteration we end on. The data we're given won't necessarily be
    % iterations 0, 1, 2, 3, 4, 5!. Maybe it's, 2, 5, 9, 11, 12!. model_fun
    % is the place to take care of that.
    
n_cycles = iteration_end-iteration_start;
    % The total number of iterations we have to do.
Temp_Iteration_Vector = iteration_start:iteration_end;
    
x_pred = zeros(n_cycles+1, 1);
x_pred(1) = initial_x;

% Calculate the predicted x values just like we calculated the real data.
for step = 1:n_cycles
   x_pred(step+1) = Generate_Next_Point(parameters, x_pred(step)); 
end

% Now, our output must only be the x values for the iteration numbers
% specified in iteration_vector. In the Huang implementation, we need to
% control for the independent vector time using interpolation (why? what
% does the ODE solver do to time???). In this case, we need to compare the
% input iteration_vector with the iterations we actualy did,
% Temp_Iteration_Vector, and take only the x values corresponding to
% iteration_vector.

[Keep_Indicies] = ismember(Temp_Iteration_Vector, iteration_vector);
    % See the documentation for ismember if you're interested in how this
    % works. It's not necessary for this homework or the next one. It may
    % be useful in 20.440!
x_final = x_pred(Keep_Indicies);

end