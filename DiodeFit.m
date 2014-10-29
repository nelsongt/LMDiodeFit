clear

% Constants
q = 1;                  % Elementary charge. Unit: e
kb = 8.617 * 10^-5;     % Boltzmann constant. Unit: eV/K
T = 300;                % Temperature. Unit: T
Vth = kb*T/q;           % Thermal voltage. Unit: V


% Plot boundaries & accuracy (used for plotting the model)
Vmin = 0;
Vmax = 1.4;
Vstep = 0.01;


% Load data to be fitted - comma delimited with single header row
datafile = 'Sample_IV.csv';
M = csvread(datafile,1,0);
data_V = M(:,1).';
data_I = M(:,2).';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for known parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Known "correct" parameters - from NRL multibands
J01 = 1.268 * 10^-5;
n1 = 1.62;
Rs = 2.14*10^-5;
Rsh = 1 * 10^10;

% Solve for known currents
i=0;
for V = data_V
    i = i+1;
    known_I(i) = fzero(@(I) (J01 * exp((V - I*Rs)/(n1*Vth)) + ((V - I*Rs)/Rsh) - I), 1);
end

% plot for visual inspection
%semilogy(data_V,known_I,'r')
%hold on
%semilogy(M(:,1),M(:,2),'b')

% Calculate the sum of square difference for comparison to my fit
i=0;
for V = data_V
    i = i+1;
    dev_I(i) = (known_I(i) - data_I(i))^2;
end

sq_sum_known = sum(dev_I);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% The following is my attempt at a solution using a recursive LMF
% least-squares fit


% Initial guesses
x(1) = 1 * 10^-5;     % Dark current
x(2) = 1;           % Ideality factor
x(3) = 2 * 10^-5;     % Series resistance
x(4) = 100000;      % Shunt resistance
%x4 = 100;
x(5) = 1 * 10^-5;     % Diode 2 dark current
x(6) = 2;           % Diode 2 ideality
%x6 = 2.0;


% Main loop -
% Loop finds the calculated current values based on the initial guesses of
% the parameters. Next, based on these (wrong) current values, the
% (incomplete) parameters are calculated using a
% Levenberg-Marquardt-Fletcher algorithm. Using these (more-correct)
% parameters, the values for current are again calculated, which should be
% converging towards the correct values. Once the correct values for
% current and parameters have both converged, the fit will be complete.
for n = 1:20
    
    % Solve for guessed currents - updates guessed currents every main loop
    % iteration
    i=0;
    for V = data_V
        i = i+1;
        %guess_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + ((V - I*x(3))/x(4)) - I), 1);
        %guess_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + ((V - I*x(3))/x4) - I), 1);
        guess_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + x(5) * exp((V - I*x(3))/(x(6)*Vth)) + ((V - I*x(3))/x(4)) - I), 1);
        %guess_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + x(4) * exp((V - I*x(3))/(x6*Vth)) + ((V - I*x(3))/x4) - I), 1);
    end

    % Logarithmic form of residual is required to properly weigh the errors
    % (else fit will neglect low x- and y-values)
    % This "linearizes" the data if there is negligible series resistance
    % More thorough weighting may provide a better fit for high series resistance. 
    %res = @(c) log(c(1)) + data_V.'/(c(2)*Vth) - guess_I.'*c(3)/(c(2)*Vth) - log(data_I.' - (data_V.' - guess_I.'*c(3))/c(4));
    %res = @(c) log(c(1)) + data_V.'/(c(2)*Vth) - guess_I.'*c(3)/(c(2)*Vth) - log(data_I.' - (data_V.' - guess_I.'*c(3))/x4);
    res = @(c) [log(c(1) * exp((data_V.' - guess_I.'*c(3))/(c(2)*Vth)) + c(5) * exp((data_V.' - guess_I.'*c(3))/(c(6)*Vth))) - log(data_I.' - (data_V.' - guess_I.'*c(3))/c(4))
];
    %res = @(c) log(c(1) * exp((data_V.' - guess_I.'*c(3))/(c(2)*Vth)) + c(4) * exp((data_V.' - guess_I.'*c(3))/(x6*Vth))) - log(data_I.' - (data_V.' - guess_I.'*c(3))/x4);

    % LMF nonlinear least squares fit to the data with the guessed currents
    % for this loop iteration.  This function is itself iterative (nested).
    % MaxIter will improve accuracy at cost of compute time (default is 400)
    x = LMFnlsq(res,x,'Display',0,'MaxIter',400);
    
end

% Calculate the sum of square difference for comparison to NRL fit
i=0;
for V = data_V
    i = i+1;
    dev_I(i) = (guess_I(i) - data_I(i))^2;
end

sq_sum_new = sum(dev_I);


% Solve for final parameters
V1 = Vmin:Vstep:Vmax;
i=0;
for V = V1
    i = i+1;
    %final_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + ((V - I*x(3))/x(4)) - I), 1);
    %final_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + ((V - I*x(3))/x4) - I), 1);
    final_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + x(5) * exp((V - I*x(3))/(x(6)*Vth)) + ((V - I*x(3))/x(4)) - I), 1);
    %final_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + x(4) * exp((V - I*x(3))/(x6*Vth)) + ((V - I*x(3))/x4) - I), 1);
end

% plot for visual inspection
semilogy(V1,final_I,'r')
hold on
semilogy(M(:,1),M(:,2),'b')
hleg1 = legend('fit','data');
