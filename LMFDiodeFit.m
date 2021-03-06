clear

% Constants
q = 1;                  % Elementary charge. Unit: e
kb = 8.617 * 10^-5;     % Boltzmann constant. Unit: eV/K
T = 300;                % Temperature. Unit: T
Vth = kb*T/q;           % Thermal voltage. Unit: V


% Plot boundaries & accuracy (used for plotting the model)
Vmin = -1.0;
Vmax = 1.4;
Vstep = 0.01;


% Load data to be fitted - comma delimited with single header row
datafile = 'Sample_IV.csv';
M = csvread(datafile,1,0);
data_V = M(:,1).';
data_I = M(:,2).';



% The following is my attempt at a solution using a recursive LMF
% least-squares fit


% Initial guesses
x(1) = 0.0003667;     % Dark current (mA/cm^2)
x(2) = 1.0722;           % Ideality factor
x(3) = 0.0014;     % Series resistance  (kOhm*cm^2)
x(4) = 0.019      % Shunt resistance  (kOhm*cm^2)
%x4 = 100;
x(5) = 0.0933;     % Diode 2 dark current
x(6) = 1.8025;           % Diode 2 ideality
%x6 = 2.0;


% Main -
% Loop finds the calculated current values based on the initial guesses of
% the parameters. Next, based on these (wrong) current values, the
% (incomplete) parameters are calculated using a
% Levenberg-Marquardt-Fletcher algorithm. Using these (more-correct)
% parameters, the values for current are again calculated, which should be
% converging towards the correct values. Once the correct values for
% current and parameters have both converged, the fit will be complete.


    % Logarithmic form of residual is required to properly weigh the errors
    % (else fit will neglect low x- and y-values)
    % This linearizes the data if there is negligible series resistance
    % More thorough weighting may provide a better fit for high series resistance. 
    %res = @(c) log(c(1)) + data_V.'/(c(2)*Vth) - guess_I.'*c(3)/(c(2)*Vth) - log(data_I.' - (data_V.' - guess_I.'*c(3))/c(4));
    %res = @(c) log(c(1)) + data_V.'/(c(2)*Vth) - guess_I.'*c(3)/(c(2)*Vth) - log(data_I.' - (data_V.' - guess_I.'*c(3))/x4);
    res = @(c) [log(abs(ImplicitDiodeF(c,data_V))).' - log(abs(data_I).')
        (c(1)<0)*c(1)*100000000
        (c(2)<0)*c(2)*100
        (c(3)<0)*c(3)*10000000
        (c(4)<0)*c(5)*100000000
        (c(5)<0)*c(5)*1000000
        (c(6)<0)*c(6)*100]
    %res = @(c) log(c(1) * exp((data_V.' - guess_I.'*c(3))/(c(2)*Vth)) + c(4) * exp((data_V.' - guess_I.'*c(3))/(x6*Vth))) - log(data_I.' - (data_V.' - guess_I.'*c(3))/x4);

    reslife1 = log(abs(ImplicitDiodeF(x,data_V))).' - log(abs(data_I).');
    
% LMF nonlinear least squares fit to the data with the guessed currents
% for this loop iteration.  This function is itself iterative (nested).
% MaxIter will improve accuracy at cost of compute time (default is 400)
    
%x = LMFnlsq2(@(x) (log(ImplicitDiodeF(x,data_V)).' - log(data_I.')),x,'Display',0,'MaxIter',10);
%x = LMFnlsq2(res,x,'Display',[1,0],'MaxIter',100);
%x = lsqcurvefit(@(x,dava_V) ImplicitDiodeF(x,data_V),x,data_V,data_I);

reslife2 = log(abs(ImplicitDiodeF(x,data_V))).' - log(abs(data_I).');

% Calculate the sum of square difference for comparison to NRL fit
i=0;
guess_I = ImplicitDiodeF(x,data_V);
for V = data_V
    i = i+1;
    dev_I(i) = (guess_I(i) - data_I(i))^2;
end

sq_sum_new = sum(dev_I);


% Solve for final parameters
V1 = Vmin:Vstep:Vmax;
final_I = ImplicitDiodeF(x,V1);


% plot for visual inspection
semilogy(V1,abs(final_I),'r')
hold on
semilogy(M(:,1),abs(M(:,2)),'b')
hleg1 = legend('fit','data');
