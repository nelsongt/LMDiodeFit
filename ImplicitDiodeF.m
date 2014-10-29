function fout = ImplicitDiodeF(params,vdata)
% Constants
q = 1;                  % Elementary charge. Unit: e
kb = 8.617 * 10^-5;     % Boltzmann constant. Unit: eV/K
T = 300;                % Temperature. Unit: 
Vth = kb*T/q;           % Thermal voltage. Unit: 

fout = zeros(size(vdata));
i=0;
for V = vdata
    i = i+1;
    %fout(i) = fzero(@(I) (params(1) * exp((V - I*params(3))/(params(2)*Vth)) + ((V - I*params(3))/params(4)) - I), 0.1,optimset('display','off'));    
    %fout(i) = fsolve(@(I) (params(1) * exp((V - I*params(3))/(params(2)*Vth)) + ((V - I*params(3))/params(4)) - I), 1);
    %guess_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + ((V - I*x(3))/x4) - I), 1);    
    fout(i) = fzero(@(I) (params(1) * exp((V - I*params(3))/(params(2)*Vth)) + params(5) * exp((V - I*params(3))/(params(6)*Vth)) + ((V - I*params(3))/params(4)) - I), 1);
    %guess_I(i) = fzero(@(I) (x(1) * exp((V - I*x(3))/(x(2)*Vth)) + x(4) * exp((V - I*x(3))/(x6*Vth)) + ((V - I*x(3))/x4) - I), 1)
end    