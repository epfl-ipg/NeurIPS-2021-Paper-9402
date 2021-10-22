function [mse] = timeEvolution(TSTEPS, r, s, lambda, res)

phi = res.phi;
psi = res.psi;
mu = res.mu;
nu = res.nu;
c = phi/psi;

Xspace = res.Xspace;

GX1 = res.GX1;
GX3 = res.GX3;
HX4 = res.HX4;
TX1 = res.TX1;

GX1_0 = res.GX1_0;
GX3_0 = res.GX3_0;
HX4_0 = res.HX4_0;
TX1_0 = res.TX1_0;

Q1 = res.Q1;
Q2 = res.Q2;
Q4 = res.Q4;

Q1_00 = res.Q1_00;
Q2_00 = res.Q2_00;
Q4_00 = res.Q4_00;

Q1dirac = res.Q1dirac;
Q2dirac = res.Q2dirac;
Q4dirac = res.Q4dirac;

Q1line = res.Q1line;
Q2line = res.Q2line;
Q4line = res.Q4line;

%% Time evolution

GT = zeros(TSTEPS,1);
LT = zeros(TSTEPS,1);
HT = zeros(TSTEPS,1);
UT = zeros(TSTEPS,1);
Tspace = logspace(-3, 10, TSTEPS)';
%Tspace = linspace(0, 80, TSTEPS)';
delta = lambda*c;

for i=1:length(Tspace)
    t = Tspace(i);
    
    %% Train 
    DIAG = (r^2)*(Xspace+delta).*exp(-2*t*(Xspace+delta)) .* imag(GX1)/pi - ...
        ((1-exp(-2*t*(Xspace + delta)))./(Xspace + delta) ) .* imag(-(s^2)*c*GX3-HX4)/pi;
    DIRAC = (r^2)*delta.*exp(-2*t*(delta)) * GX1_0 - ...
        ((1-exp(-2*t*(delta)))./(delta) ) * (-(s^2)*c*GX3_0-HX4_0);
    UT(i) = trapz(Xspace, DIAG) + DIRAC;
    
    %% Test
    
    % gt
    DIAG = ( (1-exp(-t*(Xspace + delta))) ./ (Xspace + delta)) .* imag(TX1)/pi;
    DIRAC = ( (1-exp(-t*(delta))) ./ (delta)) * TX1_0;
    GT(i) = trapz(Xspace, DIAG) + DIRAC;
    
    % lt 
    DIAG = (r^2)*exp(-2*t*(Xspace+delta)) .* imag(GX1)/pi + ...
        (((1-exp(-t*(Xspace + delta)))./(Xspace + delta) ).^2 ) .* imag(-(s^2)*c*GX3-HX4)/pi;
    DIRAC = (r^2)*exp(-2*t*(delta)) * GX1_0 + ...
        (((1-exp(-t*(delta)))./(delta) ).^2 ) * (-(s^2)*c*GX3_0-HX4_0);
    LT(i) = trapz(Xspace, DIAG) + DIRAC;
    
    
    % ht
    
    MESH = exp(-t*(Xspace + Xspace.' + 2*delta)) .* Q1;
    DIAG = exp(-t*(2*Xspace + 2*delta)) .* Q1dirac;
    DIRAC = exp(-t*(2*delta)) * Q1_00;
    LINE = exp(-t*(Xspace+2*delta)) .* Q1line;
    
    H0t = trapz(Xspace, trapz(Xspace, MESH, 2))/(2*pi^2) + trapz(Xspace, DIAG)/pi + DIRAC/2 + trapz(Xspace, LINE)/pi;
    H0t = (r^2)*H0t;
    %fprintf("LINE1 = %f, ", trapz(Xspace, LINE)/pi);
    
    
    MESH = (  ((1-exp(-t*(Xspace + delta)))./(Xspace + delta)) * ((1-exp(-t*(Xspace + delta)))./(Xspace + delta))' ) .* ((s^2)*c*Q4 + Q2);
    DIAG = (( (1-exp(-t*(Xspace + delta)))./(Xspace + delta) ).^2) .* ((s^2)*c*Q4dirac + Q2dirac);
    DIRAC = (( (1-exp(-t*(delta)))./delta )^2 ) * ((s^2)*c*Q4_00 + Q2_00);
    LINE = (  ((1-exp(-t*(Xspace + delta)))./(Xspace + delta)) .* ((1-exp(-t*delta))./delta) ) .* ((s^2)*c*Q4line + Q2line);
    
    Wt = trapz(Xspace, trapz(Xspace, MESH, 2))/(2*pi^2) + trapz(Xspace, DIAG)/pi + DIRAC/2 + trapz(Xspace, LINE)/pi;
    %fprintf("LINE0 = %f\n", trapz(Xspace, LINE)/pi);
    
    
    HT(i) = H0t + Wt;
end

MSEtest = 1 + (s^2) - 2*mu*GT + (mu^2)*HT + (nu^2)*LT; 
MSEtrain = 1 + (s^2) + (1/c)*UT;

mse = array2table([Tspace, MSEtrain, MSEtest], 'VariableNames', {'Tspace', 'MSEtrain', 'MSEtest'});


[MSEtrain_asym, MSEtest_asym] = asymptoticMSE(res.phi, res.psi, res.mu, res.nu, lambda, s);


fprintf("(%.3f,%.3f) final MSEtest: %.2f / %.2f\n", res.phi, res.psi, MSEtest(end), MSEtest_asym);
fprintf("(%.3f,%.3f) final MSEtrain: %.2f / %.2f\n", res.phi, res.psi, MSEtrain(end), MSEtrain_asym);


end
