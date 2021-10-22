% This is an interesting profile which could be used for a config
%warning('off', 'symbolic:numeric:NumericalInstability')

%mu = 1.0; 
%nu = 0.15; 
mu = 0.5;
nu = 0.3;

phi = 5; % new

delta0 = 0.0001; 
r=1.;
s= 0.2; % new


Cinv_SPACE = round(logspace(-2,2,50),2);
%Cinv_SPACE = round(linspace(0.1,2,30),2);
Cinv_LEN = length(Cinv_SPACE);

MSEtrain = zeros(Cinv_LEN, 1);
MSEtest = zeros(Cinv_LEN, 1);

for i=1:Cinv_LEN
    cinv=Cinv_SPACE(i);
    psi = cinv*phi;
    
    fprintf("start calculation (%d)\n", i);
    
    
    c = phi/psi;

    [MSEtrainTMP, MSEtestTMP] = asymptoticMSE(phi, psi, mu, nu, delta0, s);
    MSEtrain(i) = MSEtrainTMP;
    MSEtest(i) = MSEtestTMP;

    
end

semilogx(Cinv_SPACE, MSEtrain); hold on;
semilogx(Cinv_SPACE, MSEtest); hold on;

%ar = array2table([C_SPACE', MSEtrain, MSEtest], 'VariableNames',{'c', 'MSEtrain', 'MSEtest'});
%writetable(ar, "triple_descent_profile.csv");




