% it shows another the triple descent profile

mu = 0.9;
nu = 0.1;

psi = 2.;

delta0 = 0.0001; 
r=1.;
s=0.8;


%C_SPACE = [0.01, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.5, 2., 3.];
C_SPACE = round(linspace(0.1,2,30),2);
%C_SPACE = round(linspace(0.1,2,100),2);
%C_SPACE = round(logspace(-2,1,40),2);
C_LEN = length(C_SPACE);

MSEtrain = zeros(C_LEN, 1);
MSEtest = zeros(C_LEN, 1);

for i=1:C_LEN
    a=C_SPACE(i);
    
    fprintf("start calculation (%d)\n", i);
    
    %psi = a*phi;
    phi = a*psi;
    
    [MSEtrainTMP, MSEtestTMP] = asymptoticMSE(phi, psi, mu, nu, delta0, s);
    MSEtrain(i) = MSEtrainTMP;
    MSEtest(i) = MSEtestTMP;
    
end

plot(C_SPACE, MSEtrain); hold on;
plot(C_SPACE, MSEtest); hold on;

%ar = array2table([C_SPACE', MSEtrain, MSEtest], 'VariableNames',{'c', 'MSEtrain', 'MSEtest'});
%writetable(ar, "triple_descent_profile.csv");




