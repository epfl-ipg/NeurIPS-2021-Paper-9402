%%%%%%%%%%%%%%%%%%%%%%
% to compare with experimental results

mu = 0.5; 
nu = round(0.5*sqrt(1 - 2/pi), 4); 

psi = 1.8;
phi = 1.4;

PATH = "./output/config0";
mkdir(PATH);


res = finiteElement(phi, psi, mu, nu);
exportResult(PATH+"/MSE_"+0+"_log.bin", res);
res = importResult(PATH+"/MSE_"+0+"_log.bin")

r = 1.;
s = 0.;
lambda = 0.01;

mse = timeEvolution(100, r, s, lambda, res);

writetable(mse, PATH+"/mse1.csv");