function [MSEtrain, MSEtest] = asymptoticMSE(phi, psi, mu, nu, lambda, s)

c = phi/psi;

x = -lambda*c;
deltax = 1e-10;

% derivative computation
[gx1, tx1, gx3, hx4] = solveXfast2(phi, psi, mu, nu, x-deltax);
V0 = (s^2)*(phi/psi)*(1-gx3)+(c-hx4);

% repeat again
[gx1, tx1, gx3, hx4] = solveXfast2(phi, psi, mu, nu, x);
V1 = (s^2)*(phi/psi)*(1-gx3)+(c-hx4);

% now we have V
V = (V1-V0)/deltax;

[q1,q2,q4] = solveXYfast2(phi, psi, mu, nu, x, x, gx1, hx4, tx1, gx1, hx4, tx1);

K = tx1;
W = (s^2)*(phi/psi)*q4 + q2;

MSEtrain = 1 + (s^2) - (1/c)*V1;
MSEtest = 1 + (s^2) -2*mu*K + (mu^2)*W + (nu^2)*V;


end