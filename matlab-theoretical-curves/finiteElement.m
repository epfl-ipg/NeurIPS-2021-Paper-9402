function res = finiteElement(phi, psi, mu, nu)
    

c = phi/psi;

%% dirac in 0

[GX1_0, GX3_0, HX4_0, TX1_0] = diracX(phi, psi, mu, nu);

%% Interval 
[a,b] = spectrumInterval2(phi, psi, mu, nu);
[Xspace,RES] = adaptint2(phi, psi, mu, nu, a, b, 0.01); 

err = abs(GX1_0+trapz(Xspace, imag(RES/pi))-1)
% Workaround for wrong err
if err>0.01
    nullVals = find(RES==0)
    if length(nullVals)>0
        fprintf("(%.3f,%.3f) null vals detected\n", phi, psi);

        [aa, ~] = dichotomy(phi, psi, mu, nu, Xspace(nullVals(1)-1), Xspace(nullVals(1)), 20);
        [~, bb] = dichotomy(phi, psi, mu, nu, Xspace(nullVals(end)-1), Xspace(nullVals(end)), 20);

        [XspaceA,RESA] = adaptint3(phi, psi, mu, nu, a, aa, 0.01, 15);
        [XspaceB,RESB] = adaptint3(phi, psi, mu, nu, bb, b, 0.01, 10);

        Xspace = cat(1, XspaceA, XspaceB);
        RES = cat(1, RESA, RESB);
    else
        fprintf("(%.3f,%.3f)trying lowering precision\n", phi, psi);

        [a, b] = spectrumInterval(phi, psi, mu, nu);
        [Xspace,RES] = adaptint2(phi, psi, mu, nu, a, b, 0.01); 
        
    end
    
    err = abs(GX1_0+trapz(Xspace, imag(RES/pi))-1);
end

N = length(Xspace);

fprintf("(%.3f,%.3f) err: %f (for N=%d points)\n", phi, psi, err, N);


%% else

GX1 = zeros(N,1);
GX3 = zeros(N,1);
HX4 = zeros(N,1);
TX1 = zeros(N,1);


% DIRAC IN (0,0)
delta = 0.0000000001;
x = 0.00+1j*delta; 
y = 0.00+1j*delta;

[gx1, tx1, gx3, hx4] = solveXfast2(phi, psi, mu, nu, x);
[q1, q2, q4] = solveXYdistFast2(phi, psi, mu, nu, x, y, gx1, hx4, tx1, gx1, hx4, tx1);


Q1_00 = q1*delta*delta;
Q2_00 = q2*delta*delta;
Q4_00 = q4*delta*delta;


% DIAGONAL

Q1dirac = zeros(N,1);
Q2dirac = zeros(N,1);
Q4dirac = zeros(N,1);


for i=1:N
    x = Xspace(i)+1j*delta;
    [gx1, tx1, gx3, hx4] = solveXfast2(phi, psi, mu, nu, x);
    [q1, q2, q4] = solveXYdistFast2(phi, psi, mu, nu, x, x, gx1, hx4, tx1, gx1, hx4, tx1);
    Q1dirac(i) = q1*delta;
    Q2dirac(i) = q2*delta;
    Q4dirac(i) = q4*delta;
end


% THE REST
for i=1:N
    [gx1, tx1, gx3, hx4] = solveXfast2(phi, psi, mu, nu, Xspace(i));
    GX1(i) = gx1;
    GX3(i) = gx3;
    HX4(i) = hx4;
    TX1(i) = tx1;
end

fprintf("(%.3f,%.3f) the following should be equal to 1: %f (for N=%d points)\n", phi, psi, trapz(Xspace, imag(GX1)/pi) + GX1_0, N);


Q1 = zeros(N,N);
Q2 = zeros(N,N);
Q4 = zeros(N,N);


for j=2:N
    fprintf("(%.3f,%.3f) progress: %.2f\n", phi, psi, 100*j/N);
    i = 1:(j-1);
    [q1, q2, q4] = solveXYdistFast2(phi, psi, mu, nu, Xspace(i), Xspace(j), GX1(i), HX4(i), TX1(i), GX1(j), HX4(j), TX1(j));
    Q1(i,j) = q1;
    Q2(i,j) = q2;
    Q4(i,j) = q4;
end

Q1 = Q1 + Q1';
Q2 = Q2 + Q2';
Q4 = Q4 + Q4';


y = 0 + delta*1j;
[gy1, ty1, ~, hy4] = solveXfast2(phi, psi, mu, nu, y);
[Q1line, Q2line, Q4line] = solveXYdistFast2(phi, psi, mu, nu, Xspace, y, GX1, HX4, TX1, gy1, hy4, ty1);
Q1line = Q1line*delta;
Q2line = Q2line*delta;
Q4line = Q4line*delta;

%% Export

res = struct();

res.phi = phi;
res.psi = psi;
res.mu = mu;
res.nu = nu;

res.Xspace = Xspace;

res.GX1 = GX1;
res.GX3 = GX3;
res.HX4 = HX4;
res.TX1 = TX1;

res.GX1_0 = GX1_0;
res.GX3_0 = GX3_0;
res.HX4_0 = HX4_0;
res.TX1_0 = TX1_0;

res.Q1 = Q1;
res.Q2 = Q2;
res.Q4 = Q4;

res.Q1dirac = Q1dirac;
res.Q2dirac = Q2dirac;
res.Q4dirac = Q4dirac;

res.Q1_00 = Q1_00;
res.Q2_00 = Q2_00;
res.Q4_00 = Q4_00;

res.Q1line = Q1line;
res.Q2line = Q2line;
res.Q4line = Q4line;

end 






