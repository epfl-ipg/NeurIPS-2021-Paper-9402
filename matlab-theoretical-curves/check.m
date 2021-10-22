res = importResult("./output/config7/MSE_11_log.bin");
Xspace = res.Xspace;

d = 1000;
N = round(res.psi*d,0);
n = round(res.phi*d,0);

X = randn(n,d);
Theta = randn(N,d);
XThetaT = X*Theta'/(d^0.5);
Z = max(XThetaT,0)-1/(2*pi)^0.5;
ZtX = Z'*X;
ThetaThetaT = Theta*Theta';

Q1 = trace(ThetaThetaT/d)/N;
Q2 = trace(ZtX'*ThetaThetaT*ZtX/d/N/N)/d;
Q4 = trace(Z*ThetaThetaT*Z'/d/N)/n;


MESH = res.Q1;
DIAG = res.Q1dirac;
DIRAC = res.Q1_00;
A = trapz(Xspace, trapz(Xspace, MESH, 2))/(2*pi^2);
B = trapz(Xspace, DIAG)/pi;
C = DIRAC/2;
Q = A+B+C;
fprintf("Q4 = %f | %f (err=%f) (%f %f %f)\n", Q, Q1, abs(Q-Q1)/Q1, A, B, C);

MESH = res.Q2;
DIAG = res.Q2dirac;
DIRAC = res.Q2_00;
A = trapz(Xspace, trapz(Xspace, MESH, 2))/(2*pi^2);
B = trapz(Xspace, DIAG)/pi;
C = DIRAC/2;
Q = A+B+C;
fprintf("Q4 = %f | %f (err=%f) (%f %f %f)\n", Q, Q2, abs(Q-Q2)/Q2, A, B, C);

MESH = res.Q4;
DIAG = res.Q4dirac;
DIRAC = res.Q4_00;
A = trapz(Xspace, trapz(Xspace, MESH, 2))/(2*pi^2);
B = trapz(Xspace, DIAG)/pi;
C = DIRAC/2;
Q = A+B+C;
fprintf("Q4 = %f | %f (err=%f) (%f %f %f)\n", Q, Q4, abs(Q-Q4)/Q4, A, B, C);


