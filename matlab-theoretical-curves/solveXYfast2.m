function [q1, q2, q4] = solveXYfast2(phi, psi, mu, nu, X, y, GX1, HX4, TX1, gy1, hy4, ty1)
    q1 = sym("q1_", [length(GX1),1]);
    q2 = sym("q2_", [length(GX1),1]);
    q4 = sym("q4_", [length(GX1),1]);
    q5 = sym("q5_", [length(GX1),1]);
    
    HX1 = 1-mu*TX1;
    hy1 = 1-mu*ty1;

    EQxy = [
        -mu^2*gy1*q2 + mu^2*HX4.*q1 + mu*gy1*TX1 + mu*gy1*ty1 - nu^2*phi*gy1*q4/psi + nu^2*q1.*(-psi*GX1.*X + phi - psi)/psi - gy1 - q1.*X, ... 
        mu*(-psi*GX1.*X + phi - psi).*(-mu*GX1.*q2 + mu*hy4*q1 + GX1*ty1) + phi*hy1*q4/psi - q2, ...
        -mu^2*GX1.*HX1.*q4 + mu^2*q5*(phi - psi*gy1*y - psi)/(phi*psi) - nu^2*GX1.*q4 + nu^2*q1*(phi - psi*gy1*y - psi)/phi - q4, ...
        mu^2*phi*gy1*hy1*GX1.*q4 - mu^2*GX1.*q5.*(phi - psi*GX1.*X - psi)/psi + psi*GX1*gy1*hy1 + hy1.*q1 - q5/psi
    ];

    res = vpasolve(EQxy, [q1; q2; q4; q5]);
    
    q1 = subs(q1, res);
    q2 = subs(q2, res);
    q4 = subs(q4, res);
    %q5 = subs(q5, res);
    
    
end
