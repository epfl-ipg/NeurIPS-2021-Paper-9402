function [I,Q] = adaptint(phi, psi, mu, nu, a,b,Tol)
    fa = f(phi, psi, mu, nu, a); 
    fb = f(phi, psi, mu, nu, b);
    QT = 0.5*(b-a)*(fa+fb);        % Trapezoid rule

    c = (a+b)/2; 
    fc = f(phi, psi, mu, nu, c);
    %plot([c c],[-1,1]*1e-2,'b'); % plot new point

    QS = (b-a)/6*(fa+4*fc+fb);   % Simpson rule
    % for small intervals QS is much closer to I than QT
    % hence we can approximate error QT-I by QT-QS

    %disp([a,b]);
    %if (depth>=8 & (fa==0 | fb==0)) | 
    if (abs(QT-QS)<=Tol) | ((b-a)<0.001)
        I = [a;b];
        Q = [fa;fb];
    else
        [Ia,Qa] = adaptint(phi, psi, mu, nu, a,c,0.5*Tol);
        [Ib,Qb] = adaptint(phi, psi, mu, nu, c,b,0.5*Tol);
        I = cat(1, Ia(1:end-1), Ib);
        Q = cat(1, Qa(1:end-1), Qb);
    end
end


function [gx1] = f(phi, psi, mu, nu, x)
    hx1 = sym("hx1");
    hx4 = sym("hx4");
    gx1 = sym("gx1");
    tx1 = sym("tx1");

    EQx = [nu^2*(gx1.^2).*x+gx1.*(-mu^2*hx4 - nu^2*phi/psi + nu^2 + x) + 1, ...
    -mu^2*(hx1.^2)/psi + nu^2/psi + hx1.*(-mu^2*gx1.*x - mu^2 + mu^2/psi - nu^2/psi), ...
     hx1.*(phi/psi - gx1.*x - 1) - hx4, ...
    -mu^2*gx1.*hx4 - hx1/psi + psi^(-1), ...
    -mu*tx1 - hx1 + 1
    ];

    res = vpasolve(EQx, [gx1 tx1 hx1 hx4], [-1e10+1e-10i 1e10+1e10i; NaN NaN; -1e10+0i 1e10-1e10i; NaN NaN; ]);

    gx1 = res.gx1;

    if isempty(gx1) 
        gx1 = 0.;
    end
end