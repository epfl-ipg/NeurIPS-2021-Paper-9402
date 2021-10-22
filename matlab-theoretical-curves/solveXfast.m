function [gx1, hx1, hx4, tx1] = solveXfast(phi, psi, mu, nu, x)
    % old version, see solveXfast2.m
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

    %res = vpasolve(EQx, [gx1 tx1 hx1 hx4]);

    if imag(x)==0
        if real(x)<0
            res = vpasolve(EQx, [gx1 tx1 hx1 hx4], [0. 1e10; NaN NaN; NaN NaN; NaN NaN; ]);
        else
            res = vpasolve(EQx, [gx1 tx1 hx1 hx4], [-1e10+1e-10i 1e10+1e10i; NaN NaN; -1e10+0i 1e10-1e10i; NaN NaN; ]);
        end
    else
        res = vpasolve(EQx, [gx1 tx1 hx1 hx4], [-1e10+0i 1e10+1e10i; NaN NaN; -1e10+0i 1e10-1e10i; NaN NaN; ]);
    end

    gx1 = res.gx1;
    hx1 = res.hx1;
    hx4 = res.hx4;
    tx1 = res.tx1;

    if isempty(gx1) 
        gx1 = 0.;
        hx1 = 0.;
        hx4 = 0.;
        tx1 = 0.;
    end

end