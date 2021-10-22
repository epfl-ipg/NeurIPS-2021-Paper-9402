function [gx1, tx1, gx3, hx4] = solveXfast2(phi, psi, mu, nu, x)
    % new version
    tx1 = sym("tx1");
    hx4 = sym("hx4");
    gx1 = sym("gx1");
    gx3 = sym("gx3");

    EQx = [
    gx1.*(-mu.^2.*hx4 - (phi/psi).*gx3.*nu.^2 + x) + 1, ...
    gx3.*(-mu.^2.*phi.*gx1.*hx4 + (phi/psi) ) - hx4  , ...
    (phi/psi)*(1 - gx3) - gx1.*x - 1, ...
    mu.*psi*gx1.*hx4 - tx1 ...
    %mu.*phi.*gx1.*gx3.*(phi.*gx1.*gx3.*nu.^2 + phi.*gx3 - phi + 1) - tx1 ...
    ];
     
    %res = vpasolve(EQx, [gx1 tx1 gx3 hx4], [0. 1e10; NaN NaN; 0. 1e10; NaN NaN; ]);

    if imag(x)==0
        if real(x)<0
            res = vpasolve(EQx, [gx1 tx1 gx3 hx4], [0. 1e10; NaN NaN; 0. 1e10; NaN NaN; ]);
        else
            res = vpasolve(EQx, [gx1 tx1 gx3 hx4], [-1e10+1e-10i 1e10+1e10i; NaN NaN; -1e10+0i 1e10-1e10i; NaN NaN; ]);
        end
    else
        res = vpasolve(EQx, [gx1 tx1 gx3 hx4], [-1e10+0i 1e10+1e10i; NaN NaN; -1e10+0i 1e10-1e10i; NaN NaN; ]);
        %res = vpasolve(EQx, [gx1 tx1 gx3 hx4], [-1e10+0i 1e10+1e10i; NaN NaN; NaN NaN; NaN NaN; ]);
    end

    if length(res.gx1)>1
        fprintf("error! res length greater than 0 (%f %f %f %f %f)\n", phi, psi, mu, nu, x);
        disp(res);
    end
     
    gx1 = res.gx1;
    gx3 = res.gx3;
    hx4 = res.hx4;
    tx1 = res.tx1;

    if isempty(gx1) 
        gx1 = 0.;
        gx3 = 0.;
        hx4 = 0.;
        tx1 = 0.;
    end

end