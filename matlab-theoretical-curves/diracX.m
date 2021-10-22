function [GX1_0, GX3_0, HX4_0, TX1_0] = diracX(phi, psi, mu, nu)
    delta = 0.0000000001;
    [gx1, tx1, gx3, hx4] = solveXfast2(phi, psi, mu, nu, delta*1j);
    GX1_0 = imag(gx1)*delta;
    GX3_0 = imag(gx3)*delta;
    HX4_0 = imag(hx4)*delta;
    TX1_0 = imag(tx1)*delta;
end
