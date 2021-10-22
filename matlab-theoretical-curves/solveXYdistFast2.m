ffunction [q1,q2,q4] = solveXYdistFast2(phi, psi, mu, nu, x, y, GX1, HX4, TX1, gy1, hy4, ty1)

    [q1a, q2a, q4a] = solveXYfast2(phi, psi, mu, nu, x, y, GX1, HX4, TX1, gy1, hy4, ty1);
    
    y = conj(y);
    gy1 = conj(gy1);
    hy4 = conj(hy4);
    ty1 = conj(ty1);
    
    [q1b, q2b, q4b] = solveXYfast2(phi, psi, mu, nu, x, y, GX1, HX4, TX1, gy1, hy4, ty1);

    q1 = real(q1b - q1a);
    q2 = real(q2b - q2a);
    q4 = real(q4b - q4a);
end