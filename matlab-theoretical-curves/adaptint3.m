function [I,Q] = adaptint3(phi, psi, mu, nu, a,b,Tol,depth)
    fa = solveXfast2(phi, psi, mu, nu, a); 
    fb = solveXfast2(phi, psi, mu, nu, b);
    QT = 0.5*(b-a)*(fa+fb);        % Trapezoid rule

    c = (a+b)/2; 
    fc = solveXfast2(phi, psi, mu, nu, c);

    QS = (b-a)/6*(fa+4*fc+fb);   % Simpson rule
    % for small intervals QS is much closer to I than QT
    % hence we can approximate error QT-I by QT-QS

    fprintf("%d %.20f %.20f\n", depth, a, b);
    %if (depth>=8 & (fa==0 | fb==0)) | 
    if (abs(QT-QS)<=Tol) || depth<1
        I = [a;b];
        Q = [fa;fb];
    else
        [Ia,Qa] = adaptint3(phi, psi, mu, nu, a,c,0.5*Tol,depth-1);
        [Ib,Qb] = adaptint3(phi, psi, mu, nu, c,b,0.5*Tol,depth-1);
        I = cat(1, Ia(1:end-1), Ib);
        Q = cat(1, Qa(1:end-1), Qb);
    end
end
