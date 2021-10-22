function [a,b] = spectrumInterval(phi, psi, mu, nu)
    a = 0;
    b = 20;
    
    eps = 0.1;

    fx = solveXfast2(phi, psi, mu, nu, a);
    while fx==0
        a = a + eps;
        fx = solveXfast2(phi, psi, mu, nu, a);
    end
    
    b = a+eps;

    
    fx = solveXfast2(phi, psi, mu, nu, b);
    while fx~=0
        b = b + 1;
        fx = solveXfast2(phi, psi, mu, nu, b);
    end
    
    fx = solveXfast2(phi, psi, mu, nu, b);
    while fx==0
        b = b - eps;
        fx = solveXfast2(phi, psi, mu, nu, b);
    end
    
    [b,tmp] = dichotomy0(phi, psi, mu, nu, b, b+eps, 0);
    [tmp,a] = dichotomy0(phi, psi, mu, nu, a-eps, a, 0);
    
end



function [x,y] = dichotomy0(phi, psi, mu, nu, x, y, depth)

    if depth < 10
       
        xy = 0.5*(x+y);
        
        fx = solveXfast2(phi, psi, mu, nu, x);
        fxy = solveXfast2(phi, psi, mu, nu, xy);
        
        if (fx==0 & fxy == 0) | (fx ~= 0 & fxy ~= 0)
            [x,y] = dichotomy0(phi, psi, mu, nu, xy, y, depth+1);
        else
            [x,y] = dichotomy0(phi, psi, mu, nu, x, xy, depth+1);
        end
    end

end