function [x,y] = dichotomy(phi, psi, mu, nu, x, y, depth)

    if depth > 0
       
        xy = 0.5*(x+y);
        
        fx = solveXfast2(phi, psi, mu, nu, x);
        fxy = solveXfast2(phi, psi, mu, nu, xy);
        
        if (fx==0 & fxy == 0) | (fx ~= 0 & fxy ~= 0)
            [x,y] = dichotomy(phi, psi, mu, nu, xy, y, depth-1);
        else
            [x,y] = dichotomy(phi, psi, mu, nu, x, xy, depth-1);
        end
    end

end