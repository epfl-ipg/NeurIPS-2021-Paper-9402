function [a,b] = spectrumInterval2(phi, psi, mu, nu) 
    eps_a = 0.01;
    a = 0;
    
    fx = solveXfast2(phi, psi, mu, nu, a);
    while fx==0
        a = a + eps_a;
        fx = solveXfast2(phi, psi, mu, nu, a);
    end
    
    b = max(a+10,50); %a+eps;
    eps_b = 0.1;
    
    fx = solveXfast2(phi, psi, mu, nu, b);
    while fx~=0
        b = b + 5;
        fx = solveXfast2(phi, psi, mu, nu, b);
    end
    
    fx = solveXfast2(phi, psi, mu, nu, b);
    while fx==0
        b = b - eps_b;
        fx = solveXfast2(phi, psi, mu, nu, b);
    end
      
    [b,tmp] = dichotomy(phi, psi, mu, nu, b, b+eps_b, 20);
    [tmp,a] = dichotomy(phi, psi, mu, nu, a-eps_a, a, 20);
    
end

