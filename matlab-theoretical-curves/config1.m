%%%%%%%%%%%%%%%%%%%%%%
% as in config 11
Cinv_SPACE = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 5.0, 10.0, 15.0, 20., 50., 100.];
Cinv_LEN = length(Cinv_SPACE);

mu = 0.5/0.3; 
nu = 1.0; 

phi = 3.;

PATH = "./output/config1";
mkdir(PATH);

PARAMS = [(1:Cinv_LEN)', phi*ones(Cinv_LEN,1), phi*Cinv_SPACE'];


parpool(Cinv_LEN)
parfor i=1:Cinv_LEN
%for i=1:Cinv_LEN
    c_inv=Cinv_SPACE(i);
    psi = c_inv*phi;
    
    fprintf("start calculation (%d)\n", i);
    
    res = finiteElement(phi, psi, mu, nu);
    exportResult(PATH+"/MSE_"+i+"_log.bin", res);
    
    fprintf("stop calculation (%d)\n", i);
    
    if i==1 
        TMP = array2table(PARAMS, 'VariableNames',{'id', 'phi','psi'});
        writetable(TMP, PATH+"/params.csv");
    end
    
end
