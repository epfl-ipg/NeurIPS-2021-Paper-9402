%%%%%%%%%%%%%%%%%%%%%%
Cinv_SPACE = round(logspace(-1,2,31),2);
Cinv_LEN = length(Cinv_SPACE);

mu = 0.5; 
nu = 0.3; 

phi = 0.5;

PATH = "./output/config6";
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
