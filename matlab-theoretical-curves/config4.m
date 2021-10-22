%%%%%%%%%%%%%%%%%%%%%%
C_SPACE = round(linspace(0.1,2,30),2);
C_LEN = length(C_SPACE);

mu = 0.9;
nu = 0.1;

psi = 2.0;

PATH = "./output/config4";
mkdir(PATH);

PARAMS = [(1:C_LEN)', psi*C_SPACE', psi*ones(C_LEN,1)];

parpool(C_LEN)
parfor i=1:C_LEN
%for i=1:C_LEN
    c=C_SPACE(i);
    phi = c*psi;
    
    fprintf("start calculation (%d)\n", i);
    
    res = finiteElement(phi, psi, mu, nu);
    exportResult(PATH+"/MSE_"+i+"_log.bin", res);
    
    fprintf("stop calculation (%d)\n", i);

    if i==1
        TMP = array2table(PARAMS, 'VariableNames',{'id', 'phi','psi'});
        writetable(TMP, PATH+"/params.csv");
    end
    
end


