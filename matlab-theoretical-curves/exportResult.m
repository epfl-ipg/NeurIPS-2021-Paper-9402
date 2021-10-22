function exportResult(PATH, res)

fid = fopen(PATH, 'w');

fwrite(fid,double(length(res.GX1)),"int");

fwrite(fid,double(res.phi),"double");
fwrite(fid,double(res.psi),"double");
fwrite(fid,double(res.mu),"double");
fwrite(fid,double(res.nu),"double");

fwrite(fid,double(res.GX1_0),"double");
fwrite(fid,double(res.GX3_0),"double");
fwrite(fid,double(res.HX4_0),"double");
fwrite(fid,double(res.TX1_0),"double");

fwrite(fid,double(imag(res.GX1)), "double");
fwrite(fid,double(imag(res.GX3)), "double");
fwrite(fid,double(imag(res.HX4)), "double");
fwrite(fid,double(imag(res.TX1)), "double");

fwrite(fid,double(res.Q1_00),"double");
fwrite(fid,double(res.Q2_00),"double");
fwrite(fid,double(res.Q4_00),"double");

fwrite(fid,double(res.Q1dirac),"double");
fwrite(fid,double(res.Q2dirac),"double");
fwrite(fid,double(res.Q4dirac),"double");

fwrite(fid,double(res.Q1),"double");
fwrite(fid,double(res.Q2),"double");
fwrite(fid,double(res.Q4),"double");

fwrite(fid,double(res.Q1line),"double");
fwrite(fid,double(res.Q2line),"double");
fwrite(fid,double(res.Q4line),"double");

fwrite(fid,double(res.Xspace),"double");

fclose(fid);




end
