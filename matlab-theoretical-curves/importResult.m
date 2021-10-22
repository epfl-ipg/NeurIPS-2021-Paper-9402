function res = importResult(PATH)

res = struct();

fid = fopen(PATH, 'r');

N = fread(fid, 1, "int")

res.phi = fread(fid, 1, "double");
res.psi = fread(fid, 1, "double");
res.mu = fread(fid, 1, "double");
res.nu = fread(fid, 1, "double");

res.GX1_0 = fread(fid, 1, "double");
res.GX3_0 = fread(fid, 1, "double");
res.HX4_0 = fread(fid, 1, "double");
res.TX1_0 = fread(fid, 1, "double");

res.GX1 = 1j*fread(fid, [N], "double");
res.GX3 = 1j*fread(fid, [N], "double");
res.HX4 = 1j*fread(fid, [N], "double");
res.TX1 = 1j*fread(fid, [N], "double");

res.Q1_00 = fread(fid, 1, "double");
res.Q2_00 = fread(fid, 1, "double");
res.Q4_00 = fread(fid, 1, "double");

res.Q1dirac = fread(fid, [N], "double");
res.Q2dirac = fread(fid, [N], "double");
res.Q4dirac = fread(fid, [N], "double");

res.Q1 = fread(fid, [N,N], "double");
res.Q2 = fread(fid, [N,N], "double");
res.Q4 = fread(fid, [N,N], "double");

res.Q1line = fread(fid, [N], "double");
res.Q2line = fread(fid, [N], "double");
res.Q4line = fread(fid, [N], "double");

res.Xspace = fread(fid, [N], "double");

fclose(fid);



end
