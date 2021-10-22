res = importResult("./output/config2/MSE_14_log.bin")

LAMBDA_SPACE = logspace(-8,1.5,50);
i=1;
Tsteps = 100;

MSEtest = zeros(length(LAMBDA_SPACE), Tsteps);

for lambda=LAMBDA_SPACE

    mse = timeEvolution(Tsteps, 2.0, 0.5, lambda, res); 
    
    MSEtest(i,:) = mse.MSEtest;
    i=i+1;
end

Tspace = mse.Tspace;

[X,Y] = meshgrid((Tspace), (LAMBDA_SPACE));


grad = multigradient('preset', 'div.cb.spectral.10');

mmin = min(MSEtest,[],"all");
mmax = max(MSEtest,[],"all");

%% Graph    
tiledlayout(1,2)

ax1 = nexttile;
surfc(X,Y,MSEtest);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
colormap(gca, grad);
xlabel('time (t)') 
ylabel('\lambda') 
zlabel('MSE Test')
caxis([mmin mmax])
zlim([0. mmax])
colormap(gca, grad);

h = colorbar;
h.Title.String = "MSE";




%% second graph

PATH = "./output/config4";
lambda = 0.0001; r=1.; s= 0.8;


params = readtable(PATH+"/params.csv");
np = length(params.id);

if REVERSE_PSI_PHI
    Xspace = params.phi./params.psi;
else
    Xspace = params.psi./params.phi;
end

Tsteps = 100;

MSEtest = zeros(np, Tsteps);
MSEtrain = zeros(np, Tsteps);

for i=1:np
    res = importResult(PATH+"/MSE_"+i+"_log.bin");
    cinv = params.psi(i)/params.phi(i);
    if CINV
        mse = timeEvolution(100, r, s, lambda*cinv, res);
    else
        mse = timeEvolution(100, r, s, lambda, res);
    end
    
    MSEtest(i,:) = mse.MSEtest;
    MSEtrain(i,:) = mse.MSEtrain;
end

Tspace = mse.Tspace;

[X,Y] = meshgrid((Tspace), (Xspace));

%mmin = min( min(MSEtrain,[],"all"), min(MSEtest,[],"all"));
%mmax = max( max(MSEtrain,[],"all"), max(MSEtest,[],"all"));

mmin = min(MSEtest,[],"all");
mmax = max(MSEtest,[],"all");


%% plot

ax2 = nexttile;
surfc(X,Y,MSEtest);
set(gca, 'XScale', 'log')
if LOGY_SCALE
    set(gca, 'YScale', 'log')
end
xlabel('time (t)') 
ylabel("\phi/\psi") 
zlabel('MSE Test')
caxis([mmin mmax])
zlim([0. mmax])
colormap(gca, grad);

h = colorbar;
h.Title.String = "MSE";

set(gcf,'Position',[0 0 800 0.4*800])
