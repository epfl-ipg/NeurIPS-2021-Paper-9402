%  epoch-wise descent in the overparameterized regime

res = importResult("./output/config2/MSE_14_log.bin")


PATH = "./output/config2_lambda_graphs";
mkdir(PATH);

LAMBDA_SPACE = logspace(-8,1.5,50);
i=1;
Tsteps = 100;

MSEtrain = zeros(length(LAMBDA_SPACE), Tsteps);
MSEtest = zeros(length(LAMBDA_SPACE), Tsteps);

for lambda=LAMBDA_SPACE

    mse = timeEvolution(Tsteps, 2.0, 0.5, lambda, res); 
    %semilogx(mse.Tspace, mse.MSEtest); hold on;
    
    %semilogx(mse.Tspace, mse.MSEtrain); hold on;

    %writetable(mse, PATH+"/MSE_"+i+"_log.csv");
    
    MSEtest(i,:) = mse.MSEtest;
    MSEtrain(i,:) = mse.MSEtrain;
    i=i+1;
end

% TMP = array2table(LAMBDA_SPACE', 'VariableNames',{'lambda'});
% writetable(TMP, PATH+"/Lambdaspace.csv");

Tspace = mse.Tspace;

[X,Y] = meshgrid((Tspace), (LAMBDA_SPACE));

mmin = min( min(MSEtrain,[],"all"), min(MSEtest,[],"all"));
%mmax = max( max(MSEtrain,[],"all"), max(MSEtest,[],"all"));
mmax =  max(MSEtest,[],"all");
mmax = 2.
grad = multigradient('preset', 'div.cb.spectral.10');

tiledlayout(1,2)
    
ax1 = nexttile;
pcolor(X,Y,min(MSEtrain,mmax));
hold on;
shading interp;
contour(X,Y,min(MSEtrain,mmax),20,':','LineWidth',0.8,'LineColor','k')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('time (t)') 
ylabel('\lambda') 
title("Train")
caxis([mmin mmax])
colormap(gca, grad);

ax2 = nexttile;
pcolor(X,Y,MSEtest);
hold on;
shading interp;
contour(X,Y,MSEtest,40,':','LineWidth',0.8,'LineColor','k')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlabel('time (t)') 
ylabel('\lambda') 
title("Test")
caxis([mmin mmax])
colormap(gca, grad);

h = colorbar;
h.Title.String = "MSE";
%set(gcf,'Position',[0 0 2*600 0.4*800])
set(gcf,'Position',[0 0 800 0.3*800])





%% R
% pause;
% 
% r_SPACE = linspace(0.,2.,50);
% i=1;
% MSE = zeros(length(r_SPACE), Tsteps);
% for r=r_SPACE
%     mse = timeEvolution(Tsteps, r, 0.4, 0.0001, res); 
%     MSE(i,:) = mse.MSEtest;
%     i=i+1;
% end
% 
% 
% Tspace = mse.Tspace;
% 
% [X,Y] = meshgrid((Tspace), (r_SPACE));
% 
% pcolor(X,Y,MSE);
% hold on;
% shading interp;
% contour(X,Y,MSE,20,':','LineWidth',0.8,'LineColor','k')
% %set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlabel('time (t)') 
% title("Test")
% %caxis([mmin mmax])
% 
% h = colorbar;
% h.Title.String = "MSE";





