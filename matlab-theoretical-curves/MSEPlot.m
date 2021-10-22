CONFIG = 4;
PLOT_3D = false;
REVERSE_PSI_PHI = false;
LOGY_SCALE = true;
CINV = false;

if CONFIG==2
    PATH = "./output/config2";
    r = 2.; s = 0.4;lambda = 0.001;
elseif CONFIG==4
    % PERMUTE PHI/PSI
    REVERSE_PSI_PHI = true;
    LOGY_SCALE = false;
    PATH = "./output/config4";
    lambda = 0.0001; r=1.; s= 0.8;
    % try also lambda = 0.001;
elseif CONFIG==5
    PATH = "./output/config5";
    r = 1.; s = 0.1;lambda = 0.001;
elseif CONFIG == 6
    PATH = "./output/config6";
    r = 2.; s = 0.1;lambda = 0.003;
elseif CONFIG == 7
    PATH = "./output/config7";
    r = 1.; s = 0.1;lambda = 0.001;
elseif CONFIG==8
    PATH = "./output/config8";
    r = 2.; s = 0.4;lambda = 0.001;
elseif CONFIG==9
    PATH = "./output/config9";
    r = 2.; s = 0.4;lambda = 0.001;
elseif CONFIG==10
    PATH = "./output/config10";
    r = 2.; s = 0.4;lambda = 0.001;
end



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

mmin = min( min(MSEtrain,[],"all"), min(MSEtest,[],"all"));
mmax = max( max(MSEtrain,[],"all"), max(MSEtest,[],"all"));


%grad = multigradient([8 17 129; 160 56 141; 245 200 79; 243 247 88]./255.);
%grad = multigradient('preset', 'div.cb.BuYlRd.10');
grad = multigradient('preset', 'div.cb.spectral.10');
%grad = multigradient('preset', 'div.cb.GnYlRd.10');
%grad = multigradient('preset', 'div.cb.BuRd.10'); % pas mal
%grad = multigradient('preset', 'div.cb.PuOr.10'); % why not

if REVERSE_PSI_PHI
    label_ = "\phi/\psi";
else
    label_ = "\psi/\phi";
end


%% 3D PLOTS

if PLOT_3D
    tiledlayout(1,2)

    ax1 = nexttile;
    surfc(X,Y,MSEtrain);
    set(gca, 'XScale', 'log')
    if LOGY_SCALE
        set(gca, 'YScale', 'log')
    end
    xlabel('time (t)') 
    ylabel(label_) 
    zlabel('MSE Train')
    caxis([mmin mmax])
    zlim([0. mmax])
    colormap(gca, grad);

    ax2 = nexttile;
    surfc(X,Y,MSEtest);
    set(gca, 'XScale', 'log')
    if LOGY_SCALE
        set(gca, 'YScale', 'log')
    end
    xlabel('time (t)') 
    ylabel(label_) 
    zlabel('MSE Test')
    caxis([mmin mmax])
    zlim([0. mmax])
    colormap(gca, grad);
    
    h = colorbar;
    h.Title.String = "MSE";
    set(gcf,'Position',[0 0 800 0.4*800])
    
    %a =  h.Position;
end

%% 2D PLOTS

if ~PLOT_3D
    tiledlayout(1,2)
    
    ax1 = nexttile;
    pcolor(X,Y,MSEtrain);
    hold on;
    shading interp;
    contour(X,Y,MSEtrain,20,':','LineWidth',0.8,'LineColor','k')
    if LOGY_SCALE
        set(gca, 'YScale', 'log')
    end
    set(gca, 'XScale', 'log')
    xlabel('time (t)') 
    ylabel(label_) 
    title("Train")
    caxis([mmin mmax])
    colormap(gca, grad);

    ax2 = nexttile;
    pcolor(X,Y,MSEtest);
    hold on;
    shading interp;
    [C,h] = contour(X,Y,MSEtest,20,':','LineWidth',0.8,'LineColor','k');
    if LOGY_SCALE
        set(gca, 'YScale', 'log')
    end
    set(gca, 'XScale', 'log')
    xlabel('time (t)') 
    title("Test")
    caxis([mmin mmax])
    colormap(gca, grad);

    h = colorbar;
    h.Title.String = "MSE";
    set(gcf,'Position',[0 0 800 0.3*800])
end