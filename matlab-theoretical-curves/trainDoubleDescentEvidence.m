%  training regime
PATH = "./output/config2";

res = importResult(PATH+"/MSE_"+31+"_log.bin");

mse = timeEvolution(400, 2, 0.4, 0.1, res);
semilogx(mse.Tspace, mse.MSEtrain); 
writetable(mse, "./output/doubleTrainingDescent.csv")



%% test regime

res = importResult(PATH+"/MSE_"+14+"_log.bin");
mse = timeEvolution(400, 2, 0.4, 0.0001, res);

semilogx(mse.Tspace, mse.MSEtest); 
writetable(mse, "./output/doubleTestDescent.csv")
