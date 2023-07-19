function yfit = regf(Xtrain,ytrain,Xtest)
b = regress(ytrain,Xtrain);
yfit = Xtest*b;
end