% where X is predictor mx
% y is data to be predicted

link = 'identity';

[b,dev,stats] = glmfit(X,y,...
        'normal','link',link)    
    
yfit = glmval(b,[predict],link);

% get likelihood w/ actual data
mu = glmval(b,X,(link)); st = stats.s
LL = nansum(log(normpdf(y,mu,st)));
n = size(X,2);

% plot model comparisons
figure(); clear h; hold on;
[aic,bic] = aicbic(LL,n,length(y));

k = n;
correction = ((2*k).*(k+1))./(length(y)-k-1);
aicc = aic + correction;

x = 1:length(aic);
h(1) = plot(x,aic,'.b');
h(2) = plot(x,aicc,'.k');
h(3) = plot(x,bic,'.r');
set(h,'LineWidth',2)
xlabel('maximum order of terms included in model','FontSize',16)
ylabel('information criterion','FontSize',16)
set(gca,'FontSize',16)

legend('AIC','AICc','BIC');