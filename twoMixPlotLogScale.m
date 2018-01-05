function fH = twoMixPlotLogScale(times,xMax)

if nargin < 2; xMax = 40; end

x = [0:xMax];

% save the fitted mixture function:
f2 = @(x,theta) (theta(3))*(1./(1+theta(1)))*((theta(1)./(1+theta(1))).^x) + ...
    (1-theta(3))*(1./(1+theta(2)))*((theta(2)./(1+theta(2))).^x); % p 25 mT, middle eq

fH = figure('Position',[ 440   123   525   575]);
subplot(4,1,1:3); clear h;

theta = exp2mix(times);

% theta = exp2mix_plusDelta(times);

y = histc(times,x)./length(times);
% y = histc(times,x);
% y = y./trapz(x,y);

semilogy(x+.5,y,'.','Color',[.5 .5 .5],...
    'MarkerSize',30,'LineStyle','none');

hold on;
set(gca,'FontSize',16);

% fZ2, shift the plotting over, but eval starting at 0
h = semilogy(x+.5,f2(x,theta),'-','LineWidth',3); hold on; 
set(h,'Color',[.42 .72 .95])

th2 = theta;

% now the single exp, maximum likelihood version of the same
theta = nanmean(times); % mle estimate, eq 2.3 mT
% marginal distribution of times: (eq 2.1 mT)
f1 = @(x,theta) (1./(1+theta))*((theta./(1+theta)).^x);
h(2) = semilogy(x+.5,f1(x,theta),'--k','LineWidth',3);

th = plot(x+.5,th2(3)*f1(x,th2(1)),'--');
th(2) = plot(x+.5,(1-th2(3))*f1(x,th2(2)),'--');
set(th,'Color',[.42 .72 .95])

xlim([min(x) max(x)])
ylim([10^-5 1])

legend(h,'2mix','1exp')
ylabel('probability');


subplot(4,1,4); hold on;
set(gca,'FontSize',16);

theta = exp2mix(times);
if length(y) < size(y,1);
    y = y';
end
plot(x+1,y.*log(y./f2(x,theta)));
% how much does observed P deviate from predicted P?
% extra information needed to code each message if our model is wrong

h = line([min(x) max(x)],[0 0]);
set(h,'Color','k');
ylabel('K-L divergence');
xlabel('run length');

xlim([min(x) max(x)])
% for i = x
%     plot(i,sum(times==i)*log(1/f2(i,theta))); % self information of observation
% end