function fH = twoMixPlotLogScale_plusNormal(times,xMax)

if nargin < 2; xMax = 40; end

%%
% times = allTimes-1;
xMax = 40;

x = [0:xMax];

% the fitted mixture function:
f2c1 = @(x,theta) (theta(3))*(1./(1+theta(1)))*((theta(1)./(1+theta(1))).^x);
f2c2 = @(x,theta) (theta(4))*(1./(1+theta(2)))*((theta(2)./(1+theta(2))).^x);
f2c3 = @(x,theta) (1-theta(3)-theta(4))*(2*pi)^(-1/2)*(theta(6)^-1).*exp(((-1/2)*(x-theta(5)).^2)./(theta(6).^2));

f2 = @(x,theta) f2c1(x,theta) + ...
    f2c2(x,theta) + ...
    f2c3(x,theta);

fH = figure('Position',[ 440   123   525   575]);
subplot(4,1,1:3); clear h;

% theta = exp2mix(times);

theta = exp2mix_plusNormal(times);

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

% now plot each component
th = plot(x+.5,f2c1(x,theta),'--');
th(2) = plot(x+.5,f2c2(x,theta),'--');
th(3) = plot(x+.5,f2c3(x,theta),'--');
set(th,'Color',[.42 .72 .95])

xlim([min(x) max(x)])
ylim([10^-5 1])

legend([h,th(1)],'mixture','component')
ylabel('probability');

subplot(4,1,4); hold on;
set(gca,'FontSize',16);

% theta = exp2mix_plusNormal(times);
if length(y) < size(y,1);
    y = y';
end
plot(x+1,y.*log(y./f2(x,theta)'));
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