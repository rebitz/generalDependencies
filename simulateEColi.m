% simulate e.coli's path

% we'll save this as a movie:
cd('/Users/becket/Documents/MATLAB/frames')

fieldSize = [500,500]; % w, h

scale = 1;
sigma = 200;
origin = [randi(fieldSize(1)),randi(fieldSize(2))];
startPos = -(origin-mean([[0,0];fieldSize]))+mean([[0,0];fieldSize]);

gauss = @(scale,sigma,origin,evalpos)...
    scale*exp(-((((evalpos(:,1)-origin(1)).^2)./(2*sigma^2)) + ...
    (((evalpos(:,2)-origin(2)).^2)./(2*sigma^2))));

xpos = [1:fieldSize(1)];
ypos = [1:fieldSize(2)];

[xx,yy] = meshgrid(xpos,ypos);
evalpos = [xx(:),yy(:)];

gradient = gauss(scale,sigma,origin,evalpos);
Z = reshape(gradient,fieldSize(1),fieldSize(2));

%% now we'll have the little guy walk

figure(99); axis square; hold on;
set(gca,'FontSize',16);
ylabel('y position');
xlabel('x position');

contour(-Z); contourcmap('gray')
plot(origin(1),origin(2),'ok','MarkerSize',40)

nSteps = 1500;
stepSize = 1;

path = NaN(nSteps,2); % make a place to keep the path
path(1,:) = startPos; % seed the path
pOre = 0.5; pOreStep = 0.01; % how fast does pOre change?
minOre = 0.02; maxOre = 0.95;

for i = 1:nSteps
    if i == 1 || rand < pOre
        % pick a direction at random
        [xStep,yStep] = pol2cart(unifrnd(0,2*pi),stepSize);
        pOre = 0.1; % reset prob of exploring
    else % step forward in the current direction
        xStep = path(i,1)-path(i-1,1);
        yStep = path(i,2)-path(i-1,2);
        
        % and update pOre based on increasing reward gradient
        currentPos = round(path(i,:));
        currentVal = gradient(and(xx == currentPos(1), yy == currentPos(2)));
        
        lastPos = round(path(i-1,:));
        lastVal = gradient(and(xx == lastPos(1), yy == lastPos(2)));
        
        if lastVal > currentVal % if going down the gradient
            pOre = min([(pOre + pOreStep), maxOre]);
        else
            pOre = max([(pOre - pOreStep), minOre]);
        end
    end
    path(i+1,:) = path(i,:) + [xStep,yStep];
    
end

hStep = plot(path(:,1),path(:,2),'k','LineWidth',1.25);

%% now make a movie of this

for i = 1:nSteps
        
    delete(hStep);
    hStep = plot(path([1:i],1),path([1:i],2),'k','LineWidth',1.25);
    hStep(2) = plot(path(i,1),path(i,2),'.k','MarkerSize',40);
%     title(num2str(pOre))
    
    if i < 10
        zPad = '000';
    elseif i < 100
        zPad = '00';
    elseif i < 1000
        zPad = '0';
    else
        zPad = '';
    end
    fNumStr = strcat(zPad,num2str(i));
    saveas(99,strcat('frame',fNumStr,'.jpg'))
end

%% now illustrate the exponential curve being fit?

% we'll save this as a movie:
cd('/Users/becket/Documents/MATLAB/frames')

% maximum likelihood single Exp fit
f1 = @(x,theta) (1./(1+theta))*((theta./(1+theta)).^x);

figure(99); hold on;
set(gca,'FontSize',16);
ylabel('probability');
xlabel('run length');

nReps = 2000;
pStop = 0.2;
xMaxForPlot = 25;
x = [0:1:xMaxForPlot]; % for plotting
multiple = 1;
    
times = [];

for i = 1:nReps
    stopped = 0; j = 0;
    
    while ~stopped
        if rand < pStop
            times = [times; j];
            stopped = 1;
        else
            j = j+1;
        end
    end

    % if we've got 10 new obs, append these to the plot
    if length(times) > 5*multiple
        y = histc(times,x)./length(times);
        % y = y ./ trapz(x,y);
        try; delete(hBar); end
        hBar = bar(x+1,y);
        set(hBar,'LineStyle','none','FaceColor',[.5 .5 .5])
        set(hBar,'BarWidth',1);
        xlim([0 xMaxForPlot]); ylim([0 .35])

        if length(times) > 800
            % plot the single-term exponential
            xeval = [min(x):.1:max(x)];
            theta = nanmean(times); % mle estimate, eq 2.3 mT
            h = plot(xeval+1,f1(xeval,theta),'--k');
            set(h,'LineWidth',3); clear h
        end
%         pause(0.05)

        if multiple < 10
            zPad = '00';
        elseif multiple < 100
            zPad = '0';
        else
            zPad = '';
        end
        
        fNumStr = strcat(zPad,num2str(multiple));
        saveas(99,strcat('frame',fNumStr,'.jpg'))
        
        multiple = multiple+1
    end
end

%% now illustrate different hazards:

% maximum likelihood single Exp fit
f1 = @(x,theta) (1./(1+theta))*((theta./(1+theta)).^x);

figure(99); hold all;
set(gcf,'Position',[0 0 250 200])
set(gca,'FontSize',16);
ylabel('probability');
xlabel('run length');

colors = gray; colors = colors(1:20:end,:);

x = [0:1:xMaxForPlot]; % for plotting

thetaOpts = [3,9]; i = 1;

for theta = thetaOpts

    xeval = [min(x):.1:max(x)];
    h = plot(xeval+1,f1(xeval,theta),'-','Color',colors(i,:));
    set(h,'LineWidth',3); clear h;
    i = i+1;
end

colormap(99,'gray')
lh = legend(cellstr(num2str(1./(thetaOpts+1)')))

% Move the legend box to make room for title
set(lh,'position',get(lh,'position').*[0.92 0.92 1 1])
% Title text
txh = text(0.5,1,'{p(switch):}',...
                'Parent',lh,'VerticalAlign','bottom',...
                'HorizontalAlign','center',...
                'FontSize',16);

ylim([0 .3])
xlim([0 xMaxForPlot])