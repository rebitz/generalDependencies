function fH = learningByCondition(rewards,choices,condition)
%% learning during ore/oit

maxLag = 15;
lags = [1:maxLag];    

kernels = NaN(2,maxLag,length(rewards));

for k = 1:length(rewards);

    switches = deal([0, diff(choices{k})~=0]);
    switches(end+1:end+maxLag) = deal(NaN); % pad the end w/ NaNs
    
    chs = choices{k};
    chs(end+1:end+maxLag) = deal(NaN); % pad the end w/ NaNs
    
    rwds = rewards{k};
    states = condition{k};

    for lag = lags
        for i = 1:2
            switch i
                case 1
                    rightState = find(states==1);
                case 2
                    rightState = find(states~=1);
            end
            
            % NEED TO DOUBLE CHECK THIS - specifically which choice is
            % conserved
            sameChoice = chs(rightState) == chs(rightState+lag-1);
            % only look at decisions FOLLOWING same-choice future trials
            %   - this is to make sure the signs line up - rwd now should
            %   cause you to stay in the future
            rightState = rightState(sameChoice);
            a = nanmean(switches(intersect(rightState,find(rwds==1))+lag));
            b = nanmean(switches(intersect(rightState,find(rwds==0))+lag));
            c = nanmean(switches(rightState+lag));
            kernels(i,lag,k) = (b-a)/c; % (rwd-noRwd)/pSwitch
            
            % this should look the same if we reverse it, but it does not
%             sameChoice = choices(rightState) ~= choices(rightState+lag-1);
%             % look for switch backs?
%             rightState = rightState(sameChoice);
%             a = nanmean(switches(intersect(rightState,find(rwds==1))+lag));
%             b = nanmean(switches(intersect(rightState,find(rwds==0))+lag));
%             c = nanmean(switches(rightState+lag));
%             kernels(i,lag,k) = (a-b)/c; % (noRwd-rwd)/pSwitch

        end
    end
end

m = nanmean(kernels,3)';
e = nanstd(kernels,[],3)'./sqrt(size(kernels,3)-1);

fH = figure(); hold on;
set(gca,'FontSize',16);

h = plot(lags*-1,m,'.','MarkerSize',25,'LineStyle','none');
set(h(1),'Color','k'); set(h(2),'Color','b')

legend(flipud(h),'explore','exploit','Location','NorthWest');

h = errorbar([lags'*-1,lags'*-1],m,e);
for i = 1:length(h); removeEBarEnds(h(i)); end
set(h(1),'Color','k'); set(h(2),'Color','b');
set(h,'LineStyle','none');

% add exponential fits
includeOffset = 1;

% now we'll fix an exponential corve to this
xpos = -lags;

if includeOffset
    stopts = [0,2,0];
    fexp = 'a+exp(-x*(1/b))*c'; % exponential function
    % can get rid of the A (offset) & they come together
else
    stopts = [1];
    fexp = 'exp(-x*(1/b))*-1'; % exponential function
    % can get rid of the A (offset) & they come together    
end

opts = fitoptions('method','NonlinearLeastSquares',...
    'StartPoint',stopts);
ftype = fittype(fexp);

clear h;
for i = 1:2
    switch i
        case 1; cStr = 'k'; % oit
        case 2; cStr = 'b'; % ore
    end
    
    m = squeeze(kernels(i,:,:))';
    
    % we have a lot of NaN's because we don't have a lot of obs per subj
    % so, we'll replace them with the mean of the rest of the column
    means = nanmean(m);
    for col = 1:size(m,2)
        m(isnan(m(:,col)),col) = means(col);
    end
    
    tmp = repmat(xpos,size(m,1),1);

    % now fit
    [cf,gof] = fit(-1*tmp(:),m(:),ftype,opts);

    evalpos = min(xpos):range(xpos)/100:max(xpos);
    yhat = cf(-1*evalpos);

    h(i) = plot(evalpos,yhat,'color',cStr);
end

xlim([-maxLag-1 0])

set(h,'LineWidth',2,'LineStyle','--')

ylabel('effect of reward on p(swtich)')
ylabel('previous trial number')
