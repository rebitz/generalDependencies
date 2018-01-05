function [params] = fitLogisticThreeParam(x,y,plotIt)

    if nargin<3
        plotIt = 0;
    end
    
    % stretch short or long
    tstretch = [];
    if max(x)-min(x) < 1.e-4 || max(x)-min(x) > 1e5;
        tstretch = 1./(max(x) - min(x));
        x = x*tstretch;
    end

    dataN = length(x);

    %model's parameters to be fitted
    mlpar = NaN(1,7);

    options = optimset('MaxFunEvals',2000,'MaxIter',700,...
        'TolX',0.000001,'TolFun',0.000001,'Display','off');

    fvalmin = dataN; % length of data to be fitted??
% keyboard();
    for i=1:20
                % slope  intercept  scale
        initP = [rand(1) rand(1) 1];

        [guess fval exitflag output] = fminsearch(@(params) logisticCost(params,x,y),...
            initP,options);
%         [guess fval res exitflag output] = lsqcurvefit(@logistic,initP,x,y,[],[],options);
% params = lsqcurvefit(@generalLogistic,[.5,0,A,K],selection,M')

        if fval <= fvalmin %&& exitflag==1 % save output if we get improvement
            fvalmin = fval;
            mlpar(1:length(initP)) = guess(1:length(initP));
            mlpar(5) = fval;
            mlpar(6) = fval./dataN;
            mlpar(7) = output.iterations;
            mlpar(8) = exitflag;
        end
    end

    params = mlpar(1:length(initP));
    
    if ~isempty(tstretch)
        params(1) = params(1)/tstretch;
        params(2) = params(2)*tstretch;
    end  

    if plotIt == 1
        figure(); hold on;
        plot(x,y,'.k')
        xpos = [min(x):.01:max(x)];
        
        [~,id] = histc(x,xpos);
        for i = 1:length(unique(id))-1;
            h = plot(xpos(i),nanmean(y(id == i)),'x');
            set(h,'Color','r')
        end
        
        h = plot(xpos,logistic(params,xpos));
        set(h,'Color','r');
    end

    function z = logisticCost(params,x,y)
        r = params(1);% slope
        x0 = params(2);% intercept
        s = params(3);% scale

        h0 = logistic(params,x);
        
%         keyboard();
%         h = 1./(1+exp(-r*(x-x0)));
%         numerator = -dot(h,y)^2;
%         denominator = dot(h,h);
%         z = numerator/denominator;

        z = nansum((h0-y).^2);
    end

    function y = logistic(params,xdata)
        r = params(1); % slope
        x0 = params(2); % intercept
        s = params(3); % intercept
%         A = 0; % minimum
%         K = 1; % maximum

        y = ((1-s)./2) + (s./(1+exp(-r*(xdata-x0))));
       
    end 
end