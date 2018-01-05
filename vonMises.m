function [fitPut,thisF,thisGof,betterThanMean] = vonMises(thetas,firings,xpos)
% vonMises does a fit of a firing rate vector
%   [fitPut,f,gof] = vonMises(thetas,firings,xToEval)

if mean(thetas) > (4*pi)
    fprintf('\n Input appears to be in degrees, treating as such.\n')
    thetas = thetas * (pi/180);
    xpos = xpos * (pi/180);
    postConvert = 1;
else
    postConvert = 0;
end

fprintf('\n Generating Fit \n')

% gaussian fit is thisun:
fititG = fittype('rMax.*exp((-1/2)*(((x-thMax)/sigma).^2))+c',...
    'coefficients',{'rMax','thMax','sigma','c'});

% von mises fit is thisun:
fititV = fittype('(rMax.*(exp(k*(cos(x-thMax)/sigma))./exp(k)))+c',...
    'coefficients',{'rMax','thMax','sigma','k','c'});

if length(thetas) > size(thetas,1)
    thetas = thetas';
end
if length(firings) > size(firings,1)
    firings = firings';
end

minGof = Inf;

for i = 1:5
    try
    optsG = fitoptions('method','nonlinearleastsquares',...
        'StartPoint',[0,max(firings),15,min(firings)]);
    optsV = fitoptions('method','nonlinearleastsquares',...
        'StartPoint',[5+5*rand,500+100.*randn,2.*randn,3.*randn,min(firings)],...
        'Lower',[0,360,[],[],[]],...
        'MaxFunEvals',300);
    %[0,max(firings),15,1,min(firings)]

    [f,gof] = fit(thetas,firings,fititV,optsV);
    
    if gof.sse < minGof 
        thisF = f;
        thisGof = gof;
        minGof = gof.sse;
    end
    end
end

try,
fitPut = thisF(xpos);

vmPred = thisF(thetas);
vmLogl = -sum(vmPred) + sum(firings.*log(vmPred));

meanVal = nanmean(firings);
mnLogl = -sum(~isnan(firings))*meanVal + sum(firings*log(meanVal));

[aic, bic] = aicbic([vmLogl mnLogl],[5 1], sum(~isnan(firings)));
bicW = exp(-(bic - min(bic))/2)

if bic(1)<bic(2) && bicW(2) < 0.01;
    fprintf('\n SUCESS: \n')
    fprintf('\n VonMises fits better than mean FR! \n')
    fprintf('\n Mean-only model BICw = %d \n',bicW(2))
    betterThanMean = 1;
else
    fprintf('\n CAUTION: \n')
    fprintf('\n Poor fit w/ VonMises \n')
    betterThanMean = 0;
end
catch
    fprintf('\n FATAL ERROR! No function generated! \n')
end