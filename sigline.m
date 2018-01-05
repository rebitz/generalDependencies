function h = sigline
% SIGLINE Add SIGMOIDAL/LOGISTIC least-squares fit to scatter plot.
%   SIGLINE superimposes the logistic least squares line in the
%   current axes for plots made using PLOT, LINE, SCATTER, or any
%   plot based on these functions.  Any line objects with
%   LineStyles '-', '--', or '.-' are ignored.
% 
%   H = SIGLINE returns the handle to the line object(s) in H.
%   
%   Requires fitLogisticTwoParam.m
%
%   See also LSLINE, POLYFIT, POLYVAL, REFLINE.

%   Adapted from LSLINE 08/12/15, RB Ebitz
%   Copyright 1993-2010 The MathWorks, Inc.


% Find any line objects that are descendents of the axes.
AxCh = get(gca,'Children');
lh = findobj(AxCh,'Type','line');
% Ignore certain continuous lines.
if ~isempty(lh)
    style = get(lh,'LineStyle');
    if ~iscell(style)
        style = cellstr(style);
    end
    ignore = strcmp('-',style) | strcmp('--',style) | strcmp('-.',style);
    lh(ignore) = [];
end

% Find hggroups that are immediate children of the axes, such as plots made
% using SCATTER.
hgh = findobj(AxCh,'flat','Type','hggroup');
% Ignore hggroups that don't expose both XData and YData.
if ~isempty(hgh)
    ignore = ~isprop(hgh,'XData') | ~isprop(hgh,'YData');
    hgh(ignore) = [];
end

hh = [lh;hgh];
numlines = length(hh);
if numlines == 0
    warning(message('stats:lsline:NoLinesFound'));
    hlslines = [];
else
    for k = 1:length(hh)
        if isprop(hh(k),'ZData')
            zdata = get(hh(k),'ZData');
            if ~isempty(zdata) && ~all(zdata(:)==0)
                warning(message('stats:lsline:ZIgnored'));
            end
        end
        % Extract data from the points we want to fit.
        xdat = get(hh(k),'XData'); xdat = xdat(:);
        ydat = get(hh(k),'YData'); ydat = ydat(:);
        ok = ~(isnan(xdat) | isnan(ydat));
        if isprop(hh(k),'Color')
            datacolor = get(hh(k),'Color');
        else
            datacolor = [.75 .75 .75]; % Light Gray
        end
        
        % Fit the points and plot the line.
        % START CHANGE FROM ORIGINAL LSLINE
        
        [params] = fitLogisticTwoParam(xdat,ydat);

        xpos = [min(xdat):range(xdat)/100:max(xdat)];
        logistic = @(params,xdata) (1./(1+exp(-params(1)*(xdata-params(2)))));
        yhat = logistic(params,xpos);
        hlslines(k) = plot(xpos,yhat,'-');
        
        % END CHANGE

        set(hlslines(k),'Color',datacolor);
    end
    set(hlslines,'Tag','lsline');
end

if nargout == 1
    h = hlslines;
end
