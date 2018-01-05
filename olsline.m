function h = olsline
% OLSLINE Add TOTAL least-squares fit line to scatter plot.
%   OLSLINE superimposes the orthoganal least squares line in the
%   current axes for plots made using PLOT, LINE, SCATTER, or any
%   plot based on these functions.  Any line objects with
%   LineStyles '-', '--', or '.-' are ignored.
% 
%   H = OLSLINE returns the handle to the line object(s) in H.
%   
%   See also LSLINE, POLYFIT, POLYVAL, REFLINE.

%   Adapted from LSLINE 10/04/13, RB Ebitz
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
        
        % Demming regression w/ delta set to 1
        xbar = mean(xdat);
        ybar = mean(ydat);
        sXX = sum((xdat-xbar).^2)./(length(xdat)-1);
        sXY = sum((xdat-xbar).*(ydat-ybar))./(length(xdat)-1);
        sYY = sum((ydat-ybar).^2)./(length(ydat)-1);

        B1 = sYY-sXX+sqrt((sYY-sXX).^2+4*sXY.^2);
        B1 = B1./(2*sXY);
        B0 = ybar-B1*xbar;
        B = [B1,B0];
        
        % END CHANGE
        
        % plot it in the current axis
        hlslines(k) = refline(B);

        set(hlslines(k),'Color',datacolor);
    end
    set(hlslines,'Tag','lsline');
end

if nargout == 1
    h = hlslines;
end
