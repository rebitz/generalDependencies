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
        % BELOW IS CHANGE FROM ORIGINAL LSLINE
         % sum of squared perpendicular distances
        R = @(beta) sum((abs(beta(1)*xdat(ok,:)-ydat(ok,:)+beta(2))./sqrt(beta(1).^2+1)).^2);

        % initial guesses - can be anything
        beta0 = -.1; %offset
        beta1 = 0; %slope
        x0 = [beta1 beta0];

        % do minimixation
        options = optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
        [B,resid] = fminunc(R,x0,options);

        % plot it in the current axis
        hlslines(k) = refline(B);

        set(hlslines(k),'Color',datacolor);
    end
    set(hlslines,'Tag','lsline');
end

if nargout == 1
    h = hlslines;
end
