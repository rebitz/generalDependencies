

% parse Saccades
% take data struct "trials", which contains a field "eye"
% eye must be a m x 4 matrix, where columns = time, x, y, pupil
%   output from "edfProcess" and from "concatMIBsubj"
% IMPORTANT - this only works on MATLAB generated data - not CTX!!!

if exist('trials') == 1
    if isfield(trials,'eye')
        fprintf('\n passed initial checks \n')
        fprintf('\n beginning saccade identification. \n')
    end
end

errants = [];

normTo = 1000;

plotit = 0;

fixPad = 0; % pad the fixation so we skip the saccade, better est of fixation origin
smoothby = 10;
preWindow = 300; % NOT IMPLEMENTED
peakthresh = 100; % must cross this to be a peak
saccthresh = 20; % soft implementation - range of saccade starts
errorWindow = 10; %area around the trial-identified fixation to count
targVectError = 50; % acceptable range of initial saccades, deg

%%
plotit = 1;
for i = 1:length(trials)
% try,        
    % fix the offset in eye position time stamps:

    if trials(i).correct == 1 && sum(~isnan(trials(i).eye(:,1)))>10
        
        % first, operate on timing - find the block between fix acquire and
        % target onset in the eye data
        fixIdx = and(trials(i).eye(:,1) > [trials(i).fixAcq]+fixPad,trials(i).eye(:,1) < [trials(i).targOn]);
        origin = [nanmean(trials(i).eye(fixIdx,2)),nanmean(trials(i).eye(fixIdx,3))];

        % pull the eye data into vectors
        tstamps = trials(i).eye(:,1);
        xdata = trials(i).eye(:,2);
        ydata = trials(i).eye(:,3);
        % deal with unreasonable values
%         xdata(ydata>3000) = NaN; ydata(ydata>3000) = NaN;
        
        % first derivative of the data, also divide by timing
        xvelo = diff(xdata)./diff(tstamps);
        xvelo = xvelo(:); % arrange into column
        yvelo = diff(ydata)./diff(tstamps);
        yvelo = yvelo(:);

        % smooth the velocity information for easier saccade finding
        velo = sqrt(xvelo.^2 + yvelo.^2);
        smoothvelo = gsmooth(velo,smoothby);
        
        % find start of the saccade based on first point after fix to cross
        % the peak threshold
        temp = find(smoothvelo(find(fixIdx,1,'last'):end) > peakthresh,1,'first');
        tempIdx = temp+find(fixIdx,1,'last'); % get index of that point
        % then walk back to the start of the high velo period to find the
        % actual starting point
        startIdx = find(round(diff(smoothvelo(1:tempIdx))*normTo)/normTo <= 0,1,'last');
        
        % saccade end point,
        % find the point after the saccade with non-decreasing velocity
        temp = find(smoothvelo(tempIdx:end) < peakthresh,1,'first');
        tempIdx = temp+tempIdx; % get index of that point

        endIdx = find(round(diff(smoothvelo(tempIdx:end))*normTo)/normTo >= 0,1,'first');
        endIdx = endIdx+tempIdx;
        
        % now get the angle of the saccade
        x1 = origin(1); xEnd = trials(i).eye(endIdx,2);
        y1 = origin(2); yEnd = trials(i).eye(endIdx,3);
        x = [xEnd-x1]; y = [yEnd-y1]; % orient relative to the origin
        
        saccadeTheta = cart2pol(x,y); % convert to polar
        saccadeTheta = saccadeTheta*(180/pi); % convert to degrees
        if saccadeTheta<0; saccadeTheta = saccadeTheta+360; end % correct neg. output
        %saccadeTheta = 360-saccadeTheta;
        trials(i).saccadeTheta = saccadeTheta;
        
        targTheta = NaN;
        if ~exist('thetaCheck','var')
            if isfield(trials,'theta')
                % eccentric target position
                targTheta = (trials(i).theta);
                thetaCheck = 1;
            else
                thetaCheck = 0;
            end
        end
        
        if isempty(endIdx) || isempty(startIdx)
%             errants = [errants; i];
            trials(i).saccadeStart = NaN;
            trials(i).saccadeEnd = NaN;
            trials(i).errantSaccade = NaN;
            trials(i).peakVelocity = NaN;
        %then, check to make sure that the angle of the saccade matches the angle of the target    
        elseif thetaCheck && (saccadeTheta > targTheta+targVectError && ~(saccadeTheta-360 > targTheta-targVectError)) ...
                || saccadeTheta < targTheta-targVectError
            fu = [i targTheta saccadeTheta]
            errants = [errants; fu];
            trials(i).errantSaccade = 1;
            trials(i).saccadeStart = NaN;
            trials(i).saccadeEnd = NaN;
            trials(i).peakVelocity = NaN;
        elseif ~isempty(startIdx) && abs(trials(i).eye(startIdx,2)-origin(1)) < errorWindow && abs(trials(i).eye(startIdx,3)-origin(2)) < errorWindow 
            trials(i).saccadeStart = trials(i).eye(startIdx,1);
            trials(i).saccadeEnd = trials(i).eye(endIdx,1);
            trials(i).errantSaccade = 0;
            trials(i).peakVelocity = max(smoothvelo(startIdx:endIdx));
        elseif trials(i).saccadeEnd < trials(i).saccadeStart
            fprintf('\n end listed as before start, trial number %d \n',i)
            trials(i).saccadeStart = NaN;
            trials(i).saccadeEnd = NaN;
            trials(i).errantSaccade = NaN;
            trials(i).peakVelocity = NaN;
        else
            trials(i).saccadeStart = NaN;
            trials(i).saccadeEnd = NaN;
            trials(i).errantSaccade = NaN;
            trials(i).peakVelocity = NaN;
        end
        
        if find(fixIdx,1,'last') > startIdx+10;
            plotit = 1;
            trials(i).saccadeStart = NaN;
            trials(i).saccadeEnd = NaN;
            trials(i).errantSaccade = NaN;
            trials(i).peakVelocity = NaN;
        end
        
        if plotit
            figure(99);  clf; hold on; %subplot(2,1,1); 
%             set(gcf,'Position',[2400 400 350 300]);
            plot(smoothvelo);
            title(i)
            temp = xlim;
            line([temp(1),temp(2)],[peakthresh peakthresh]);
            line([temp(1),temp(2)],[saccthresh saccthresh]);
            
            % targ on marker
            tmp = find(trials(i).eye(:,1) > [trials(i).targOn],1,'first');
            h = line([tmp(1),tmp(1)],[0 100]);
            set(h,'Color','g')
            % fix on marker
            tmp = find(trials(i).eye(:,1) > [trials(i).fixOn],1,'first');
            h = line([tmp(1),tmp(1)],[0 100]);
            set(h,'Color','g')
            
            temp = ylim;
            h = line([endIdx endIdx],[temp(1) temp(2)]);
            set(h,'Color','k');
            h = line([startIdx startIdx],[temp(1) temp(2)]);
            set(h,'Color','r');
            
            h = line([find(fixIdx,1,'first') find(fixIdx,1,'last')],[temp(1)+range(temp)/2 temp(1)+range(temp)/2]);
            pause(.1)
%             subplot(2,1,2); hold on;
%             plot(diff(smoothvelo))
%             plot(diff(diff(smoothvelo)))
plotit = 0;
        end
    else
            trials(i).saccadeStart = NaN;
            trials(i).saccadeEnd = NaN;
            trials(i).saccadeTheta = NaN;
            trials(i).errantSaccade = NaN;
            trials(i).peakVelocity = NaN;
    end
%     catch
%     fprintf('\n skipping trial number %d, UNKNOWN ERROR \n',i)
%     
%             trials(i).saccadeStart = NaN;
%             trials(i).saccadeEnd = NaN;
%             trials(i).saccadeTheta = NaN;
%             trials(i).errantSaccade = NaN;
%             trials(i).peakVelocity = NaN;
% end

end

fprintf('\n saccades identified. \n')

%% additional visualization of the trial eye data
%for i = 10:20
plotIt = 1;
if ~exist('trialnum','var')
    trialnum = 1
end

if plotIt

%     trialnum = trialnum+1
% trialnum = 10
i = 20
%     list = find([trials.correct]==1);
%     i = list(trialnum);
    
    startIdx = find(trials(i).eye(:,1) > trials(i).saccadeStart,1,'first');
    endIdx = find(trials(i).eye(:,1) > trials(i).saccadeEnd,1,'first');
    
    xdata = trials(i).eye(:,2);
    ydata = trials(i).eye(:,3);
    
    xdata(ydata>3000) = NaN; ydata(ydata>3000) = NaN;
    
    figure(i); clf;
    % set(gcf,'Position',[1900 400 350 600]);
    subplot(2,1,1); hold on;
    plot(xdata,ydata)
    plot(xdata(1),ydata(1),'.b','MarkerSize',20)
    plot(xdata(startIdx),ydata(startIdx),'.r')
    plot(xdata(endIdx),ydata(endIdx),'.k')
    title(i)
%     ylim([-200 200]); xlim([-200 200])
    subplot(2,1,2); hold on;
    plot(trials(i).eye(:,1),xdata,'-b')
    plot(trials(i).eye(:,1),ydata,'-g')
    try
        x = trials(i).eye(startIdx,1);
        y = ylim;
        h = line([x x],[y(1) y(2)]);
        set(h,'Color','r')
    end
    try
        x = trials(i).eye(endIdx,1);
        y = ylim;
        h = line([x x],[y(1) y(2)]);
        set(h,'Color','k')
    end
    legend({'xpos','ypos','start','end'})
    title('over T')
end
%end