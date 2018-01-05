% TO DO:
% 3. adapt this same thing to work w/ the MGS data
% 4. re-process all existing MIB and MGS files

%%

% parse Saccades version 3:
%   - improved alignment of offset timestamps relative to v2
%
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

% set some parameters
normTo = 1000;
plotIt = 0;

smoothby = 3;
preWindow = 50; % time b4/after points to avg across for saccade ID
peakthresh = 100; % must cross this to be a peak
saccthresh = 20; % soft implementation - range of saccade starts
brokeFixThresh = 80; % velo to count as a broken fixation
errorWindow = 10; %area around the trial-identified fixation to count
targVectError = 50; % acceptable range of initial saccades, deg

errScale = 1.25; % multiply spatial errors by this value
fixLengthCheat = 60; % if nbins in fix is off by < this, don't worry

% now pull some stuff from our task_data
origin = px2deg(task_data.fix_loc,env);
fixErr = px2deg(task_data.fix_err,env);
tErr = px2deg(task_data.targ_err,env);

% pre/postpadding for each tr:
if ~exist('MGS') || ~MGS
    pres = [trials.fixAcq]-[trials.trialstart];
    posts = [trials.ITIend]-([trials.targAcq]+[task_data.targHoldTime]); % PAY ATTNS! ITI may come at start rather than end of tr in some sessions
    superPres = ([trials.targAcq]+[task_data.targHoldTime])-[trials.trialstart];
else
    pres = [trials.fixAcq]-[trials.trialStart];
    posts = [trials.ITIend]-([trials.targAcq]+[task_data.targHoldTime]); % PAY ATTNS! ITI may come at start rather than end of tr in some sessions
    superPres = ([trials.targAcq]+[task_data.targHoldTime])-[trials.trialStart];    
end
    
% initialize output of n bins we had to shift by
binsOff = NaN(1,length(trials));

%%

for i = 1:length(trials)
try
plotIt = 1;
%%
    if trials(i).correct == 1 && sum(~isnan(trials(i).eye(:,1)))>10

        i % echo where we're at
    % fix the offset in eye position time stamps:
    
        % first, find the t elapsed btw fixAcq and targ Acq
        fixLen = [trials(i).fixOff]-[trials(i).fixAcq];
        fixToTarg = [trials(i).targAcq]-[trials(i).fixAcq];
        
        if ~exist('MGS') || ~MGS
            % then chosen target locations
            switch trials(i).choice
                case 1; tLoc = px2deg(task_data.t1origin,env);
                case 2; tLoc = px2deg(task_data.t2origin,env);
                case 3; tLoc = px2deg(task_data.t3origin,env);
            end
        else % for MGS data
            % RectLeft=1, RectTop=2, RectRight=3, RectBottom=4.
            tLoc = [mean(trials(1).targRect([1,3])),...
                mean(trials(1).targRect([2,4]))];
            tLoc = px2deg(tLoc,env);
        end
        
        % pull the eye data into vectors
        tstamps = trials(i).eye(:,1);
        xdata = trials(i).eye(:,2);
        ydata = trials(i).eye(:,3);
        % deal with unreasonable values
        % xdata(ydata>3000) = NaN; ydata(ydata>3000) = NaN;

        % get first derivative of the data, also divide by timing
        xvelo = diff(xdata)./diff(tstamps);
        xvelo = xvelo(:); % arrange into column
        yvelo = diff(ydata)./diff(tstamps);
        yvelo = yvelo(:);

        % smooth the velocity information for easier saccade finding
        velo = sqrt(xvelo.^2 + yvelo.^2);
        smoothvelo = gsmooth(velo,smoothby);
        
        % now some logical-ing
        tmpVelo = [0, smoothvelo];
        veloIdx = tmpVelo < brokeFixThresh;
        
        % all the points that could be fixation fixations
        xFixIdx = and(xdata <= origin(1)+(fixErr)*errScale, xdata >= origin(1)-(fixErr)*errScale);
        yFixIdx = and(ydata <= origin(2)+(fixErr)*errScale, ydata >= origin(2)-(fixErr)*errScale);
        fixIdx = and(and(xFixIdx,yFixIdx),veloIdx');
        
        % skip the pre-period:
        fixIdx(1:ceil(pres(i)/nanmean(diff(tstamps)))) = false;
        
        % now cycle through all the fixations
        tmp = diff(fixIdx); % if we find more than one epoch
        if sum(tmp==1) > 1;
            fixNum = zeros(1,length(fixIdx)); % preallocate
            running = 1;
            for q = 2:length(fixIdx)
                if fixIdx(q) && fixIdx(q-1);
                    fixNum(q) = running; % if no change
                elseif fixIdx(q)
                    running = running+1;
                    fixNum(q) = running; % if a change
                end
            end
            % only keep the fixation that was long enough to be our guy
            tmp = unique(fixNum); tmp = tmp(2:end); % skip the 0's
            goodFix = NaN;
            for q = tmp;
                if sum(fixNum==q) > (fixLen/nanmean(diff(tstamps))-fixLengthCheat)
                    goodFix = q;
                    break % stop with the FIRST good fixation
                end
            end
            if ~isnan(goodFix); fixIdx = [fixNum==goodFix]; end
        else
            fixNum = fixIdx;
        end
        
        % get rid of any fixation points that are ridiculously late
        fixIdx(floor(superPres/nanmean(diff(tstamps))):end) = false;
        
        % all the points that could be target fixations
        xTargIdx = and(xdata <= tLoc(1)+(tErr)*errScale, xdata >= tLoc(1)-(tErr)*errScale);
        yTargIdx = and(ydata <= tLoc(2)+(tErr)*errScale, ydata >= tLoc(2)-(tErr)*errScale);
        targIdx = and(and(xTargIdx,yTargIdx),veloIdx');
        
        % require that these occur after the ID'ed fixation
        targIdx(1:find(fixIdx,1,'last')) = false;
        
        % but before the end points of the trial
        targIdx(end-ceil(posts(i)/nanmean(diff(tstamps))):end) = false;
        
        % now take a look at the damage
        if plotIt
            figure(98);
            subplot(2,1,1); hold on; plot(xdata)
            tmp = xlim; title('xpos')
            h = line([tmp(1),tmp(2)],[origin(1)+(fixErr)*errScale origin(1)+(fixErr)*errScale])
            h(2) = line([tmp(1),tmp(2)],[origin(1)-(fixErr)*errScale origin(1)-(fixErr)*errScale])
            set(h,'Color','g','LineStyle','--');
%             h = line([tmp(1),tmp(2)],[tLoc(1)+(tErr)*errScale tLoc(1)+(tErr)*errScale])
%             h(2) = line([tmp(1),tmp(2)],[tLoc(1)-(tErr)*errScale tLoc(1)-(tErr)*errScale])
%             set(h,'Color','m');

            subplot(2,1,2); hold on; plot(ydata)
            tmp = xlim; title('ypos')
            line([tmp(1),tmp(2)],[origin(2)+(fixErr)*errScale origin(2)+(fixErr)*errScale])
            line([tmp(1),tmp(2)],[origin(2)-(fixErr)*errScale origin(2)-(fixErr)*errScale])
        end
        
        % check to see if this is WAY off from the targ Acq point, if so,
        % then shift the timepoints again to make this congruent:
        evtBin = find(trials(i).eye(:,1) > [trials(i).targAcq],1,'first');
        binOffset = find(targIdx,1,'first')-evtBin;
        timeOffset = binOffset*nanmean(diff(tstamps));
        
        % adjust the eye field time stamps accordingly
        if abs(binOffset)>5; % if we seem to be > 1 bin off
            lateShift = 1;
            %fprintf('\n targAcq is off; adjusting timestamps by %d',timeOffset)
            trials(i).eye(:,1) = tstamps-timeOffset;
%             startIdx = startIdx+binOffset;
%             endIdx = endIdx+binOffset;
        end
        tstamps = trials(i).eye(:,1);
        
        % save this info for later so we can check our shifts:
        try, binsOff(i) = binOffset; end
        
        % quick plot to check this out        
        if plotIt;
            figure(99); clf;
            subplot(3,1,1:3); hold on;

            plot(smoothvelo); title(i);
            temp = xlim;
            
            % thresholds
            line([temp(1),temp(2)],[peakthresh peakthresh]);
            line([temp(1),temp(2)],[saccthresh saccthresh]);
            
            % targ on marker
            tmp = find(trials(i).eye(:,1) > [trials(i).targOn],1,'first');
            a = line([tmp(1),tmp(1)],[0 100]);
            set(a,'Color','m','LineStyle','--')
            % targ acq marker
            tmp = find(trials(i).eye(:,1) > [trials(i).targAcq],1,'first');
            q = line([tmp(1),tmp(1)],[0 100]);
            set(q,'Color','m','LineWidth',3)
            % fix on marker
            tmp = find(trials(i).eye(:,1) > [trials(i).fixOn],1,'first');
            b = line([tmp(1),tmp(1)],[0 100]);
            set(b,'Color','g','LineStyle','--')
            % fix acq marker
            tmp = find(trials(i).eye(:,1) > [trials(i).fixAcq],1,'first');
            e = line([tmp(1),tmp(1)],[0 100]);
            set(e,'Color','g','LineWidth',3)
            
            j = plot(fixIdx*150);
            set(j,'Color','g')
            
            k = plot(targIdx*200);
            set(k,'Color','m')
            
            % identified fix w/ length
            temp = ylim;
            fu = find(fixIdx,1,'first');
            r = line([fu fu+fixLen/nanmean(diff(tstamps))],[temp(1)+range(temp)/2 temp(1)+range(temp)/2]);
            set(r,'Color','k');
            
            legend([a,q,b,e,j,k],'targOn','targAcq','fixOn','fixAcq','fixFix','targFix')
            
        end
        
        %% now for the saccade
        %  we'll find the start and end and traj. of the saccade
                
        % find start of the saccade based on first point after fix to cross
        % the peak threshold
        temp = find(smoothvelo(find(fixIdx,1,'last'):end) > peakthresh,1,'first');
        tempIdx = temp+find(fixIdx,1,'last'); % get index of that point
        % then walk back to the start of the high velo period to find the
        % actual starting point
        % startIdx = find(round(diff(smoothvelo(1:tempIdx))*normTo)/normTo <= 0,1,'last');
        % instead, we're going to use low-threshold crossing
        startIdx = find(smoothvelo(1:tempIdx) <= saccthresh,1,'last');
        
        % saccade end point,
        % find the point after the saccade with non-decreasing velocity
        temp = find(smoothvelo(tempIdx:end) < peakthresh,1,'first');
        tempIdx = temp+tempIdx; % get index of that point

        endIdx = find(round(diff(smoothvelo(tempIdx:end))*normTo)/normTo >= 0,1,'first');
        endIdx = endIdx+tempIdx;

        endIdx = find(smoothvelo(tempIdx:endIdx) >= saccthresh,1,'last');
        endIdx = endIdx+tempIdx;

        if plotIt
            for fig = 1:3
                switch fig
                    case 1; figure(99)
                    case 2; figure(98); subplot(2,1,1);
                    case 3; figure(98); subplot(2,1,2);
                end
                
                temp = ylim;
                c = line([endIdx endIdx],[temp(1) temp(2)]);
                set(c,'Color','k');
                d = line([startIdx startIdx],[temp(1) temp(2)]);
                set(d,'Color','r');
                drawnow;
            end
                pause(0.1);
        end
        
        % now for the trajectory of the saccade:
        % we'll need windows instead of points
        startPts = [startIdx-preWindow:startIdx];
        endPts = [endIdx:endIdx+preWindow];
        
        % now get the angle of the saccade
        x1 = nanmedian(trials(i).eye(startPts,2));
        y1 = nanmedian(trials(i).eye(startPts,3));
        
        xEnd = [nanmedian(trials(i).eye(endPts,2))];
        yEnd = [nanmedian(trials(i).eye(endPts,3))];
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
        
    else
            trials(i).saccadeStart = NaN;
            trials(i).saccadeEnd = NaN;
            trials(i).saccadeTheta = NaN;
            trials(i).errantSaccade = NaN;
            trials(i).peakVelocity = NaN;
    end
catch
    trials(i).saccadeStart = NaN;
    trials(i).saccadeEnd = NaN;
    trials(i).saccadeTheta = NaN;
    trials(i).errantSaccade = NaN;
    trials(i).peakVelocity = NaN;
    fprintf('UNKNOWN ERROR, trial number %i',i)
end
end

fprintf('\n saccades identified. \n')