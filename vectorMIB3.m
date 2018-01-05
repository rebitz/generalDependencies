%%
% this really only works if the direction of motion is
% perpendicular to the saccade - which it's not if the 
% same stimulus is used for all thetas

MATLAB = 1; % set to 1 if data comes from MATLAB
meanRegister = 1; % align to some within condition factor
    % (else, default to register to veridical positions)
    toSaccades = 0; % to center of saccades for that condition
        runningSaccadeMean = 0; filterLength = 150; % DO NOT USE!!
    toCenterLocation = 1; % to center targ loc. for that condition
rotate = 0; % don't use this! deals with actual saccade position now

plotit = 1;

[trials.xChangeRaw] = deal(NaN);
[trials.yChangeRaw] = deal(NaN);

% formula for vector magnitude.
getMag = @(v) (sqrt(v(1,:).^2+v(2,:).^2));

for i = 1:length(trials);
    
    i
    
    % first grab timestamps
    startIdx = trials(i).saccadeStart;
    endIdx = trials(i).saccadeEnd;
    
    if trials(i).correct == 1 && ~isnan(startIdx) && ~isnan(endIdx)

        % EYE POSITION IS FLIPPED UP/DOWN!!!
        eye = trials(i).eye;

        % convert to indices
        startIdx = find(eye(:,1) >= startIdx,1,'first');
        endIdx = find(eye(:,1) <= endIdx,1,'last');

        % starts and ends
        yPosStart = eye(startIdx,3);
        yPosEnd = eye(endIdx,3);
        xPosStart = eye(startIdx,2); % deal with the L/R flipping
        xPosEnd = eye(endIdx,2);

        % get saccade vector
        xPos = [xPosEnd - xPosStart];% xPos = xPosEnd;
        yPos = [yPosEnd - yPosStart];% yPos = yPosEnd;

        if rotate
            % extract target position and align to matlab CCW space, rather than
            % ctx's CW convention
            realTh = trials(i).theta+180;

            % convert to radians
            realTh = realTh/(180/pi);

            % rotate the coordinate system so that all saccades go straight R
            % where that is the actual position of the target in space
            rotPt = @(x,y,th) [x*cos(th) - y*sin(th), x*sin(th) + y*cos(th)];
            temp = rotPt(xPos,yPos,realTh);
            xPos = temp(1); yPos = temp(2);
        end
        
        trials(i).yChangeRaw = yPos;
        trials(i).xChangeRaw = xPos;
    end
     
end

v = [[trials.xChangeRaw];[trials.yChangeRaw]];
stepsPerSaccade = nanmedian(getMag(v));

if isfield(trials,'targEcc')
    degreesPerSaccade = nanmean([trials.targEcc]);
else
    degreesPerSaccade = 8; % SET ME CORRECTLY
end
degPerStep = degreesPerSaccade / stepsPerSaccade;

% xChange and yChange are in degrees(ish)
temp = num2cell([trials.yChangeRaw] .* degPerStep);
[trials.yChangeRaw] = deal(temp{:});

temp = ([trials.xChangeRaw] .* degPerStep);
temp = num2cell(temp);
[trials.xChangeRaw] = deal(temp{:});

% get rid of truly ridiculuous values
idx = or(abs([trials.xChangeRaw])>3000,abs([trials.yChangeRaw])>3000);
[trials(idx).xChangeRaw] = deal(NaN);
[trials(idx).yChangeRaw] = deal(NaN);

% first get ready to run the analysis
[trials.projection] = deal(NaN);
holdUnits = NaN(2,length(trials));

if ~MATLAB
    thetas = 180+[trials.theta];
else
    thetas = 180+[trials.theta];
end
thetas = ((thetas)/(180/pi));

% get the target locations
[x,y] = pol2cart(thetas,degreesPerSaccade);

if meanRegister
    % ALTERNATIVELY: get the middle part of saccades to that actual position
    tmp = unique([trials.condition]);
    for i = 1:length(tmp)
        
        if toSaccades
            idx = [trials.condition] == tmp(i);
            
            if ~runningSaccadeMean
                x1 = nanmean([trials(idx).xChangeRaw]);
                y1 = nanmean([trials(idx).yChangeRaw]);
            else
                
                disp('DOES NOT PLAY WELL WITH OUTLIERS!!!!');
                x1 = gsmooth([trials(idx).xChangeRaw],filterLength);
                y1 = gsmooth([trials(idx).yChangeRaw],filterLength);
            end
            
            x(idx) = x1;
            y(idx) = y1;
        elseif toCenterLocation
            idx = [trials.condition] == tmp(i);
            th = round(nanmean([trials(idx).theta]));
            if ~MATLAB; th = 180+th; end
            th = ((th)/(180/pi));
            
            thetas(idx) = deal(th);
            [x,y] = pol2cart(thetas,degreesPerSaccade);
        end
    end
end

for i = 1:2 % for the two directions of motion

    switch i
        case 1
            % select only upward going trials
            idx = [trials.up] == 1;
            locations = [y;-x]; % perpendicular vectors
        case 2
            % now for downward deflections
            idx = [trials.up] == 0;
            locations = [-y;x]; % perpendicular vectors
    end
    
    % unit vectors are the mag == 1 vectors perpendicular to 
    %   a straight saccade to the target and congruent with motion
    mag = getMag(locations); mag = [mag;mag];
    unitVectors = locations./mag;

    % sanity check for perpendicularity
    nanmean(dot([x;y],unitVectors)) % should == 0

    % pull out the locations of the actual saccades
    x1 = [trials(idx).xChangeRaw];
    y1 = [trials(idx).yChangeRaw];

    % register to the ideal (center of mass) saccade
    xOff = x1-x(idx); yOff = y1-y(idx);

    % project relevant saccade offsets to motion direction unit vectors
    devs = dot(unitVectors(:,idx),[xOff;yOff]);

    % save these to the trial struct
    devs = num2cell(devs);
    [trials(idx).projection] = deal(devs{:});

    % and hang onto the unit vectors we generated for plotting
    holdUnits(:,idx) = unitVectors(:,idx);
    
end

%% deal with undershot saccades, defined as < mu+3std of the rest

x = [trials.xChangeRaw];
y = [trials.yChangeRaw];

[t,r] = cart2pol(x,y);
badIdx = or(r>nanmean(r)+(3*nanstd(r)),r<nanmean(r)-(3*nanstd(r)));

[trials(badIdx).xChangeRaw] = deal(NaN);
[trials(badIdx).yChangeRaw] = deal(NaN);
[trials(badIdx).projection] = deal(NaN);


% if removeOutliers
%     tmp = unique([trials.condition]);
%     tmp = tmp(~isnan(tmp));
%     
%     for i = 1:length(tmp)
%         idx = [trials.condition] == tmp(i);
%         xM = nanmean([trials(idx).xChangeRaw]);
%         yM = nanmean([trials(idx).yChangeRaw]);
%         
%     end
% end
%% plot an example trial
plotit = 0
if plotit
which = 78;
figure(); hold on; axis square

theta = [trials(which).theta]+180; theta = theta/(180/pi);
[x,y] = pol2cart(theta,degreesPerSaccade);
perpVec = [y;-x]; % perpendicular vectors

% mag = getMag(perpVec); mag = [mag;mag];
% unitVectors = perpVec./mag;

h = line([x x+holdUnits(1,which)],[y,y+holdUnits(2,which)]);
set(h,'Color','k','LineWidth',2)
    
% if [trials(which).up] == 1
%     h = line([x x+unitVectors(1)],[y,y+unitVectors(2)]);
%     set(h,'Color','k','LineWidth',5)
% else
%     h = line([x x-unitVectors(1)],[y,y-unitVectors(2)]);
%     set(h,'Color','r','LineWidth',2)
% end

line([0 x],[0,y])
h = line([0 perpVec(1)],[0,perpVec(2)]);
set(h,'LineStyle','--')

x = [trials.xChangeRaw]; y = [trials.yChangeRaw];
h = line([0 x(which)],[0 y(which)])
set(h,'LineStyle','--','Color','r')

[x1,y1] = pol2cart(theta,degreesPerSaccade);
x = [trials.xChangeRaw]; y = [trials.yChangeRaw];
h = line([x1 x(which)],[y1 y(which)])
set(h,'LineStyle','-','Color','r')

ba = [trials(which).projection];
% h = line([x x+holdUnits(1,which)],[y,y+holdUnits(2,which)]);

% ylim([-10 10])
% xlim([-10 10])

title({num2str([trials(which).up]);num2str(ba)});
end
%% visualize the whole session
% 
% [trials.projection] = deal(NaN);
plotit = 1;
figure(); hold on; axis square

idx = [trials.up] == 1;

x = [trials(idx).xChangeRaw];
y = [trials(idx).yChangeRaw];

plot(x,y,'.k')

idx = [trials.up] == 0;

x = [trials(idx).xChangeRaw];
y = [trials(idx).yChangeRaw];

plot(x,y,'ok')

% fixation square
plot(0,0,'sk','MarkerFaceColor',[.5 .5 .5],'MarkerSize',15)

colors = {'r','b','g'};
tmp = unique([trials.condition]);
tmp = tmp(~isnan(tmp));
for i = 1:length(tmp)
    
    if i == 1
        selex = [trials.condition] == tmp(i);
        idx = and(selex,[trials.up] == 1);

        x = [trials(idx).xChangeRaw];
        y = [trials(idx).yChangeRaw];

        plot(x,y,'.r')

        idx = and(selex,[trials.up] == 0);

        x = [trials(idx).xChangeRaw];
        y = [trials(idx).yChangeRaw];

        plot(x,y,'or')

    end
    
    idx = [trials.condition] == tmp(i);
    x = [trials(idx).xChangeRaw];
    y = [trials(idx).yChangeRaw];
    plot(nanmean(x),nanmean(y),'s','MarkerFaceColor',colors{i});

    % plot the mean vector in that condition
    h = line([0 nanmean(x)],[0 nanmean(y)]);
    set(h,'LineStyle','--','Color',colors{i});

    % now the vector to the actual saccade position
    thetas = [trials(idx).theta];
    
    if ~MATLAB
        thetas = nanmean(thetas+180);
    else
        thetas = nanmean(thetas);
    end
    
    thetas = thetas/(180/pi);
    [x,y] = pol2cart(thetas,degreesPerSaccade);
    h = line([0 x],[0 y]);
    set(h,'Color',colors{i});
    
end

title({'solid line is actual position';'dotted line is mean vector'})

% xlim([-10 10]); ylim([-10 10])