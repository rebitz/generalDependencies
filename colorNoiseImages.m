nColors = 32;
imSize = [100,100];

orientations = [1:360/nColors:360]*(pi/180);

cmx = pickColors(nColors);

% we'll have to rescale to make the actual movie:
% cmx = cmx.*255;

% make a pretend movie for trying stuff out
nFrames = 50;
% and keep the 
keepW = NaN(nFrames,nColors);

% weights
% uniform:
w = 1/nColors;
unif_w = repmat(w,1,nColors);

for frame = 1:nFrames
    % biased
%     vonmisespdf = @(x,mu,K) exp(K*cos(x-mu)) ./ (2*pi * besseli(0,K));
    vonmises = @(x,mu,K) exp((-1/2)*(((x-mu)/K).^2));
    K = pi/4;
    mu = orientations(randi(nColors,1));
    bias_w = vonmises(orientations,mu,K);
    bias_w = bias_w ./ sum(bias_w);

    % illustrate this kernel
    % figure(); plot(orientations,w)
    
    % combine the biased and uniform in some proportion
    bias = (rand/2)+0.5; % keep the bias pretty high - this term makes all the stimuli more correlated with eachother
    w = bias*bias_w + (1-bias)*unif_w;
    keepW(frame,:) = w;

    % get index to the color we want to put in there
    imIndx = randsample(nColors,prod(imSize),true,w);
    im = NaN(imSize(1),imSize(2),3);
    for i = 1:3 % rgb
        im(:,:,i) = reshape(cmx(imIndx,i),imSize);
    end

    figure(99);
    imagesc(im); axis square
    
    pause(0.2);
end

% xlim([0 imSize(1)]);
% ylim([0 imSize(2)]);

%% 

% mean representation of each color ch
figure();
bar(nanmean(keepW))

% time series of same:
figure();
plot(keepW)