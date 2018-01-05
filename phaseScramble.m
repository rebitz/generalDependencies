function PhaseScrable(filename)

% Phase scrambles images using random noise in the spatial frequency
% domain. Uses the same random noise in all three color channels to
% preserve some the color percept.
%
% Produces 3 different images:
%   SCR: scrambled
%   LSF: low spatial frequency preserved
%   HSF: high spatial frequency preserved
%
% 6/13/2008: created, SVS/RAE
% 6/13/2008: added masks, RAE
% 4/20/2009: added high SF preservation, RAE

% For testing:
% filename = 'test.jpg';

% Input image, convert to double
Im = imread(filename);
Im = double(Im);

% Get image size
ImSize = size(Im);

% Generate random noise map in frequency space.
RandomPhase = angle(fft2(rand(ImSize(1), ImSize(2))));
SFMask = zeros(ImSize(1), ImSize(2));
Msk1=(ImSize(1)-1);
Msk2=(ImSize(2)-1);
SFMask(1:(Msk1),1:(Msk2)) = 1;
RandomPhase = RandomPhase .* SFMask;

% Fourier transform image, dissasociate amplitude and phase, and masked
% noise to the phase domain.
ImFourier = fft2(Im);
Amp = abs(ImFourier);
Phase = angle(ImFourier);
Phase = Phase + cat(3,RandomPhase,RandomPhase,RandomPhase);

%Back Transform
ImScr = ifft2(Amp.*exp(i*Phase));

% Get rid of machine error. If there's anything weird, double check to make
% sure imaginary part is small.
ImScr = real(ImScr);


% Create a mask of the original image, use it to make sharp boundaries in
% the phase scrambled output.
% Optional - Stephen doesn't like this.
%mask_gray = Im(:,:,1);
%mask = zeros(320,320);
%index = find(mask_gray > 0);
%mask(index) = ones(size(index));
%mask = cat(3,mask,mask,mask);

%x = 1;
%y = 1;
%layer = 1;

%for layer = 1:ImSize(3)
%    for x = 1:ImSize(1)
%        for y = 1:ImSize(2)
%            Masked(x,y,layer) = ImScr(x,y,layer)*mask(x,y,layer);
%            y = y+1;
%        end
%        x = x+1;
%    end
%    layer = layer+1;
%end


% Preserve low spatial frequency space
RandomPhase = angle(fft2(rand(ImSize(1), ImSize(2))));
SFMask = zeros(ImSize(1), ImSize(2));
Msk1=(ImSize(1)-2);
Msk2=(ImSize(2)-2);
SFMask(2:(Msk1),2:(Msk2)) = 1;

RandomPhase = RandomPhase .* SFMask;
ImFourier = fft2(Im);
Amp = abs(ImFourier);
Phase = angle(ImFourier);
Phase = Phase + cat(3,RandomPhase,RandomPhase,RandomPhase);

%Back Transform
ImScrLSF = ifft2(Amp.*exp(i*Phase));

% Get rid of machine error. If there's anything weird, double check to make
% sure imaginary part is small.
ImScrLSF = real(ImScrLSF);

% Preserve high spatial frequency space
RandomPhase = angle(fft2(rand(ImSize(1), ImSize(2))));
SFMask = ones(ImSize(1), ImSize(2));
Msk1=(ImSize(1)-7);
Msk2=(ImSize(2)-7);
SFMask(7:(Msk1),7:(Msk2)) = 0;

RandomPhase = RandomPhase .* SFMask;
ImFourier = fft2(Im);
Amp = abs(ImFourier);
Phase = angle(ImFourier);
Phase = Phase + cat(3,RandomPhase,RandomPhase,RandomPhase);

%Back Transform
ImScrHSF = ifft2(Amp.*exp(i*Phase));

% Get rid of machine error. If there's anything weird, double check to make
% sure imaginary part is small.
ImScrHSF = real(ImScrHSF);


% % Display Images.
% subplot(2,2,1);
% imshow(uint8(Im));
% axis image ij;
% subplot(2,2,2);
% imshow(uint8(ImScr));
% axis image ij;
% subplot(2,2,3);
% imshow(uint8(ImScrLSF));
% axis image ij;
% subplot(2,2,4);
% imshow(uint8(ImScrHSF));
% axis image ij;


% Save output, both phase scrambled and masked images

cd scram

imwrite(uint8(ImScr),['SCR_' filename],'png');
imwrite(uint8(ImScrLSF),['LSF_' filename],'png');
imwrite(uint8(ImScrHSF),['HSF_' filename],'png');

cd ..