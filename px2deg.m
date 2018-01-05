function degrees = px2deg(pixels,env)

% firs convert px to cm
convF = env.physicalWidth./env.width;
inputCm = pixels * convF;
rads = atan(inputCm./env.distance);
degrees = rads*(180/pi);
