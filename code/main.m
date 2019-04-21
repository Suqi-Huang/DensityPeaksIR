clc; clear; close all;

%% input image
srcImgRGB = imread( '..\imageSamples\6.bmp' );

% rgb converts to gray.
[~, ~, channel] = size(srcImgRGB);
if ( channel == 3 )
    srcImg = rgb2gray(srcImgRGB);
else
    srcImg = srcImgRGB;
end
srcImg = double(srcImg);



%% Directly produce the final detection results, ignoring the intermediate process.
% if using this block of code, you can comment out the next block of code.

% detWay = DensityPeaksIR();
% [tarPos, tarCon] = detWay.finalDetect(srcImg);
   
%% Details of the detection method are presented here.
tic
detWay = DensityPeaksIR();
rhoMat = srcImg;
m = size(rhoMat, 1);
[rho, delta] = iterationElection( detWay, rhoMat );
[ classInitial ] = singularFind( detWay, rho, delta );
singularIndex = find( classInitial ~=  0 );
classCenterRows = mod( singularIndex, m );
classCenterRows(classCenterRows == 0) = m;
classCenterCols = ceil( singularIndex / m );
seedPos = [ classCenterCols, classCenterRows ];
gvr = regionGrow( detWay, rhoMat, seedPos );

% Threashold operation
confidence = confidenceCal( detWay, gvr );
posIndex = confidence > detWay.thdQuatile;
tarPos = seedPos(posIndex, :);
tarCon = confidence(posIndex);
toc

% Results show
figure
imshow(srcImgRGB);

figure
plot(rho, delta, 'LineStyle', 'none', ...
    'Color', 'k', 'Marker', '.', 'MarkerSize', 16 ); hold on;
plot(rho(singularIndex), delta(singularIndex), 'LineStyle', 'none', ...
    'Color', 'b', 'Marker', '.', 'MarkerSize', 24 );
grid on;
set(gca,'FontSize',16,'GridLineStyle',':','GridColor','k','GridAlpha',1);
xlabel('\rho','FontSize',20); ylabel('\delta','FontSize',20);
title('\rho-\delta space');

figure
imshow(srcImgRGB);
hold on;
plot(classCenterCols, classCenterRows, 'LineStyle', 'none', ...
    'LineWidth', 1.5, 'Color', 'b', 'Marker', 'o', 'MarkerSize', 8 );
hold on;
plot( tarPos(:, 1), tarPos(:, 2), 'LineStyle', 'none', ...
    'LineWidth', 1.5, 'Color', 'r', 'Marker', 'o', 'MarkerSize', 8 );
title('Candidate targets (blue) and detected targets (red)');

