clc; clear; close all;

% Import files for edge detection
folderPath = '1.0s_Exposure_S1_2/1.0s_Exposure/';
files = dir(folderPath);
fileNames = {files(~[files.isdir]).name};
fileNames = string(fullfile(folderPath, fileNames));


%Initialize for center and edge
background = rgb2gray(imread(fileNames(end)));
maxGFPImg = rgb2gray(imread(fileNames(end-1))-background);
maxGFP = max(maxGFPImg, [], 'all');
disp('Click on the center of the plate.');
figure;
imshow(imread(fileNames(end)));
[x_center, y_center] = ginput(1);

disp('Click on the edge of the plate.');
[x_edge, y_edge] = ginput(1);

plateRadius = sqrt((x_center-x_edge)^2 + (y_center-y_edge)^2);



radiusList = [];
fileNameList = fileNames.';
for i = 1:size(fileNames, 2)-1
    img = imread(fileNames(i));
    gray = rgb2gray(img-background);
    maxThresh = max(gray(:, :));

    figure;
    hold on

    %Pick edge that is above 20% of max threshold at each timepoint
    edge = (img-background) >= 0.2*maxThresh;
    radius = sqrt((x_center-xGFP)^2 + (y_center-yGFP)^2);
    radiusList = [radiusList; radius];

    viscircles([x_center, y_center], radius);
end