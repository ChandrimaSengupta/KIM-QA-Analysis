% Static localization tests
function staticTest_LARK(value_LR, value_SI, value_AP, ParentPath, coord)

clc
close all

global KIM;

% GO TO THE FOLDER CALLED "Static" and copy paste the path into
% "ParentPath"

% Inputs
shiftMag = 0; % In mm
frameAverage = 1;

%parentPath =  strcat('C:\LARK_QA\QA codes\Robot_traces\vrt_3_long_4_lat_5\'); 
%outputPath = 'C:\LARK_QA\QA codes\Robot_traces\vrt_3_long_4_lat_5\';
outputPath = ParentPath;

%Obtain coordinate data
% First row x, second row y, third row z

fid = fopen(coord);
coordData = fscanf(fid, '%f %f %f');
fclose(fid);

Avg_marker_x = (coordData(1) + coordData(4) + coordData(7))/3
Avg_marker_y = (coordData(2) + coordData(5) + coordData(8))/3
Avg_marker_z = (coordData(3) + coordData(6) + coordData(9))/3

Avg_marker_x_iso = 10*(Avg_marker_x - coordData(10));
Avg_marker_y_iso = 10*(Avg_marker_y - coordData(11));
Avg_marker_z_iso = 10*(Avg_marker_z - coordData(12));

Avg_marker_x = Avg_marker_x_iso;
Avg_marker_y = Avg_marker_z_iso;
Avg_marker_z = -Avg_marker_y_iso;

%Avg_marker_x = -2.2;
%Avg_marker_y = 1.7;
%Avg_marker_z= -1.6;
%Avg_marker_x = -7.2;
%Avg_marker_y = -3.4;
%Avg_marker_z = 3.2; 

data_no_displacement = readKIMData([ParentPath], frameAverage, Avg_marker_x, Avg_marker_y, Avg_marker_z);

% Compute metrics
[data1] = computeMetrics('LR', data_no_displacement, value_LR, value_SI, value_AP);


data = [data1];

metricsFile = strcat([ParentPath],'\','Metrics.txt')

fid = fopen(metricsFile, 'w');

fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\r', ' ', 'Mean (mm)', ' ', ' ', 'Std (mm)', ' ', ' ', 'Percentile (5, 95)');
fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\r', ' ','LR', 'SI', 'AP', 'LR', 'SI', 'AP', 'LR', 'SI', 'AP');

%str_header = {['Zero'],  ['+5mmLR'], ['-5mmLR'], ['+5mmSI'], ['-5mmSI'] , ['+5mmAP'], ['-5mmAP']  }; 

for i = 1:1
    fprintf(fid, '\t%1.2f\t %1.2f\t %1.2f\t %1.2f\t %1.2f\t %1.2f\t (%1.2f, %1.2f) \t(%1.2f, %1.2f) \t(%1.2f, %1.2f) \t', data(i,:));
    fprintf(fid,'\n');
end 

fclose('all')


end

function [data1] = computeMetrics(testDirection, data_no_displacement, value_LR, value_SI, value_AP)

%value_LR = 5;
%value_SI = 4;
%value_AP = -3;

initial = 10; 

data_no_displacement.xCent = data_no_displacement.xCent(initial:end); 
data_no_displacement.yCent = data_no_displacement.yCent(initial:end); 
data_no_displacement.zCent = data_no_displacement.zCent(initial:end); 

% Compute mean

meanLRZero = mean(data_no_displacement.xCent - str2double(value_LR));
meanSIZero = mean(data_no_displacement.yCent - str2double(value_SI));
meanAPZero = mean(data_no_displacement.zCent - str2double(value_AP));


% Compute stds
stdLRZero = std(data_no_displacement.xCent - str2double(value_LR));
stdSIZero = std(data_no_displacement.yCent - str2double(value_SI));
stdAPZero = std(data_no_displacement.zCent - str2double(value_AP));



%compute percentiles
pctLRZero = tsprctile((data_no_displacement.xCent - str2double(value_LR)), [5 95]);
pctSIZero = tsprctile((data_no_displacement.yCent - str2double(value_SI)), [5 95]);
pctAPZero = tsprctile((data_no_displacement.zCent - str2double(value_AP)), [5 95]);

data1 = [meanLRZero, meanSIZero, meanAPZero, stdLRZero, stdSIZero, stdAPZero, pctLRZero, pctSIZero, pctAPZero];

end


function KIMData = readKIMData(kVFolder, frameAverage, Avg_marker_x, Avg_marker_y, Avg_marker_z)
% Gets trajectory data from the trajectory file in folderPath

% Copy the trajectory and covOutput file obtained online and add the tag
% (OL)
%covFile = '\covOutput (OL).txt';
%fileKIM = '\MarkerLocationsGA_CouchShift_0.txt';

%fid = fopen([kVFolder, '\MarkerLocationsGA_CouchShift_0 (OL).txt']);

%if exist(strcat(kVFolder, fileKIM), 'file') == 0
 %   copyfile(strcat(kVFolder, '\MarkerLocationsGA_CouchShift_0.txt'), strcat(kVFolder ,fileKIM));
%end

%if exist([kVFolder covFile], 'file') == 0
%    copyfile(strcat(kVFolder ,'\covOutput.txt'), strcat(kVFolder ,covFile));
%end

% Read the file
fid=fopen([kVFolder, '\MarkerLocationsGA_CouchShift_2.txt']);
rawKIMData = textscan(fid, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f, %s %*[^\n]', 'headerLines', 1);
%rawKIMData = textscan(fid, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f', 'headerLines', 1);
fclose(fid);

KIMData.kVFrameNo = rawKIMData{1};
KIMData.timestamps = rawKIMData{2};
KIMData.timestamps = KIMData.timestamps - KIMData.timestamps(1);
KIMData.kVSourceAngle = rawKIMData{3};


if frameAverage == 1
    KIMData.kVFrameNo = KIMData.kVFrameNo + 1;
else
    KIMData.kVFrameNo = KIMData.kVFrameNo * frameAverage;
    KIMData.kVFrameNo(1) = 1;
end

% K. 3D trajectories for KIM data
% Index the markers by SI position where 1 is the most cranial and 3 the most caudal
array = [rawKIMData{6}(1) rawKIMData{9}(1) rawKIMData{12}(1)];
sortedArray = sort(array, 'descend');
indexes = [find(sortedArray(1) == array) find(sortedArray(2) == array) find(sortedArray(3) == array)];

for n = 1:3
    if indexes(n) == 1
        eval(['KIMData.x' num2str(n) '= rawKIMData{5};']);
        eval(['KIMData.y' num2str(n) '= rawKIMData{6};']);
        eval(['KIMData.z' num2str(n) '= rawKIMData{4};']);
    elseif indexes(n) == 2
        eval(['KIMData.x' num2str(n) '= rawKIMData{8};']);
        eval(['KIMData.y' num2str(n) '= rawKIMData{9};']);
        eval(['KIMData.z' num2str(n) '= rawKIMData{7};']);
    elseif indexes(n) == 3
        eval(['KIMData.x' num2str(n) '= rawKIMData{11};']);
        eval(['KIMData.y' num2str(n) '= rawKIMData{12};']);
        eval(['KIMData.z' num2str(n) '= rawKIMData{10};']);
    end
end

KIMData.r1 = sqrt(KIMData.x1.^2 + KIMData.y1.^2 + KIMData.z1.^2);
KIMData.r2 = sqrt(KIMData.x2.^2 + KIMData.y2.^2 + KIMData.z2.^2);
KIMData.r3 = sqrt(KIMData.x3.^2 + KIMData.y3.^2 + KIMData.z3.^2);

% K. Compute centroid 3D trajectories for KIM data
KIMData.xCent = (KIMData.x1 + KIMData.x2 + KIMData.x3)/3 - Avg_marker_x;
KIMData.yCent = (KIMData.y1 + KIMData.y2 + KIMData.y3)/3 - Avg_marker_y;
KIMData.zCent = (KIMData.z1 + KIMData.z2 + KIMData.z3)/3 - Avg_marker_z;
KIMData.rCent = sqrt(KIMData.xCent.^2 + KIMData.yCent.^2 + KIMData.zCent.^2);

disp(mean(KIMData.xCent))
disp(mean(KIMData.yCent))
disp(mean(KIMData.zCent))

plot(KIMData.timestamps, KIMData.xCent,'gx', KIMData.timestamps, KIMData.yCent, KIMData.timestamps, KIMData.zCent)
legend('LAT (KIM)', 'LONG (KIM)', 'VRT (KIM)');

KIMData.xCentOff = KIMData.xCent - KIMData.xCent(1);
KIMData.yCentOff = KIMData.yCent - KIMData.yCent(1);
KIMData.zCentOff = KIMData.zCent - KIMData.zCent(1);
KIMData.rCentOff = sqrt(KIMData.xCentOff.^2 + KIMData.yCentOff.^2 + KIMData.zCentOff.^2);

end





