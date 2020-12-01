
% Dynamic localization tests
%function dynamicTest(traj_type, folder_KIM_root)

function dynamicTest(traj_type, folder_KIM_root, coord_file, paramfile)

%close all
global deepestFolder;
global shift;
global latency;
global plotOffline
global KIM
global paramData



% ***************  Input section *****************************************
% ************************************************************************ 

%-Input Robot trajectory file

%fileHexa = strcat('LiverTraj_BreathHold1_robot', '.txt');
%fileHexa = input('Please enter the robot trajectory name within apostrophe (with the extension e.g. .txt):');
fileHexa=KIM.KIMRobotFile
%- Input KIM trajectory folder name 
%folderKIM = strcat('C:\LARK_QA\QA codes\Robot_traces', '\' , 'Transient');
folderKIM = KIM.KIMTrajFolder
%input('Please enter the folder name within apostrophe:')

coordFile = KIM.KIMcoordFile

paramfile = KIM.KIMparamFile
%- Input KIM trajectory file name 

fileKIM = '\MarkerLocationsGA_CouchShift_0.txt';
%fileKIM =input('enter file name')
fileKIMOff = '\MarkerLocationsGA_CouchShift_0.txt';


%- Input cov filename
[~, deepestFolder] = fileparts(folderKIM);
%covFile = '\covOutput (OL).txt' ;

%- Copy file
dirHexa = [fileHexa]; %[folderHexa ];
dirKIM = [folderKIM fileKIM];
dirKIMOff = [folderKIM fileKIMOff];

% Copy the trajectory and covOutput file obtained online and add the tag (OL) 
%if exist([folderKIM fileKIM], 'file') == 0
%    copyfile([folderKIM '\MarkerLocationsGA_CouchShift_0.txt'], [folderKIM fileKIM]);
%end

%if exist([folderKIM covFile], 'file') == 0
 %   copyfile([folderKIM '\covOutput.txt'], [folderKIM covFile]);
%end


%- Input parameters 

shiftHexa = 0; % time in seconds by which Hexamotion trajectory was shifted
noOfArcs = 1;

% Average of the three marker positions with respect to the isocenter
%Avg_marker_x = -7.2;
%Avg_marker_y = -3.3;
%Avg_marker_z= 3.4;

%Avg_marker_x = -2.2;
%Avg_marker_y = 1.7;
%Avg_marker_z= -1.6;


plotOffline = true;
latency = 0.2; %0.350;

% ************************************************************************
% ************************** End Input ***********************************

%*************************************************************************
%************************** Output File **********************************


%************************** Analysis *************************************
%*************************************************************************

% Obtain data from the files 
%Obtain coordinate data
% First row x, second row y, third row z
fid = fopen(coordFile);
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

% Obtain parameters 
fid = fopen(paramfile);
paramData = fscanf(fid, '%f %f %f');
disp(paramData(1));
disp(paramData(2));
disp(paramData(3));

% Obtain Hexa data
fid=fopen(dirHexa);
rawDataHexa = textscan(fid, '%f %f %f %f %f %f %f');
fclose(fid);

% Obtain KIM data
fid=fopen(dirKIM);
rawDataKIM = textscan(fid,    '%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%u,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s', 'headerLines', 1);
%celldisp(rawDataKIM(1,end));
%rawDataKIM = textscan(fid,    '%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%u,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s %*[^\n]', 'headerLines', 1);
fclose(fid);

% Obtain Offline KIM data
fid=fopen(dirKIMOff);
rawDataKIMOff = textscan(fid, '%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%u,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s %*[^\n]', 'headerLines', 1);
%rawDataKIMOff = textscan(fid, '%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%u,%u,%f,%f,%f,%f,%f,%f', 'headerLines', 1);
%fprintf("\n%f\t%f\n", rawDataKIMOff{1}, rawDataKIMOff{2});
fclose(fid);

% H. Obtain Hexamotion trajectory data
% 3D trajectories for Hexa data
dataHexa.x = (1).*rawDataHexa{2};
dataHexa.y = rawDataHexa{3};
dataHexa.z = (1).*rawDataHexa{4};
dataHexa.r = sqrt(dataHexa.x.^2 + dataHexa.y.^2 + dataHexa.z.^2);

dataHexa.xOff = dataHexa.x - dataHexa.x(1);
dataHexa.yOff = dataHexa.y - dataHexa.y(1);
dataHexa.zOff = dataHexa.z - dataHexa.z(1);
dataHexa.rOff = sqrt(dataHexa.xOff.^2 + dataHexa.yOff.^2 + dataHexa.zOff.^2);

dataHexa.timestamps = rawDataHexa{1};
%fprintf("\n%f\n",rawDataHexa{1});
% If Hexa was shifted forward, start the Hexa timestamp at the time of the
% shift
indexAtStart = shiftHexa/0.02 + 1;

dataHexa.x = dataHexa.x(indexAtStart:end);
dataHexa.y = dataHexa.y(indexAtStart:end);
dataHexa.z = dataHexa.z(indexAtStart:end);
dataHexa.r = sqrt(dataHexa.x.^2 + dataHexa.y.^2 + dataHexa.z.^2);

dataHexa.xOff = dataHexa.xOff(indexAtStart:end);
dataHexa.yOff = dataHexa.yOff(indexAtStart:end);
dataHexa.zOff = dataHexa.zOff(indexAtStart:end);
dataHexa.rOff = sqrt(dataHexa.xOff.^2 + dataHexa.yOff.^2 + dataHexa.zOff.^2);

offset = 0; 
dataHexa.timestamps = rawDataHexa{1} + offset;


% K. Obtain KIM trajectory data
% Timestamps, Correlations and Segmentations from MarkerLocations_CouchShift_0.txt
dataKIM.timestamps = rawDataKIM{2};

%disp(dataKIM.timestamps);
dataKIM.timestamps = dataKIM.timestamps - dataKIM.timestamps(1);
dataKIM.gantry = rawDataKIM{3};
dataKIM.index = rawDataKIM{1};



% KO. Obtain offline KIM trajectory data
% Timestamps, Correlations and Segmentations from MarkerLocations_CouchShift_0.txt
dataKIMOff.timestamps = rawDataKIMOff{2};
dataKIMOff.timestamps = dataKIMOff.timestamps - dataKIMOff.timestamps(1);
dataKIMOff.gantry = rawDataKIMOff{3};
dataKIMOff.index = rawDataKIMOff{1};

indexForOff = zeros(length(dataKIM.index),1);
for n = 1:length(dataKIM.index)
    
    indexForOff(n) = find(dataKIMOff.index == dataKIM.index(n));
    
end
%{
% 2D
% 2D trajectories for KIM data
% Index the markers by y position where 1 is the most cranial and 3 the most caudal
% C# indexes from 0 to N-1 so a + 1 is added to each 2D trajectory for
% equivalent KIMmparison to MATLAB
array = [rawDataKIM{14}(1) rawDataKIM{16}(1) rawDataKIM{18}(1)];
sortedArray = sort(array, 'descend');
indexes = [find(sortedArray(1) == array) find(sortedArray(2) == array) find(sortedArray(3) == array)];

for n = 1:3
    if indexes(n) == 1
        eval(['dataKIM.xp' num2str(n) '= rawDataKIM{13} + 1;']);
        eval(['dataKIM.yp' num2str(n) '= rawDataKIM{14} + 1;']);
    elseif indexes(n) == 2
        eval(['dataKIM.xp' num2str(n) '= rawDataKIM{15} + 1;']);
        eval(['dataKIM.yp' num2str(n) '= rawDataKIM{16} + 1;']);
    elseif indexes(n) == 3
        eval(['dataKIM.xp' num2str(n) '= rawDataKIM{17} + 1;']);
        eval(['dataKIM.yp' num2str(n) '= rawDataKIM{18} + 1;']);
    end
end

% 2D
% Compute centroids for the 2D coordinates
dataKIM.xpCent = (dataKIM.xp1 + dataKIM.xp2 + dataKIM.xp3) / 3 ;
dataKIM.ypCent = (dataKIM.yp1 + dataKIM.yp2 + dataKIM.yp3) / 3 ;
%}


% K. 3D trajectories for KIM data
% Index the markers by SI position where 1 is the most cranial and 3 the most caudal
array = [rawDataKIM{6}(1) rawDataKIM{9}(1) rawDataKIM{12}(1)];
sortedArray = sort(array, 'descend');
indexes = [find(sortedArray(1) == array) find(sortedArray(2) == array) find(sortedArray(3) == array)];

for n = 1:3
    if indexes(n) == 1
        eval(['dataKIM.x' num2str(n) '= rawDataKIM{5};']);
        eval(['dataKIM.y' num2str(n) '= rawDataKIM{6};']);
        eval(['dataKIM.z' num2str(n) '= rawDataKIM{4};']);
    elseif indexes(n) == 2
        eval(['dataKIM.x' num2str(n) '= rawDataKIM{8};']);
        eval(['dataKIM.y' num2str(n) '= rawDataKIM{9};']);
        eval(['dataKIM.z' num2str(n) '= rawDataKIM{7};']);
    elseif indexes(n) == 3
        eval(['dataKIM.x' num2str(n) '= rawDataKIM{11};']);
        eval(['dataKIM.y' num2str(n) '= rawDataKIM{12};']);
        eval(['dataKIM.z' num2str(n) '= rawDataKIM{10};']);
    end
end

dataKIM.r1 = sqrt(dataKIM.x1.^2 + dataKIM.y1.^2 + dataKIM.z1.^2);
dataKIM.r2 = sqrt(dataKIM.x2.^2 + dataKIM.y2.^2 + dataKIM.z2.^2);
dataKIM.r3 = sqrt(dataKIM.x3.^2 + dataKIM.y3.^2 + dataKIM.z3.^2);

% K. Compute centroid 3D trajectories for KIM data % 
% robot centroid = [-2.5885    1.3947 -0.9651 ]
%dataKIM.xCent = (dataKIM.x1 + dataKIM.x2 + dataKIM.x3)/3 + 2.6 - 0.8; %2.6;
%dataKIM.yCent = (dataKIM.y1 + dataKIM.y2 + dataKIM.y3)/3 - 2.5 + 1.467; %1.4;
%dataKIM.zCent = (dataKIM.z1 + dataKIM.z2 + dataKIM.z3)/3 + 1.0 + 0.9;
% I need x=-1.8, y=+1.0, z=-1.9 to subtract

dataKIM.xCent = (dataKIM.x1 + dataKIM.x2 + dataKIM.x3)/3 - Avg_marker_x; %2.2;
dataKIM.yCent = (dataKIM.y1 + dataKIM.y2 + dataKIM.y3)/3 - Avg_marker_y; %1.7;
dataKIM.zCent = (dataKIM.z1 + dataKIM.z2 + dataKIM.z3)/3 - Avg_marker_z; %1.6;
dataKIM.rCent = sqrt(dataKIM.xCent.^2 + dataKIM.yCent.^2 + dataKIM.zCent.^2);

dataKIM.xCentOff = dataKIM.xCent - dataKIM.xCent(1);
dataKIM.yCentOff = dataKIM.yCent - dataKIM.yCent(1);
dataKIM.zCentOff = dataKIM.zCent - dataKIM.zCent(1);
dataKIM.rCentOff = sqrt(dataKIM.xCentOff.^2 + dataKIM.yCentOff.^2 + dataKIM.zCentOff.^2);

% KO. 3D trajectories for offline KIM data
% Index the markers by SI position where 1 is the most cranial and 3 the most caudal
array = [rawDataKIMOff{6}(1) rawDataKIMOff{9}(1) rawDataKIMOff{12}(1)];
sortedArray = sort(array, 'descend');
indexes = [find(sortedArray(1) == array) find(sortedArray(2) == array) find(sortedArray(3) == array)];

for n = 1:3
    if indexes(n) == 1
        eval(['dataKIMOff.x' num2str(n) '= rawDataKIMOff{5};']);
        eval(['dataKIMOff.y' num2str(n) '= rawDataKIMOff{6};']);
        eval(['dataKIMOff.z' num2str(n) '= rawDataKIMOff{4};']);
    elseif indexes(n) == 2
        eval(['dataKIMOff.x' num2str(n) '= rawDataKIMOff{8};']);
        eval(['dataKIMOff.y' num2str(n) '= rawDataKIMOff{9};']);
        eval(['dataKIMOff.z' num2str(n) '= rawDataKIMOff{7};']);
    elseif indexes(n) == 3
        eval(['dataKIMOff.x' num2str(n) '= rawDataKIMOff{11};']);
        eval(['dataKIMOff.y' num2str(n) '= rawDataKIMOff{12};']);
        eval(['dataKIMOff.z' num2str(n) '= rawDataKIMOff{10};']);
    end
end

dataKIMOff.r1 = sqrt(dataKIMOff.x1.^2 + dataKIMOff.y1.^2 + dataKIMOff.z1.^2);
dataKIMOff.r2 = sqrt(dataKIMOff.x2.^2 + dataKIMOff.y2.^2 + dataKIMOff.z2.^2);
dataKIMOff.r3 = sqrt(dataKIMOff.x3.^2 + dataKIMOff.y3.^2 + dataKIMOff.z3.^2);

% KO. Compute centroid 3D trajectories for offline KIMOff data
dataKIMOff.xCent = (dataKIMOff.x1 + dataKIMOff.x2 + dataKIMOff.x3)/3;
dataKIMOff.yCent = (dataKIMOff.y1 + dataKIMOff.y2 + dataKIMOff.y3)/3;
dataKIMOff.zCent = (dataKIMOff.z1 + dataKIMOff.z2 + dataKIMOff.z3)/3;
dataKIMOff.rCent = sqrt(dataKIMOff.xCent.^2 + dataKIMOff.yCent.^2 + dataKIMOff.zCent.^2);

%{
% Obtain covariance data
id=fopen([folderKIM covFile]);
dataCov=textscan(id, '%f, %f, %f, %f, %f, %f, %f, %f', 'headerlines', 1);
fclose(id);

% Remove first 6 rows of dataCov
% for n = 1:length(dataCov)
%
%     dataCov{n} = dataCov{n}(10:end);
%
% end
%
% frameNo = dataCov{1};
% covYZ = dataCov{8};
%
% covData.frameNo = frameNo(1:3:end);
% covData.covYZMarker1 = covYZ(1:3:end);
% covData.covYZMarker2 = covYZ(2:3:end);
% covData.covYZMarker3 = covYZ(3:3:end);
%
% placeHolder = zeros(length(covData.frameNo),1);
% for n = 1 : length(placeHolder)
%         placeHolder(n) = find(dataKIM.index == covData.frameNo(n));
% end

%[placeHolder dataKIM.index(placeHolder) covData.frameNo]
%}


% Find closest match between Hexa and KIM trajectory using the SI component
shift = findClosestSI(dataHexa,dataKIM);

%disp(shift)

% Plotting
totalShift = shift + latency;
dataKIM.timestamps = dataKIM.timestamps + totalShift;

%covData.timestamps = dataKIM.timestamps(placeHolder);

% To find when the treatment beam was on
d= diff(dataKIM.timestamps);
[value, index] = sort(d,'descend');
indexOfTreatStart = min(index(1:noOfArcs)) + 1;
dataKIM.timestamps(indexOfTreatStart);
%dataKIM.indexOfTreatStart = 100; %indexOfTreatStart;
dataKIM.indexOfTreatStart = indexOfTreatStart;
plotKIMAndHexa(dataHexa,dataKIM, dataKIMOff, indexForOff, folderKIM, deepestFolder);


end

function shift = findClosestSI(dataHexa,dataKIM)
global paramData
global KIM
shiftValues = -30:0.01:30;
%rmseValues = ones(1,length(shiftValues));
sigma = 2; 
if (dataKIM.timestamps(end) > dataHexa.timestamps(end))
   index = find( (dataKIM.timestamps < dataHexa.timestamps(end) - sigma) & (dataKIM.timestamps > dataHexa.timestamps(end) - 2*sigma));  
   shiftValues = paramData(1):paramData(2):paramData(3);
   
else
    index = length(dataKIM.timestamps); 
    shiftValues = paramData(1):paramData(2):paramData(3);
    
end

rmseValues = ones(1,length(shiftValues));
for n=1:length(shiftValues)
    
    %interpHexa.z = interp1(dataHexa.timestamps, dataHexa.z, dataKIM.timestamps(1:index(1)) + shiftValues(n));
    interpHexa.y = interp1(dataHexa.timestamps, dataHexa.y, dataKIM.timestamps(1:index(1)) + shiftValues(n));
    %figure(2); plot(dataHexa.timestamps, dataHexa.y, 'r', dataKIM.timestamps, dataKIM.yCent, 'g', dataKIM.timestamps + shiftValues(n), interpHexa.y, 'b');
    rmseValues(n) = rmse(dataKIM.yCent(1:index(1)),interpHexa.y);
    %rmseValues(n) = rmse(dataKIM.zCent(1:index(1)),interpHexa.z);
    
end

shift = shiftValues(rmseValues == min(rmseValues));
disp(shift)
end

function result = rmse(x,y)
result = sum((x-y).^2)/length(x);
end


function plotKIMAndHexa(dataHexa,dataKIM, dataKIMOff, indexForOff, folderKIM, deepestFolder)
global KIM
global plotOffline;

% Auto determine plot limits
maxKIM = max(max([dataKIM.xCent dataKIM.yCent dataKIM.zCent]));
maxHexa = max(max([dataHexa.x dataHexa.y dataHexa.z]));
yUp = round(max([maxKIM maxHexa]) + 1);


minKIM = min(min([dataKIM.xCent dataKIM.yCent dataKIM.zCent]));
minHexa = min(min([dataHexa.x dataHexa.y dataHexa.z]));
yLow = round(min([minKIM minHexa]) - 1);

yUp = 10;
yLow = -7;

if dataKIM.timestamps(end) <60
    xUp = 80;
elseif dataKIM.timestamps(end)>60 && dataKIM.timestamps(end)<80
    xUp = 80;
elseif dataKIM.timestamps(end)>80 && dataKIM.timestamps(end)<100
    xUp = 100;
elseif dataKIM.timestamps(end)>100
    xUp = round(dataKIM.timestamps(end) + 20);
end

xLow = 0;

indexOfTreatStart = dataKIM.indexOfTreatStart;
dataKIM.timeTreat = dataKIM.timestamps(indexOfTreatStart:end);
dataKIM.xCentTreat = dataKIM.xCent(indexOfTreatStart:end);
dataKIM.yCentTreat = dataKIM.yCent(indexOfTreatStart:end);
dataKIM.zCentTreat = dataKIM.zCent(indexOfTreatStart:end);


% Plot the offline trajectory using the real-time timestamps
if (0) %plotOffline == true
    dataKIMOff.xCent = dataKIMOff.xCent(indexForOff);
    dataKIMOff.yCent = dataKIMOff.yCent(indexForOff);
    dataKIMOff.zCent = dataKIMOff.zCent(indexForOff);
    
    
    %figure(1);
    %axes(KIM.handles.axes1)
    %rectangle('Position',[dataKIM.timestamps(dataKIM.indexOfTreatStart),yLow+0.01,dataKIM.timestamps(end)-dataKIM.timestamps(dataKIM.indexOfTreatStart),yUp-yLow],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.9 0.9 0.9])
    %hold on
    % Plot the actual data
    
    %plot(dataKIM.timestamps, dataKIMOff.xCent,'bx', dataKIM.timestamps, dataKIMOff.yCent,'gx', dataKIM.timestamps, dataKIMOff.zCent,'rx')
    %plot(dataHexa.timestamps, dataHexa.x,'b-', dataHexa.timestamps, dataHexa.y,'g-', dataHexa.timestamps, dataHexa.z,'k-')
    
    
    %xlim([xLow xUp]);
    %ylim([yLow yUp]);
    
    %legend('LR (KIM)', 'SI (KIM)', 'AP (KIM)', 'LR (Actual)', 'SI (Actual)', 'AP (Actual)', 'Location','NorthEastOutside');
    %ylabel('Position (mm)', 'fontsize',16);
    %xlabel('Time (s)', 'fontsize',16);
    %title([deepestFolder ' (Offline)'], 'fontsize', 16);
    %set(gca,'fontsize',16);
    %hold off
    
end

% Plot
figure(2);
%axes(KIM.handles.axes1)
hold on

rectangle('Position',[dataKIM.timestamps(dataKIM.indexOfTreatStart),yLow+0.01,dataKIM.timestamps(end)-dataKIM.timestamps(dataKIM.indexOfTreatStart),yUp-yLow],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.9 0.9 0.9]);
plot(dataKIM.timestamps, dataKIM.yCent,'gx', dataKIM.timestamps, dataKIM.zCent,'rx', dataKIM.timestamps, dataKIM.xCent,'bx', 'linewidth', 3);
plot(dataHexa.timestamps, dataHexa.y,'g-', dataHexa.timestamps, dataHexa.z,'r-', dataHexa.timestamps, dataHexa.x,'b-');

%plot(dataKIM.timestamps, dataKIM.yCent);
%plot(dataHexa.timestamps, dataHexa.y);


%xlim([xLow xUp]);
%ylim([-20 20]);
% h = legend('SI (KIM)', 'AP (KIM)', 'LR (KIM)', 'SI (Actual)', 'AP (Actual)', 'LR (Actual)')
% M = findobj(h,'type','line');
% set(M,'linewidth',1.5)

ylabel('Position (mm)', 'fontsize',16);
xlabel('Time (s)', 'fontsize',16);
title(deepestFolder, 'fontsize', 16);
legend( 'SI (KIM)', 'AP (KIM)', 'LR (KIM)', 'SI (Actual)', 'AP (Actual)', 'LR (Actual)','Location','NorthEastOutside' );
set(gca,'fontsize',16)
hold off



computeStats(dataHexa,dataKIM,dataKIMOff, folderKIM, deepestFolder);

%figure(3);
%close
%figure(4);
%close

end

% Compute some statistics
function computeStats(dataHexa,dataKIM,dataKIMOff, folderKIM, deepestFolder)
global KIM

%***************** This part writes output in a file ********************

str1 = deepestFolder
disp(deepestFolder)
str2 = '_output.txt'
file_output = strcat(str1, str2)
%fileID1 = fopen(file_output, 'w');
%fileID1 = fopen('C:\LARK_QA\QA codes\Robot_traces\', file_output,  'w');
outputPath = strcat('C:\LARK_QA\QA codes\Robot_traces\');
%file_output = fullfile('C:\LARK_QA\QA codes\Robot_traces\', str1, str2);
currentfolder = pwd
fileID1 = fopen(file_output, 'w');
%************************************************************************


sigma = 0.5; 
if (dataKIM.timestamps(end) > dataHexa.timestamps(end))
    
    index = find( (dataKIM.timestamps < dataHexa.timestamps(end)) & (dataKIM.timestamps > dataHexa.timestamps(end) - sigma));  
else
    index = length(dataKIM.timestamps); 
end


% KO. KIM online metrics (treatment only and including pre arc)
indexOfTreatStart = dataKIM.indexOfTreatStart;

interpHexa.x = interp1(dataHexa.timestamps, dataHexa.x, dataKIM.timestamps(1:index(1)));
interpHexa.y = interp1(dataHexa.timestamps, dataHexa.y, dataKIM.timestamps(1:index(1)));
interpHexa.z = interp1(dataHexa.timestamps, dataHexa.z, dataKIM.timestamps(1:index(1)));

 interpHexa.xTreat = interp1(dataHexa.timestamps, dataHexa.x, dataKIM.timeTreat);
 interpHexa.yTreat = interp1(dataHexa.timestamps, dataHexa.y, dataKIM.timeTreat);
 interpHexa.zTreat = interp1(dataHexa.timestamps, dataHexa.z, dataKIM.timeTreat);


% Get rid of offset - due to imperfect phantom alignment
% dataKIM.xCentTreat = dataKIM.xCentTreat;  %- mean(dataKIM.xCentTreat(10:100)); 
% dataKIM.yCentTreat = dataKIM.yCentTreat; %- mean(dataKIM.yCentTreat(10:100)); 
% dataKIM.zCentTreat = dataKIM.zCentTreat; %- mean(dataKIM.zCentTreat(1:10)); 

% Get rid of offset - due to imperfect phantom alignment
%  dataKIM.xCent = dataKIM.xCent - mean(dataKIM.xCent(1:100)); 
%  dataKIM.yCent = dataKIM.yCent - mean(dataKIM.yCent(1:100)); 
%  dataKIM.zCent = dataKIM.zCent - mean(dataKIM.zCent(200:250)); 

% KO. Treatment only
meanLRTreat = mean(dataKIM.xCentTreat - interpHexa.xTreat);
meanSITreat = mean(dataKIM.yCentTreat - interpHexa.yTreat);
meanAPTreat = mean(dataKIM.zCentTreat - interpHexa.zTreat);
stdLRTreat = std(dataKIM.xCentTreat - interpHexa.xTreat);
stdSITreat = std(dataKIM.yCentTreat - interpHexa.yTreat);
stdAPTreat = std(dataKIM.zCentTreat - interpHexa.zTreat);
pctLRTreat = tsprctile((dataKIM.xCentTreat - interpHexa.xTreat), [5 95]);
pctSITreat = tsprctile((dataKIM.yCentTreat - interpHexa.yTreat), [5 95]);
pctAPTreat = tsprctile((dataKIM.zCentTreat - interpHexa.zTreat), [5 95]);

% KO. Pre arc included
meanLR = mean(dataKIM.xCent(1:index(1)) - interpHexa.x);
meanSI = mean(dataKIM.yCent(1:index(1)) - interpHexa.y);
meanAP = mean(dataKIM.zCent(1:index(1)) - interpHexa.z);
stdLR = std(dataKIM.xCent(1:index(1)) - interpHexa.x);
stdSI = std(dataKIM.yCent(1:index(1)) - interpHexa.y);
stdAP = std(dataKIM.zCent(1:index(1)) - interpHexa.z);
pctLR = tsprctile((dataKIM.xCent(1:index(1)) - interpHexa.x), [5 95]);
pctSI = tsprctile((dataKIM.yCent(1:index(1)) - interpHexa.y), [5 95]);
pctAP = tsprctile((dataKIM.zCent(1:index(1)) - interpHexa.z), [5 95]);
figure(12);

hold on

%rectangle('Position',[dataKIM.timestamps(dataKIM.indexOfTreatStart),yLow+0.01,dataKIM.timestamps(end)-dataKIM.timestamps(dataKIM.indexOfTreatStart),yUp-yLow],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.9 0.9 0.9]);
plot(dataKIM.timestamps, dataKIM.yCent(1:index(1)) - interpHexa.y,'k-', 'linewidth', 1);

disp(dataKIM.timestamps(2))
disp(length(dataKIM.timestamps))
% KOFF. KIM offline metrics (treatment only and including pre arc)
dataKIMOff.xCentTreat = dataKIMOff.xCent(indexOfTreatStart:end);
dataKIMOff.yCentTreat = dataKIMOff.yCent(indexOfTreatStart:end);
dataKIMOff.zCentTreat = dataKIMOff.zCent(indexOfTreatStart:end);

% KOFF. Treatment only
meanLRTreatOff = mean(dataKIMOff.xCentTreat - interpHexa.xTreat);
meanSITreatOff = mean(dataKIMOff.yCentTreat - interpHexa.yTreat);
meanAPTreatOff = mean(dataKIMOff.zCentTreat - interpHexa.zTreat);
stdLRTreatOff = std(dataKIMOff.xCentTreat - interpHexa.xTreat);
stdSITreatOff = std(dataKIMOff.yCentTreat - interpHexa.yTreat);
stdAPTreatOff = std(dataKIMOff.zCentTreat - interpHexa.zTreat);
pctLRTreatOff = tsprctile((dataKIMOff.xCentTreat - interpHexa.xTreat), [5 95]);
pctSITreatOff = tsprctile((dataKIMOff.yCentTreat - interpHexa.yTreat), [5 95]);
pctAPTreatOff = tsprctile((dataKIMOff.zCentTreat - interpHexa.zTreat), [5 95]);

% KOFF. Pre arc included
meanLROff = mean(dataKIMOff.xCent(1:index(1)) - interpHexa.x);
meanSIOff = mean(dataKIMOff.yCent(1:index(1)) - interpHexa.y);
meanAPOff = mean(dataKIMOff.zCent(1:index(1)) - interpHexa.z);
stdLROff = std(dataKIMOff.xCent(1:index(1)) - interpHexa.x);
stdSIOff = std(dataKIMOff.yCent(1:index(1)) - interpHexa.y);
stdAPOff = std(dataKIMOff.zCent(1:index(1)) - interpHexa.z);
pctLROff = tsprctile((dataKIMOff.xCent(1:index(1)) - interpHexa.x), [5 95]);
pctSIOff = tsprctile((dataKIMOff.yCent(1:index(1)) - interpHexa.y), [5 95]);
pctAPOff = tsprctile((dataKIMOff.zCent(1:index(1)) - interpHexa.z), [5 95]);

%figure(3);
%axes(KIM.handles.axes2)
%plot(dataKIM.xCent,'bx')
%hold on
%plot(dataKIMOff.xCent,'bo')
%plot(interpHexa.x,'b')
%hold off

% figure(4);
% plot(dataKIM.yCent,'bx')
% hold on
% plot(dataKIMOff.yCent,'bo')
% plot(interpHexa.y,'b')
% hold off



% length(dataKIM.xCent)
% length(dataKIMOff.xCent)
% length(interpHexa.x)




% disp('KIM Online (Pre arc included)')
% sprintf('%1.2f\n%1.2f\n%1.2f\n%1.2f\n%1.2f\n%1.2f\n', meanLR, meanSI, meanAP, stdLR, stdSI, stdAP)
%
% disp('KIM Offline (Pre arc included)')
% sprintf('%1.2f\n%1.2f\n%1.2f\n%1.2f\n%1.2f\n%1.2f\n', meanLROff, meanSIOff, meanAPOff, stdLROff, stdSIOff, stdAPOff)

%%% Check if KIM pass this dynamic test

if 0 %(indexOfTreatStart > (index(1) + 10))
   
    All_mean = [meanLRTreat, meanSITreat, meanAPTreat];
    All_std = [stdLRTreat, stdSITreat, stdAPTreat]; 
    
else
    
    All_mean = [meanLR, meanSI, meanAP];
    All_std = [stdLR, stdSI, stdAP];
end
Failname = [' LR';' SI';' AP']; 

any_mean_fail = find(abs(All_mean) > 1); 
any_std_fail = find(abs(All_std > 2)); 


if (~isempty(any_mean_fail))
    %set(KIM.handles.text4,'string','FAIL')
    %set(KIM.handles.text4,'BackgroundColor',[1 0 0])
    line1 = ['QA result: KIM FAILED in Dynamic test with trajectory %s', deepestFolder, ': mean difference of',]; 
    fprintf(fileID1,'\n%s\n','QA result: KIM FAILED in Dynamic test');
    for i = 1: length(any_mean_fail)
        if i == length(any_mean_fail)
            line1 = [line1, Failname(any_mean_fail(i), :), ' > or = 2 mm'];
        else
            line1 = [line1, Failname(any_mean_fail(i), :), ','];
        end 
    end
    
elseif (~isempty(any_std_fail))
    %set(KIM.handles.text4,'string','FAIL')
    %set(KIM.handles.text4,'BackgroundColor',[1 0 0])
    line1 = ['QA result: KIM FAILED in Dynamic test with trajectory ', deepestFolder, ': standard deviation of difference of'];
    fprintf(fileID1,'\n%s\n','QA result: KIM FAILED in Dynamic test');
    for i = 1: length(any_std_fail)
        if i == length(any_std_fail)
            line1 = [line1, Failname(any_std_fail(i), :), ' > or = 2 mm'];
        else
            line1 = [line1, Failname(any_std_fail(i), :), ','];
        end 
    end
else
    %set(KIM.handles.text4,'string','PASS')
    %set(KIM.handles.text4,'BackgroundColor',[0 1 0])
    line1 = ['QA result: KIM PASSED in Dynamic test with trajectory ', deepestFolder];
    fprintf(fileID1,'QA result: KIM PASSED in Dynamic test');
end 
%%% Printing results
proTime1 = sprintf('Processing time per image (Online): %1.2f seconds \n', mean(diff(dataKIM.timestamps)));
proTime2 = sprintf('Processing time per image (Offline): %1.2f seconds \n', mean(diff(dataKIMOff.timestamps)));

dataLine1 = sprintf('Online\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t', ...
    All_mean(1), All_mean(2), All_mean(3), All_std(1), All_std(2), All_std(3), pctLRTreat, pctSITreat, pctAPTreat);



%KIM.Results.means={All_mean(1),All_mean(2),All_mean(3)};
%KIM.Results.std={All_std(1),All_std(2),All_std(3)};
%KIM.Results.pct={[num2str(round(pctLRTreat(1),2)) '/' num2str(round(pctLRTreat(2),2))],[num2str(round(pctSITreat(1),2)) '/' num2str(round(pctSITreat(2),2))],[num2str(round(pctAPTreat(1),2)) '/' num2str(round(pctAPTreat(2),2))]};
%KIM.Results.DisplayData=[KIM.Results.means;KIM.Results.std;KIM.Results.pct]
%set(KIM.handles.uitable1,'data',KIM.Results.DisplayData)



dataLine2 = sprintf('Offline\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t', ...
    meanLRTreatOff, meanSITreatOff, meanAPTreatOff, stdLRTreatOff, stdSITreatOff, stdAPTreatOff, pctLRTreatOff, pctSITreatOff, pctAPTreatOff);

%line1 = sprintf('First row: Online\nSecond row: Offline');
line2 = sprintf('\tMean\t\t\tStd\t\t\tPercentile(5,95)');
line3 = sprintf('\tLR\tSI\tAP\tLR\tSI\tAP\tLR\t\tSI\t\tAP');

line1 = sprintf('%s \n', line1); 

disp(line1)
disp(proTime1)
%disp(proTime2)
disp(line2)
disp(line3)
disp(dataLine1)
%disp(dataLine2)


% Write to the output file
fprintf(fileID1, '\n\n\n%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\r', ' ', 'Mean (mm)', ' ', ' ', 'Std (mm)', ' ', ' ', 'Percentile (5, 95)');
fprintf(fileID1, '%s\t %s\t %s\t %s\t %s\t %s\t %s\t\t %s\t\t %s\t\t %s\r', ' ','LR', 'SI', 'AP', 'LR', 'SI', 'AP', 'LR', 'SI', 'AP');
fprintf(fileID1, '\t%1.2f\t %1.2f\t %1.2f\t %1.2f\t %1.2f\t %1.2f\t\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\r', All_mean(1), All_mean(2), All_mean(3), All_std(1), All_std(2), All_std(3), pctLRTreat, pctSITreat, pctAPTreat);

fclose(fileID1);
fclose('all')

% Display metrics in a sexy table

% pctLR = ['(' num2str(pctLRTreat(1)) ',' num2str(pctLRTreat(2)) ')'];
% pctSI = ['(' num2str(pctSITreat(1)) ',' num2str(pctSITreat(2)) ')'];
% pctAP = ['(' num2str(pctAPTreat(1)) ',' num2str(pctAPTreat(2)) ')'];
% 
% dataLine1 = {meanLRTreat, meanSITreat, meanAPTreat, stdLRTreat, stdSITreat, stdAPTreat, pctLR, pctSI, pctAP};
% table = figure('Position',[200 300 1000 200]);
% cnames = {'Mean|LR', 'Mean|SI', 'Mean|AP', 'Std|LR', 'Std|SI', 'Std|AP' , 'Percentile (5,95)|LR', 'Percentile (5,95)|SI', 'Percentile (5,95)|AP', 'Pass/Fail'};
% table= uitable(table,'Data',dataLine1,'ColumnName',cnames,'Position',[0 0 1000 200]);
% cWidth = 80;
% set(table,'ColumnWidth',{cWidth cWidth cWidth cWidth cWidth cWidth 120 120 120 120})

end
