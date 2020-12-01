% Treatment interruption tests
function treatInterruptTest(type, folder_KIM_root, coord_file, paramfile)

%close all

global shift;
global latency;
global noOfShifts;
global deepestFolder;
global plotOffline;
global KIM;
global paramData;
% Inputs
shiftHexa = 0;
latency = 0.350;


%type = 'Breath'; 
%folderHexa = 'C:\Users\dngu2802\Documents\KIM_QA\KIM QA\QA Trajectories\Validated29Oct2014\';
%fileHexa =  strcat('C:\LARK_QA_test\treat_int_new\Sitched_traces\LiverTraj_LargeSIandAPWithBreathHold_robot_20min.txt');
fileHexa = KIM.KIMRobotFile;

%folderKIM =  strcat('C:\LARK_QA_test\treat_int_new\Large_SI_AP_Breathhold');
folderKIM = KIM.KIMTrajFolder;

coordFile = KIM.KIMcoordFile
%coordFile = strcat('C:\LARK_QA_test\treat_int_new\Input_files\coord.txt');
%paramfile = strcat('C:\LARK_QA_test\treat_int_new\Input_files\param.txt');
paramfile = KIM.KIMparamFile

[~, deepestFolder] = fileparts(folderKIM);
folderKIM = [folderKIM '\'];
dirHexa = fileHexa; %[folderHexa fileHexa];
%covFile = 'covOutput_OL.txt';

% Create the file covOutput_OL.txt if it does not exist
%if exist([folderKIM covFile], 'file') == 0
 %   copyfile([folderKIM '\covOutput.txt'], [folderKIM covFile])
%end


%Obtain coordinate data
% First row x, second row y, third row z
fid = fopen(coordFile);
coordData = fscanf(fid, '%f %f %f');
fclose(fid);

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
paramData = fscanf(fid, '%f %f %f %f');
disp(paramData(1));
disp(paramData(2));
disp(paramData(3));
%paramData(1) = -400;
%paramData(2) = 0.01;
%paramData(3) = 20;
%paramData(4) = 0.1;
% List the number of files tagged with OL
listOfOLTrajFiles = ls([folderKIM '\*GA*OL*.txt']);
noOfOLTrajFiles = size(listOfOLTrajFiles,1);

listOfTrajFiles = ls([folderKIM '\*GA*.txt']);
noOfTrajFiles = size(listOfTrajFiles,1);
    

% If the trajectory files tagged OL don't exist, create them
if noOfOLTrajFiles == 0
    
    for n = 1:noOfTrajFiles
        trajFileName  = listOfTrajFiles(n,:);
        [~,trajFileNameOL,~] = fileparts(trajFileName);
        trajFileNameOL = [trajFileNameOL '_OL.txt'];
        copyfile([folderKIM '\MarkerLocationsGA_CouchShift_' num2str(n-1) '.txt'], [folderKIM trajFileNameOL])
    end
    
end

listOfOLTrajFiles = ls([folderKIM '\*GA*OL*.txt']);
noOfOLTrajFiles = size(listOfOLTrajFiles,1);

% Specify the number of trajectory files
for n = 1:noOfOLTrajFiles
    eval(['fileKIM' num2str(n) ' = ''MarkerLocationsGA_CouchShift_' num2str(n-1) '_OL.txt'';']);
    eval(['dirKIM' num2str(n) ' = [folderKIM fileKIM' num2str(n) '];']);
end



% Obtain KIM data
for n = 1:noOfOLTrajFiles
    
    eval(['fid=fopen(dirKIM' num2str(n) ');']);
    %rawDataCurrent = textscan(fid, '%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%u,%u,%f,%f,%f,%f,%f,%f', 'headerLines', 1);
    rawDataCurrent = textscan(fid,    '%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%u,%u,%f,%f,%f,%f,%f,%f, %s %*[^\n]', 'headerLines', 1);
    eval(['rawDataKIM' num2str(n) '=rawDataCurrent;']);
    fclose(fid);
end


%disp("check");
%fprintf("\n%f\t%f", [rawDataCurrent{2}, rawDataCurrent{4}]');
   

% Read the couchShifts.txt file the image directory (units in cm)
fid=fopen([folderKIM '\couchShifts.txt']);
couchPositions = textscan(fid, '%f,%f,%f\r', 'headerlines', 1);
fclose(fid);

vrt = couchPositions{1};
lng = couchPositions{2};
lat = couchPositions{3};

for n = 1:length(lat);
    if lat(n) > 950
    lat(n) = lat(n) - 1000;    
    end   
end

% In mm
%This is for varian.
shiftsAP = -diff(vrt) * 10;
shiftsSI = diff(lng) * 10; 
shiftsLR = diff(lat) * 10; 

%shiftsAP = diff(vrt) * 10;
%shiftsSI = diff(lng) * 10; 
%shiftsLR = diff(lat) * 10;
noOfShifts = length(shiftsAP);  

% Obtain the index number where shifts occur
shiftIndex = ones(noOfShifts,1);
%fprintf("\n%f\n", shiftIndex);
%disp("test");
for n = 1:noOfShifts
    
    if (n == 1)
    shiftIndex(n) = length( eval(['rawDataKIM' num2str(n) '{1}'])) + 1;
   
    else
    shiftIndex(n) = shiftIndex(n-1) + length( eval(['rawDataKIM' num2str(n) '{1}']));
    end
    
end

% Merge the KIM data from multiple files
%disp("check");
%disp(length(rawDataKIM1));
for n = 1:length(rawDataKIM1)
    for k = 1:noOfShifts + 1                
    
        if (k==1)
            rawDataKIM{n} = [eval(['rawDataKIM' num2str(k) '{n};'])];
        else
            rawDataKIM{n} = [rawDataKIM{n}; eval(['rawDataKIM' num2str(k) '{n};']) ];
        end       
    end
    
end

% K1. Obtain KIM trajectory data
% Timestamps, Correlations and Segmentations from MarkerLocations_CouchShift_0.txt
dataKIM.timestamps = rawDataKIM{2};
% Insert shifts and latency here
dataKIM.timestamps = dataKIM.timestamps  - dataKIM.timestamps(1); 
dataKIM.gantry = rawDataKIM{3};
dataKIM.index = rawDataKIM{1};
%fprintf("\n%f", length(rawDataKIM{2}));
%fprintf("\n%f", length(rawDataKIM{3}));
%fprintf("\n%f", length(rawDataKIM{4}));
%fprintf("\n%f", length(rawDataKIM{5}));
%fprintf("\n%f", length(rawDataKIM{6}));
%fprintf("\n%f", length(rawDataKIM{7}));
%fprintf("\n%f", length(rawDataKIM{8}));
%fprintf("\n%f", length(rawDataKIM{9}));
%fprintf("\n%f", length(rawDataKIM{10}));
%fprintf("\n%f", length(rawDataKIM{11}));
%fprintf("\n%f", length(rawDataKIM{12}));

d = diff(dataKIM.timestamps);
list = find(d>0.5) + 1;
dataKIM.indexOfTreatStart = list(1);
indexOfTreatStart = dataKIM.indexOfTreatStart;

% 2D
% 2D trajectories for KIM data
% Index the markers by y position where 1 is the most cranial and 3 the most caudal
% C# indexes from 0 to N-1 so a + 1 is added to each 2D trajectory for
% equivalent comparison to MATLAB
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

% 3D
% 3D trajectories for KIM data
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
%fprintf("\n%f", length(dataKIM.x1));
%fprintf("\n%f", length(dataKIM.y1));
%fprintf("\n%f", length(dataKIM.z1));
%fprintf("\n%f", length(dataKIM.x2));
%fprintf("\n%f", length(dataKIM.y2));
%fprintf("\n%f", length(dataKIM.z2));
%fprintf("\n%f", length(dataKIM.x3));
%fprintf("\n%f", length(dataKIM.y3));
%fprintf("\n%f", length(dataKIM.z3));
%fprintf("\n%f", length(dataKIM.y));
%fprintf("\n%f", length(dataKIM.z));
%fprintf("\n%f", length(rawDataKIM{7}));
%fprintf("\n%f", length(rawDataKIM{8}));
%fprintf("\n%f", length(rawDataKIM{9}));
%fprintf("\n%f", length(rawDataKIM{10}));
%fprintf("\n%f", length(rawDataKIM{11}));
%fprintf("\n%f", length(rawDataKIM{12}));



dataKIM.r1 = sqrt(dataKIM.x1.^2 + dataKIM.y1.^2 + dataKIM.z1.^2);
dataKIM.r2 = sqrt(dataKIM.x2.^2 + dataKIM.y2.^2 + dataKIM.z2.^2);
dataKIM.r3 = sqrt(dataKIM.x3.^2 + dataKIM.y3.^2 + dataKIM.z3.^2);

% D2. Compute centroid 3D trajectories for KIM data

dataKIM.xCent = (dataKIM.x1 + dataKIM.x2 + dataKIM.x3)/3 - Avg_marker_x;
dataKIM.yCent = (dataKIM.y1 + dataKIM.y2 + dataKIM.y3)/3 - Avg_marker_y;
dataKIM.zCent = (dataKIM.z1 + dataKIM.z2 + dataKIM.z3)/3 - Avg_marker_z;
dataKIM.rCent = sqrt(dataKIM.xCent.^2 + dataKIM.yCent.^2 + dataKIM.zCent.^2);

figure(1)
plot(dataKIM.timestamps, dataKIM.xCent, 'ro', dataKIM.timestamps, dataKIM.yCent, 'gs', dataKIM.timestamps, dataKIM.zCent, 'kd')

%fprintf("\n%f", length(dataKIM.xCent));
%fprintf("\n%f", length(dataKIM.yCent));
%fprintf("\n%f", length(dataKIM.zCent));
%fprintf("\n%f", length(dataKIM.timestamps));

dataKIM.xCentOff = dataKIM.xCent - dataKIM.xCent(1);
dataKIM.yCentOff = dataKIM.yCent - dataKIM.yCent(1);
dataKIM.zCentOff = dataKIM.zCent - dataKIM.zCent(1);
dataKIM.rCentOff = sqrt(dataKIM.xCentOff.^2 + dataKIM.yCentOff.^2 + dataKIM.zCentOff.^2);


% Obtain Hexa data
fid=fopen(dirHexa);
rawDataHexa = textscan(fid, '%f %f %f %f %f %f %f');
fclose(fid);
%disp(rawDataHexa{1}));
%disp(length(rawDataHexa{3}));
%disp(length(rawDataHexa{4}));


% H1. Obtain Hexamotion trajectory data
% 3D trajectories for Hexa data
dataHexa.x = (1).*rawDataHexa{2};
dataHexa.y = rawDataHexa{3};
dataHexa.z = (1).*rawDataHexa{4};
dataHexa.r = sqrt(dataHexa.x.^2 + dataHexa.y.^2 + dataHexa.z.^2);

dataHexa.xOff = dataHexa.x - dataHexa.x(1);
dataHexa.yOff = dataHexa.y - dataHexa.y(1);
dataHexa.zOff = dataHexa.z - dataHexa.z(1);
dataHexa.rOff = sqrt(dataHexa.xOff.^2 + dataHexa.yOff.^2 + dataHexa.zOff.^2);

%dataHexa.timestamps = 0:0.02:(length(dataHexa.x)-1)*0.02;
dataHexa.timestamps = rawDataHexa{1};
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

%dataHexa.timestamps = 0:0.02:(length(dataHexa.x)-1)*0.02;


% Step 1
% Plot KIM SI with couch shifts
figure(1)
plot(dataKIM.yCent, 'g.', 'linewidth', 3)
xlabel('Index', 'fontsize',16)
ylabel('SI position (mm)', 'fontsize',16)
title('Step 1: KIM with couch shifts', 'fontsize', 16)
set(gca,'fontsize',16)

% Step 2
% Plot KIM SI after undoing couch shifts
% This is so comparison can be made to Hexa 

for n = 1:noOfShifts
    dataKIM.yCent(shiftIndex(n):end) = dataKIM.yCent(shiftIndex(n):end) - shiftsSI(n);
    
    dataKIM.xCent(shiftIndex(n):end) = dataKIM.xCent(shiftIndex(n):end) - shiftsLR(n);
    dataKIM.zCent(shiftIndex(n):end) = dataKIM.zCent(shiftIndex(n):end) - shiftsAP(n);
end

figure(2)
plot(dataKIM.yCent, 'g.', 'linewidth', 3)
xlabel('Index', 'fontsize',16)
ylabel('SI position (mm)', 'fontsize',16)
title('Step 2: KIM with couch shifts undone', 'fontsize', 16)
set(gca,'fontsize',16)

% Step 3
% Find closest match between Hexa and KIM trajectory using the SI component
shift = findClosestSI(dataHexa,dataKIM);

% Apply the time shift to the KIM trajectory
totalShift = shift + latency;
dataKIM.timestamps = dataKIM.timestamps + totalShift;
%fprintf("\n%f", length(dataKIM.yCent));



% Step 4
% Plot KIM SI after time shift
figure(4)
interpHexa.y = interp1(dataHexa.timestamps, dataHexa.y, dataKIM.timestamps);
plot(dataHexa.timestamps, dataHexa.y, 'k-', dataKIM.timestamps, dataKIM.yCent, 'g.')
%plot(dataKIM.timestamps, dataKIM.yCent, 'g.')

xlabel('Index', 'fontsize',16)
ylabel('SI position (mm)', 'fontsize',16)
title('Step 4: KIM with couch shifts undone and with time shift', 'fontsize', 16)
set(gca,'fontsize',16)


% Reapply the couch shifts to KIM trajectories and apply the couch shifts to the Hexa trajectories
hexaShiftIndex = ones(noOfShifts,1);
for n = 1:length(shiftIndex)
    dataKIM.yCent(shiftIndex(n):end) = dataKIM.yCent(shiftIndex(n):end) + shiftsSI(n);   
    dataKIM.xCent(shiftIndex(n):end) = dataKIM.xCent(shiftIndex(n):end) + shiftsLR(n);
    dataKIM.zCent(shiftIndex(n):end) = dataKIM.zCent(shiftIndex(n):end) + shiftsAP(n);
    
    
    hexaShiftIndex_temp = find(abs((dataHexa.timestamps - dataKIM.timestamps(shiftIndex(n)))) < paramData(4));
%     
    if(~isempty(hexaShiftIndex_temp))
        hexaShiftIndex(n) = hexaShiftIndex_temp(1); 
   % else
%         break; 
    end
%     
    dataHexa.x(hexaShiftIndex(n):end) = dataHexa.x(hexaShiftIndex(n):end) + shiftsLR(n);
    dataHexa.y(hexaShiftIndex(n):end) = dataHexa.y(hexaShiftIndex(n):end) + shiftsSI(n);
    dataHexa.z(hexaShiftIndex(n):end) = dataHexa.z(hexaShiftIndex(n):end) + shiftsAP(n);
end

dataKIM.timeTreat = dataKIM.timestamps(indexOfTreatStart:end);
dataKIM.xCentTreat = dataKIM.xCent(indexOfTreatStart:end);
dataKIM.yCentTreat = dataKIM.yCent(indexOfTreatStart:end);
dataKIM.zCentTreat = dataKIM.zCent(indexOfTreatStart:end);




% Compute stats. Stats for couch unshifted KIM wrt unshifted Hexa
computeStats(dataHexa,dataKIM, deepestFolder, folderKIM, shift);


% Plot
plotKIMAndHexa(dataKIM,dataHexa,deepestFolder)
%findTimeGaps(dataKIM, dataHexa)

figure(1)
close
figure(2)
close
figure(3)
close
figure(4)
close


end


function shift = findClosestSI(dataHexa,dataKIM)
global paramData;
% Plot before time shift
figure(3)
%fprintf("%f\n", length(dataHexa.timestamps));
%fprintf("%f\n", length(dataHexa.y));

interpHexa.y = interp1(dataHexa.timestamps, dataHexa.y, dataKIM.timestamps);
%fprintf("%f\n", dataHexa.y);
plot(dataHexa.timestamps, dataHexa.y, 'k-', dataKIM.timestamps, dataKIM.yCent, 'g.')
xlabel('Index', 'fontsize',16)
ylabel('SI position (mm)', 'fontsize',16)
title('Step 3: KIM before time shift', 'fontsize', 16)
set(gca,'fontsize',16)

shiftValues = paramData(1):paramData(2):paramData(3);
rmseValues = ones(1,length(shiftValues));

% check that the max timestamps of dataKIM is within max time stamp of
% dataHexa, 

if (dataHexa.timestamps(end) < dataKIM.timestamps(end)) %Hexa end time smaller than KIM end time 
    % query the time in hexa time closest to KIM end time
   
    endTime_indx = find(dataKIM.timestamps < dataHexa.timestamps(end) );   
    dataKIM.timestamps = dataKIM.timestamps(1:endTime_indx(end)); 
else 
    endTime_indx = length(dataKIM.timestamps); 
end 
for n=1:length(shiftValues)
    
    interpHexa.y = interp1(dataHexa.timestamps, dataHexa.y, dataKIM.timestamps + shiftValues(n));
    rmseValues(n) = rmse(dataKIM.yCent(1:endTime_indx(end)),interpHexa.y);
    %fprintf("%f\n", interpHexa.y);
end

shift = shiftValues(rmseValues == min(rmseValues));


end

function result = rmse(x,y)
result = sum((x-y).^2)/length(x);
end


% Compute some statistics
function computeStats(dataHexa,dataKIM, deepestFolder, folderKIM, shift)
global noOfShifts
global KIM
%***************** This part writes output in a file ********************

str1 = deepestFolder
disp(deepestFolder)
str2 = '_output.txt'
file_output = strcat(folderKIM, str1, str2);
currentfolder = pwd

fileID1 = fopen(file_output, 'w');
%************************************************************************

indexOfTreatStart = dataKIM.indexOfTreatStart;

%%%%% The only evaluate up to trajectory time approach 
% if (dataHexa.timestamps(end) < dataKIM.timestamps(end)) %Hexa end time smaller than KIM end time 
%     % query the time in hexa time closest to KIM end time
%    
%     endTime_indx = find(dataKIM.timestamps < dataHexa.timestamps(end) );   
%     endTimeTreat_index = find(dataKIM.timeTreat < dataHexa.timestamps(end) );
%     dataKIM.timestamps = dataKIM.timestamps(1:endTime_indx(end)); 
%     dataKIM.timeTreat = dataKIM.timeTreat(1:endTimeTreat_index(end)); 
%     dataKIM.xCentTreat = dataKIM.xCentTreat(1:endTimeTreat_index(end)); 
%     dataKIM.yCentTreat = dataKIM.yCentTreat(1:endTimeTreat_index(end)); 
%     dataKIM.zCentTreat = dataKIM.zCentTreat(1:endTimeTreat_index(end)); 
% end 




interpHexa.x = interp1(dataHexa.timestamps, dataHexa.x, dataKIM.timestamps);
interpHexa.y = interp1(dataHexa.timestamps, dataHexa.y, dataKIM.timestamps);
interpHexa.z = interp1(dataHexa.timestamps, dataHexa.z, dataKIM.timestamps);

interpHexa.xTreat = interp1(dataHexa.timestamps, dataHexa.x, dataKIM.timeTreat);
interpHexa.yTreat = interp1(dataHexa.timestamps, dataHexa.y, dataKIM.timeTreat);
interpHexa.zTreat = interp1(dataHexa.timestamps, dataHexa.z, dataKIM.timeTreat);
%%%%% The fill in static time approach 

% rescursively fill in the NaN value with the last know value 
interpHexa.xTreat = fill_in(interpHexa.xTreat); 
interpHexa.yTreat = fill_in(interpHexa.yTreat); 
interpHexa.zTreat = fill_in(interpHexa.zTreat); 

% Get rid of offset - due to imperfect phantom alignment
%dataKIM.xCentTreat = dataKIM.xCentTreat - mean(dataKIM.xCentTreat(1:10)); 
%dataKIM.yCentTreat = dataKIM.yCentTreat - mean(dataKIM.yCentTreat(1:10)); 
%dataKIM.zCentTreat = dataKIM.zCentTreat - mean(dataKIM.zCentTreat(1:10)); 


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
 meanLR = mean(dataKIM.xCent - interpHexa.x);
 meanSI = mean(dataKIM.yCent - interpHexa.y);
 meanAP = mean(dataKIM.zCent - interpHexa.z);
 stdLR = std(dataKIM.xCent - interpHexa.x);
 stdSI = std(dataKIM.yCent - interpHexa.y);
 stdAP = std(dataKIM.zCent - interpHexa.z);

% disp('KIM Online (Pre arc included)')
% sprintf('%1.2f\n%1.2f\n%1.2f\n%1.2f\n%1.2f\n%1.2f\n', meanLR, meanSI, meanAP, stdLR, stdSI, stdAP)



%%% Check if KIM pass this dynamic test
All_mean = [meanLR, meanSI, meanAP];
All_std = [stdLR, stdSI, stdAP]; 
Failname = [' LR';' SI';' AP']; 

any_mean_fail = find(abs(All_mean) > 1); 
any_std_fail = find(abs(All_std > 2)); 
hold on

%rectangle('Position',[dataKIM.timestamps(dataKIM.indexOfTreatStart),yLow+0.01,dataKIM.timestamps(end)-dataKIM.timestamps(dataKIM.indexOfTreatStart),yUp-yLow],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.9 0.9 0.9]);
%plot(dataKIM.timestamps, dataKIM.yCent(1:index(1)) - interpHexa.y,'k-', 'linewidth', 1);

if (~isempty(any_mean_fail))
    %set(KIM.handles.text4,'string','FAIL')
    %set(KIM.handles.text4,'BackgroundColor',[1 0 0])
    line1 = ['QA result: KIM FAILED in Dyanmic test with trajectory %s', deepestFolder, ': mean difference of',]; 
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
    line1 = ['QA result: KIM FAILED in Dyanmic test with trajectory ', deepestFolder, ': standard deviation of difference of'];
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
    line1 = ['QA result: KIM PASSED in Dyanmic test with trajectory ', deepestFolder];
    fprintf(fileID1,'QA result: KIM PASSED in Dynamic test');
end 

couchLine = sprintf('No. of couch shifts %u \n', noOfShifts );

proTime1 = sprintf('Processing time per image (Online): %1.2f \n', mean(diff(dataKIM.timestamps)));

dataLine1 = sprintf('%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t', ...
    meanLRTreat, meanSITreat, meanAPTreat, stdLRTreat, stdSITreat, stdAPTreat, pctLRTreat, pctSITreat, pctAPTreat);

%%%%%% Added %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%KIM.Results.means={All_mean(1),All_mean(2),All_mean(3)};
%KIM.Results.std={All_std(1),All_std(2),All_std(3)};
%KIM.Results.pct={[num2str(round(pctLRTreat(1),2)) '/' num2str(round(pctLRTreat(2),2))],[num2str(round(pctSITreat(1),2)) '/' num2str(round(pctSITreat(2),2))],[num2str(round(pctAPTreat(1),2)) '/' num2str(round(pctAPTreat(2),2))]};
%KIM.Results.DisplayData=[KIM.Results.means;KIM.Results.std;KIM.Results.pct]
%set(KIM.handles.uitable1,'data',KIM.Results.DisplayData)

line1 = sprintf('%s \n', line1);
line2 = sprintf('Mean\t\t\tStd\t\t\tPercentile(5,95)');
line3 = sprintf('LR\tSI\tAP\tLR\tSI\tAP\tLR\t\tSI\t\tAP');

disp(line1)
disp(couchLine)
disp(proTime1)
disp(line2)
disp(line3)
disp(dataLine1)

% Write to the output file
fprintf(fileID1, '\n\n\n%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\r', ' ', 'Mean (mm)', ' ', ' ', 'Std (mm)', ' ', ' ', 'Percentile (5, 95)');
fprintf(fileID1, '%s\t %s\t %s\t %s\t %s\t %s\t %s\t\t %s\t\t %s\t\t %s\r', ' ','LR', 'SI', 'AP', 'LR', 'SI', 'AP', 'LR', 'SI', 'AP');
fprintf(fileID1, '\t%1.2f\t %1.2f\t %1.2f\t %1.2f\t %1.2f\t %1.2f\t\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\t(%1.2f, %1.2f)\r', All_mean(1), All_mean(2), All_mean(3), All_std(1), All_std(2), All_std(3), pctLRTreat, pctSITreat, pctAPTreat);
fprintf(fileID1, '\n\n\n%s\t %f\r', 'Shift is=', shift);
fclose(fileID1);
fclose('all')

end

function findTimeGaps(dataKIM, dataHexa)
global xLow;
global xUp;
global yLow;
global yUp;

timestamps = dataKIM.timestamps;
timeDiff = diff(timestamps);

length(timestamps);
length(timeDiff);

numberOfGaps = find(timeDiff > 1);
sortedTimeDiff = flipud(sort(timeDiff));
for n = 1:length(numberOfGaps) + 1
    
    if n == 1
        eval(['timestamps' num2str(n) ' = timestamps(1:numberOfGaps(n))']);
        eval(['dataKIM.xCent' num2str(n) ' = dataKIM.xCent(1:numberOfGaps(n))']);
        eval(['dataKIM.yCent' num2str(n) ' = dataKIM.yCent(1:numberOfGaps(n))']);
        eval(['dataKIM.zCent' num2str(n) ' = dataKIM.zCent(1:numberOfGaps(n))']);
    elseif n == length(numberOfGaps) + 1
        eval(['timestamps' num2str(n) ' = timestamps(numberOfGaps(n-1)+1:end)']);
        eval(['dataKIM.xCent' num2str(n) ' = dataKIM.xCent(numberOfGaps(n-1)+1:end)']);
        eval(['dataKIM.yCent' num2str(n) ' = dataKIM.yCent(numberOfGaps(n-1)+1:end)']);
        eval(['dataKIM.zCent' num2str(n) ' = dataKIM.zCent(numberOfGaps(n-1)+1:end)']);
    else
        eval(['timestamps' num2str(n) ' = timestamps(numberOfGaps(n-1)+1:numberOfGaps(n))']);
        eval(['dataKIM.xCent' num2str(n) ' = dataKIM.xCent(numberOfGaps(n-1)+1:numberOfGaps(n))']);
        eval(['dataKIM.yCent' num2str(n) ' = dataKIM.yCent(numberOfGaps(n-1)+1:numberOfGaps(n))']);
        eval(['dataKIM.zCent' num2str(n) ' = dataKIM.zCent(numberOfGaps(n-1)+1:numberOfGaps(n))']);
    end
end

figure(9);
hold on
% Plot rectangles for when KIM monitors
d = diff(timestamps);
list = find(d>1) + 1;

for n = 1:1:length(list)
    if n==length(list)
       rectangle('Position',[dataKIM.timestamps(list(n)),yLow+0.01,dataKIM.timestamps(end)-dataKIM.timestamps(list(n)),yUp-yLow],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.9 0.9 0.9])   
    else
    rectangle('Position',[dataKIM.timestamps(list(n)),yLow+0.01,dataKIM.timestamps(list(n+1)-1)-dataKIM.timestamps(list(n)),yUp-yLow],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.9 0.9 0.9])   
    end        
end

% Plot the Hexamotion traces
plot(dataHexa.timestamps, dataHexa.x,'b-', dataHexa.timestamps, dataHexa.y,'g-', dataHexa.timestamps, dataHexa.z,'r-')

% Plot the KIM trajectory for each section
for n = 1:length(numberOfGaps) + 1

plot(eval(['timestamps' num2str(n)]), eval(['dataKIM.xCent' num2str(n)]), 'b-', eval(['timestamps' num2str(n)]), eval(['dataKIM.yCent' num2str(n)]), 'g-', eval(['timestamps' num2str(n)]), eval(['dataKIM.zCent' num2str(n)]), 'r-', 'linewidth', 4)

end

% Plot the tolerance lines
Plus3 = 3*ones(length(dataHexa.timestamps),1);
Minus3 = -3*ones(length(dataHexa.timestamps),1);
plot(dataHexa.timestamps, Plus3,'k--');
plot(dataHexa.timestamps, Minus3,'k--'); hold off 
xlim([xLow xUp]);
ylim([-8 8]);

ylabel('Position (mm)', 'fontsize',16);
xlabel('Time (s)', 'fontsize',16);

global deepestFolder
title(deepestFolder, 'fontsize', 16)
set(gca,'fontsize',16)

legend('LR (Actual)', 'SI (Actual)', 'AP (Actual)', 'LR (KIM)', 'SI (KIM)', 'AP (KIM)', 'Location','NorthEastOutside' );
end

function plotKIMAndHexa(dataKIM, dataHexa,deepestFolder)
% Auto determine plot limits
global KIM
global plotOffline;
maxKIM = max(max([dataKIM.xCent dataKIM.yCent dataKIM.zCent]));
maxHexa = max(max([dataHexa.x dataHexa.y dataHexa.z]));
yUp = round(max([maxKIM maxHexa]) + 1);

minKIM = min(min([dataKIM.xCent dataKIM.yCent dataKIM.zCent]));
minHexa = min(min([dataHexa.x dataHexa.y dataHexa.z]));
yLow = round(min([minKIM minHexa]) - 1);

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

% Plot
figure(10)
%hold on

% Plot rectangles for when KIM monitors
d = diff(dataKIM.timestamps);
list = find(d>1) + 1;
%axes(KIM.handles.axes1)
hold on
%{
for n = 1:1:length(list)
    if n==length(list)
       rectangle('Position',[dataKIM.timestamps(list(n)),yLow+0.01,dataKIM.timestamps(end)-dataKIM.timestamps(list(n)),yUp-yLow],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.9 0.9 0.9])   
    else
    rectangle('Position',[dataKIM.timestamps(list(n)),yLow+0.01,dataKIM.timestamps(list(n+1)-1)-dataKIM.timestamps(list(n)),yUp-yLow],'EdgeColor',[0.8 0.8 0.8],'FaceColor',[0.9 0.9 0.9])   
    end        
end
%}


plot(dataKIM.timestamps, dataKIM.yCent,'gx', dataKIM.timestamps, dataKIM.zCent,'rx', dataKIM.timestamps, dataKIM.xCent,'bx', 'linewidth', 3)
plot(dataHexa.timestamps, dataHexa.y,'g-', dataHexa.timestamps, dataHexa.z,'r-', dataHexa.timestamps, dataHexa.x,'b-')

%xlim([0 xUp])
%ylim([-20 20])
% h = legend('SI (KIM)', 'AP (KIM)', 'LR (KIM)', 'SI (Actual)', 'AP (Actual)', 'LR (Actual)')
% M = findobj(h,'type','line');
% set(M,'linewidth',1.5)

ylabel('Position (mm)', 'fontsize',16)
xlabel('Time (s)', 'fontsize',16)

global deepestFolder
title(deepestFolder, 'fontsize', 16)
set(gca,'fontsize',16)
legend( 'SI (KIM)', 'AP (KIM)', 'LR (KIM)', 'SI (Actual)', 'AP (Actual)', 'LR (Actual)','Location','NorthEastOutside' );
%hold off

end

function NewArray = fill_in(OldArray)

    NewArray = OldArray; 
    NaN_idx = find(isnan(OldArray) == 1); 
   for i = 1:length(NaN_idx)
        idx = NaN_idx(i); 
        NewArray(idx) = NewArray(idx - 1); 
   end
end