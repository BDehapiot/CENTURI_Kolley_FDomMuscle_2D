clearvars
%%% PatScanMuscle_2D
%%% Benoit Dehapiot, PhD
%%% benoit.dehapiot@univ-amu.fr
%%% CENTURI (Turing Center for Living Systems)
%%% Multi-Engineering Platform
%%% Aix Marseille Université

%% Options
    
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

RawName = 'MAX_Trial8_D12_C1(488-TTNrb)_C2(633-MHCall)_C3(DAPI)_C4(568-Phalloidin)_100x_01_stiched';
RootPath = 'E:\3-GitHub_BDehapiot\PatScanMuscle_2D\data';
pixSize = 0.0830266; % pixel size in µm
% pixSize = 0.0706140; % pixel size in µm

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Mask
Thresh = 100; % lower threshold
ROICutOff = 0.25; % min mean int. for ROI
SizeMin = 100; % min size for a segmented object
SizeMax = 50000; % max size for a semented object

% ROIs
gSize = 10; 
Pad1 = 15; %%% Orientations
gPad1 = gSize+(gSize*(Pad1-1))*2;
Pad2 = 5; %%% Correlations 
gPad2 = gSize+(gSize*(Pad2-1))*2;
% Pad1 must be > to Pad2

% Steerable Filter
order = 2; sigma = 1;

% Filters
minProm = 0.25; % max prominence for pattern recognition
minLoc = 1; maxLoc = 4; % min/max size for pattern recognition (µm)
minLoc = minLoc/pixSize; maxLoc = maxLoc/pixSize;
    
%% Initialize

% Set Paths
RawPath = strcat(RootPath,'\',RawName,'.tif');

% MIJ Paths
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\mij.jar'.
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\ij.jar'

% Open Data
Raw = uint16(loadtiff(RawPath,1));

% Get Variable
nY = size(Raw,1);
nX = size(Raw,2);   
nGridY = floor(nY/gSize);
nGridX = floor(nX/gSize);

% Crop Raw
Raw(nGridY*gSize+1:end,:) = [];
Raw(:,nGridX*gSize+1:end) = [];
nYCrop = size(Raw,1);
nXCrop = size(Raw,2);

%% Image Processing & Mask

choice = 1;
while choice > 0   
    
% Substract Background
tempBlur = imgaussfilt(Raw,20);
Process = Raw-tempBlur;
Process = imgaussfilt(Process,1);

% Make Binary Mask
Mask = Process;
Mask(Mask<Thresh) = 0;
Mask(Mask>=Thresh) = 1;

% Make ROIsMask
ROIsMask = zeros(nGridY,nGridX);
parfor i=1:nGridY
    for j=1:nGridX
        temp = mean(Mask(gSize*i-(gSize-1):gSize*i,gSize*j-(gSize-1):gSize*j));
        if mean(temp(:)) > ROICutOff 
            if i >= Pad1 && i <= nGridY-(Pad1-1) && j >= Pad1 && j <= nGridX-(Pad1-1)
                ROIsMask(i,j) = 1;
            end
        end
    end
end

% Apply Size Filter
ROIsMask = logical(ROIsMask);
ROIsMask = bwareafilt(ROIsMask,[SizeMin SizeMax]);

% Make a Display
% tempROIsMask = imresize(ROIsMask,[nGridY*gSize nGridX*gSize],'nearest');
% tempDisplay = vertcat(im2uint8(Process), uint8(Mask)*255,uint8(tempROIsMask)*255); 
% imshow(tempDisplay);

quant = quantile(Process,0.8, 'all');
tempProcess = Process*(255/quant);
test = im2uint8(Process*255);
imshow(test);


% Make a Display
tempROIsMask = imresize(ROIsMask,[nGridY*gSize nGridX*gSize],'nearest');
tempDisplay = vertcat(uint8(Mask)*255,uint8(tempROIsMask)*255);        
MIJ.start('C:\Program Files\ImageJ_1.52t_lite\plugins');
MIJ.createImage('Process',Process,true); MIJ.run('8-bit');
MIJ.createImage('tempDisplay',tempDisplay,true);
MIJ.run('Combine...', 'stack1=Process stack2=tempDisplay combine');
clear temp tempBlur tempROIsMask tempDisplay

% Dialog Box
choice = questdlg('What next?', ...
    'Menu', ...
    'Modify Parameters','Proceed','Proceed');
switch choice
    case 'Modify Parameters'
        choice = 1;
    case 'Proceed'
        choice = 0;
end

if choice == 1
    prompt = {'Thresh :','CutOff :','SizeMin :','SizeMax :'};
    dlgtitle = 'Input'; dims = 1;
    definput = {num2str(Thresh),num2str(ROICutOff),num2str(SizeMin),num2str(SizeMax)};
    answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
    Thresh = answer(1,1); ROICutOff = answer(2,1); SizeMin = answer(3,1); SizeMax = answer(4,1);
    MIJ.run('CloseAllWindows');
    clear prompt dlgtitle dims definput answer
end

if choice == 0
    MIJ.run('CloseAllWindows');
    clear tempt prompt dlgtitle dims definput answer
end
end

% TopHat Filtering (ImageJ)
MIJ.createImage('Process',Process,true);
MIJ.run('Gray Scale Attribute Filtering', 'operation=[Top Hat] attribute=Area minimum=300 connectivity=4 operation=[Top Hat] attribute=Area minimum=300 connectivity=4');
ProcessTopHat = MIJ.getImage('Process-attrFilt');
ProcessTopHat = double(ProcessTopHat);
MIJ.run('CloseAllWindows'); MIJ.exit