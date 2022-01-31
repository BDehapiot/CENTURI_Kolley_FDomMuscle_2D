clearvars
%%% PatScanMuscle_2D
%%% Benoit Dehapiot, PhD
%%% benoit.dehapiot@univ-amu.fr
%%% CENTURI (Turing Center for Living Systems)
%%% Multi-Engineering Platform
%%% Aix Marseille Université

%% Options
    
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

RawName = 'MAX_Trial8_D15_C1(488-TTNrb)_C2(633-MHCall)_C3(DAPI)_C4(568-Phalloidin)_100x_01_stiched';
RootPath = 'C:\Datas\3-GitHub_BDehapiot\PatScanMuscle_2D\data';

% .........................................................................

pixSize = 0.0830266; % pixel size (µm)

% .........................................................................

ScreenX = 2560; % screen horz. size (pixels)
ScreenY = 1440; % screen vert. size (pixels)

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Tophat 
FiltSize = 6; % disk size for tophat filtering (pixels)

% .........................................................................

% Mask
Thresh = 100; % lower threshold (A.U.)
SizeMin = 100; % min size for a segmented object (pixels)
SizeMax = 50000; % max size for a semented object (pixels)

% .........................................................................

% ROIs
gSize = 10; % ROIs size (pixels)
gThresh = 0.5; % ROIs thresh (mean of Mask)

% Orientations
Pad1 = 15; 
gPad1 = gSize+(gSize*(Pad1-1))*2;

% Correlations (Pad1 must be > to Pad2)
Pad2 = 5; 
gPad2 = gSize+(gSize*(Pad2-1))*2;

% .........................................................................

% Steerable Filter
order = 2; sigma = 10;

% .........................................................................

% Filters
minProm = 0.25; % max prominence for pattern recognition
minLoc = 1; maxLoc = 4; % min/max size for pattern recognition (µm)
minLoc = minLoc/pixSize; maxLoc = maxLoc/pixSize;
    
%% Initialize

% Set Paths
RawPath = strcat(RootPath,'\',RawName,'.tif');

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

%% Image processing & mask

choice = 1;
while choice > 0   
    
    % Substract background
    tempBlur = imgaussfilt(Raw,20);
    Process = Raw-tempBlur;
    Process = imgaussfilt(Process,1);

    % Tophat filtering
    Process_tophat = imtophat(Process,strel('disk',FiltSize));

    % Make binary mask
    Mask = Process;
    Mask(Mask<Thresh) = 0;
    Mask(Mask>=Thresh) = 1;

    % Make ROIsMask
    ROIsMask = zeros(nGridY,nGridX);
    parfor i=1:nGridY
        for j=1:nGridX
            temp = mean(Mask(gSize*i-(gSize-1):gSize*i,gSize*j-(gSize-1):gSize*j));
            if mean(temp(:)) > gActiv
                if i >= Pad1 && i <= nGridY-(Pad1-1) && j >= Pad1 && j <= nGridX-(Pad1-1)
                    ROIsMask(i,j) = 1;
                end
            end
        end
    end

    % Apply size filter
    ROIsMask = logical(ROIsMask);
    ROIsMask = bwareafilt(ROIsMask,[SizeMin SizeMax]);
    
    % Make a display %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 
    quantLow = quantile(Process,0.001, 'all');
    quantHigh = quantile(Process,0.999, 'all');
    
    % Plot
    subplot(4,1,1) 
    imshow(Process,[quantLow quantHigh])
    title('Raw')
    
    subplot(4,1,2) 
    imshow(Process_tophat,[quantLow quantHigh])
    title(strcat(...
        'Tophat (filt. size =',{' '},num2str(FiltSize),{' '},'pix.)'));
    
    subplot(4,1,3) 
    imshow(Mask,[0 1])
    title(strcat(...
        'Mask (thresh. =',{' '},num2str(Thresh),{' '},'A.U.)'));
    
    subplot(4,1,4) 
    imshow(ROIsMask,[0 1])
    title(strcat(...
        'ROIsMask (size = ',{' '},num2str(gSize),{' '},'pix. ;',...
        {' '},'activ =',{' '},num2str(gActiv),{' '},'%)'));
    
   	set(gcf,'Position',[100 100 900 900])
    
%     set(gcf,'Position',[0 0 FiberSize/2 900])
%     title(strcat('Day#',num2str(SelectDay),' Cat#',num2str(SelectCat),' Fiber#',num2str(randFiber(i,1))))
%     
    
%     tempROIsMask = imresize(ROIsMask,[nGridY*gSize nGridX*gSize],'nearest');
%     tempDisplay = vertcat(Process, Process_tophat, uint16(Mask)*quantHigh, uint16(tempROIsMask)*quantHigh); 
%     imshow(tempDisplay,[quantLow quantHigh]);

    % Dialog box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        prompt = {'FiltSize :','Thresh :','SizeMin :','SizeMax :'};
        dlgtitle = 'Input'; dims = 1;
        definput = {num2str(FiltSize),num2str(Thresh),num2str(SizeMin),num2str(SizeMax)};
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        FiltSize = answer(1,1); Thresh = answer(2,1); SizeMin = answer(3,1); SizeMax = answer(4,1);
        
        close
    end

    if choice == 0
        
        close
    end
    
end

%% Determine local orientation

% choice = 1;
% while choice > 0   
% 
%     % Steerable Filter
%     RawRSize = double(imresize(Raw,[nGridY nGridX],'nearest'));
%     MaskRSize = double(imresize(Mask,[nGridY nGridX],'nearest'));
%     [res,~,~,rotations] = steerableDetector(RawRSize,order,sigma,180);
%     for i=1:size(rotations,3)
%         temp = rotations(:,:,i);
%         temp(MaskRSize==0) = NaN;
%         rotations(:,:,i) = temp;
%     end
% 
%     %     % Make a Display
%     %     MIJ.start('C:\Program Files\ImageJ_1.52t_lite\plugins');
%     %     MIJ.createImage('rotations',rotations,true);  
% 
%     % Make AngleMap
%     AngleMap = zeros(nGridY,nGridX);
%     parfor i=1:nGridY
%         for j=1:nGridX
%             if ROIsMask(i,j) == 1  
% 
%                 Crop = rotations(i-(Pad1-1):i+(Pad1-1),j-(Pad1-1):j+(Pad1-1),:);  
%                 tempMax = NaN(size(Crop,3),1);
%                 for k=1:size(Crop,3)
%                     temp = Crop(:,:,k);
%                     tempMax(k,1) = nanmean(temp(:));
%                 end
%                 [M,I] = max(tempMax);
%                 AngleMap(i,j) = I;
% 
%             end
%         end
%     end
% 
%     % Make a Display
%     quantLow = quantile(AngleMap(AngleMap~=0),0.05, 'all');
%     quantHigh = quantile(AngleMap(AngleMap~=0),0.95, 'all');
%     tempAngleMap = imresize(AngleMap,[nGridY*gSize nGridX*gSize],'nearest');
%     tempDisplay = vertcat(tempAngleMap); 
%     imshow(tempDisplay,[quantLow quantHigh], 'Colormap', jet);
%     colorbar
%     
% %     MIJ.start('C:\Program Files\ImageJ_1.52t_lite\plugins');
% %     MIJ.createImage('AngleMap',AngleMap,true); MIJ.run('Fire');  
% 
%     % Dialog Box
%     choice = questdlg('What next?', ...
%         'Menu', ...
%         'Modify Parameters','Proceed','Proceed');
%     switch choice
%         case 'Modify Parameters'
%             choice = 1;
%         case 'Proceed'
%             choice = 0;
%     end
% 
%     if choice == 1
%         prompt = {'sigma :'};
%         dlgtitle = 'Input'; dims = 1;
%         definput = {num2str(sigma)};
%         answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
%         sigma = answer(1,1);
%         MIJ.run('CloseAllWindows');
%         clear prompt dlgtitle dims definput answer
%     end
% 
%     if choice == 0
%         MIJ.run('CloseAllWindows');
%         clear tempt prompt dlgtitle dims definput answer
%     end
%     
% end
