clearvars
%%% PatScanMuscle_2D
%%% Benoit Dehapiot, PhD
%%% benoit.dehapiot@univ-amu.fr
%%% CENTURI (Turing Center for Living Systems)
%%% Multi-Engineering Platform
%%% Aix Marseille Université

%% Inputs

RawName = 'MAX_Trial8_D15_C1(488-TTNrb)_C2(633-MHCall)_C3(DAPI)_C4(568-Phalloidin)_100x_01_stiched';
RootPath = 'E:\3-GitHub_BDehapiot\PatScanMuscle_2D\data';
pixSize = 0.0830266; % pixel size (µm)

%% Parameters 

% Tophat filtering
TophatSize = 6; % disk size for tophat filtering (pixels)

% Mask  (obtained from BGSub)
mThresh = 100; % lower threshold for Mask (A.U.)

% ROIsMask (obtained from Mask)
rSize = 10; % ROIs size for ROIsMask (pixels)
rThresh = 0.5; % lower threshold for ROIsMask (A.U.)
rMinSize = 100; % min size for ROIsMask's objects (pixels)
rMaxSize = 50000; % max size for ROIsMask's objects (pixels)

Pad1 = 15; % Orientations
rPad1 = rSize+(rSize*(Pad1-1))*2;
Pad2 = 5; % Correlations (Pad1 must be > to Pad2)
rPad2 = rSize+(rSize*(Pad2-1))*2;

% .........................................................................

% Steerable Filter
order = 2; sigma = 10;

% .........................................................................

% Filters
minProm = 0.25; % max prominence for pattern recognition
minLoc = 1; maxLoc = 4; % min/max size for pattern recognition (µm)
minLoc = minLoc/pixSize; maxLoc = maxLoc/pixSize;
    
%% Initialize

% Set paths
RawPath = strcat(RootPath,'\',RawName,'.tif');

% Open data
Raw = uint16(loadtiff(RawPath,1));

% Get variables
nY = size(Raw,1);
nX = size(Raw,2);  
qLow = quantile(Raw,0.001,'all');
qHigh = quantile(Raw,0.999,'all');
tempScreen = get(0,'screensize');
ScreenX = tempScreen(1,3);
ScreenY = tempScreen(1,4);

%% Create Mask & ROIsMask

choice = 1;
while choice > 0   

% Process .................................................................

    % Crop Raw
    nGridY = floor(nY/rSize);
    nGridX = floor(nX/rSize);
    Raw(nGridY*rSize+1:end,:) = [];
    Raw(:,nGridX*rSize+1:end) = [];
    nYCrop = size(Raw,1);
    nXCrop = size(Raw,2);
    
    % Substract background
    BGSub = Raw-imgaussfilt(Raw,20); % Gaussian blur #1 
    BGSub = imgaussfilt(BGSub,1); % Gaussian blur #2
    
    % Tophat filtering
    Tophat = imtophat(Raw,strel('disk',TophatSize));

    % Create Mask
    Mask = BGSub;
    Mask(Mask<mThresh) = 0;
    Mask(Mask>=mThresh) = 1;

    % Create ROIsMask
    ROIsMask = zeros(nGridY,nGridX);
    parfor i=1:nGridY
        for j=1:nGridX
            temp = mean(Mask(rSize*i-(rSize-1):rSize*i,rSize*j-(rSize-1):rSize*j));
            if mean(temp(:)) > rThresh
                if i >= Pad1 && i <= nGridY-(Pad1-1) && j >= Pad1 && j <= nGridX-(Pad1-1)
                    ROIsMask(i,j) = 1;
                end
            end
        end
    end

    % Apply size filter (ROIsMask)
    ROIsMask = logical(ROIsMask);
    ROIsMask = bwareafilt(ROIsMask,[rMinSize rMaxSize]);
    
% Display .................................................................

    % BGSub
    subplot(2,1,1) 
    imshow(BGSub,[qLow qHigh])
    title('BGSub')
    
    % Tophat
    subplot(2,2,1) 
    imshow(Tophat,[qLow qHigh])
    title('Tophat')

    % Mask
    subplot(2,1,2) 
    imshow(Mask,[0 1])
    title(strcat(...
        'Mask (thresh. =',{' '},num2str(mThresh),')'));
    
    % ROIsMask
    subplot(2,2,2) 
    imshow(ROIsMask,[0 1])
    title(strcat(...
        'ROIsMask (size = ',{' '},num2str(rSize),{' '},'pix. ;',...
        {' '},'thresh =',{' '},num2str(rThresh),{' '},';',...
        {' '},'min size =',{' '},num2str(rMinSize),{' '},'pix. ;',...
        {' '},'max size =',{' '},num2str(rMaxSize),{' '},'pix.)'));
    
   	set(gcf,'Position',[20 20 ScreenX*0.8 ScreenY*0.8])
    
% Dialog box ..............................................................
    
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

        prompt = {
            'Mask thresh. :',...
            'ROIs size :',...
            'ROIs thresh. :',...
            'ROIs min. size :',...
            'ROIs max. size :'
            };

        definput = {
            num2str(mThresh),...
            num2str(rSize),....
            num2str(rThresh),...
            num2str(rMinSize),...
            num2str(rMaxSize)
            };
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        mThresh = answer(1,1); 
        rSize = answer(2,1); 
        rThresh = answer(3,1); 
        rMinSize = answer(4,1); 
        rMaxSize = answer(5,1);
        
        close
    end

    if choice == 0
        
        close
    end
    
end

%% Tophat filtering

choice = 1;
while choice > 0   

% Process .................................................................    
    
    % Tophat filtering
    Tophat = imtophat(Raw,strel('disk',TophatSize));
    
% Display .................................................................

    % Raw
    subplot(2,1,1) 
    imshow(Raw,[qLow qHigh])
    title('Raw')

    % Process
    subplot(2,1,2) 
    imshow(Tophat,[qLow qHigh])
    title(strcat(...
        'Tophat (size =',{' '},num2str(TophatSize),{' '},'pix.)'));
    
   	set(gcf,'Position',[20 20 ScreenX*0.8 ScreenY*0.8])
    
% Dialog box ..............................................................
    
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
        
        prompt = {'Tophat size :'};
        
        definput = {num2str(TophatSize)};
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        TophatSize = answer(1,1);
        
        close
    end

    if choice == 0
        
        close
    end
    
end

%% Determine local orientation

    
%     % Substract background
%     tempBlur = imgaussfilt(Raw,20);
%     Process = Raw-tempBlur;
%     Process = imgaussfilt(Process,1);
% 
%     % Tophat filtering
%     Process_tophat = imtophat(Process,strel('disk',FiltSize));

%     subplot(4,1,2) 
%     imshow(Process_tophat,[quantLow quantHigh])
%     title(strcat(...
%         'Tophat (filt. size =',{' '},num2str(FiltSize),{' '},'pix.)'));

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
