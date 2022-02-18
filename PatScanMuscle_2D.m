clearvars
%%% PatScanMuscle_2D
%%% Benoit Dehapiot, PhD
%%% benoit.dehapiot@univ-amu.fr
%%% CENTURI (Turing Center for Living Systems)
%%% Multi-Engineering Platform
%%% Aix Marseille Université

%% Inputs

RawName = '29H_IFM_Nano2SLS647_MHC488_Actin561_63x_3,4x_Zstack3-1_z9_RSize';
RootPath = 'C:\Datas\3-GitHub_BDehapiot\PatScanMuscle_2D\data';
pixSize = 0.0830266; % pixel size (µm)
% pixSize = 0.0276060; % pixel size (µm)
nChannel = 3; % select channel to process

%% Parameters 

% Mask  (obtained from BGSub)
Mask_Thresh = 125; % lower threshold for Mask (A.U.)

% ROIsMask (obtained from Mask)
ROI_Size = 10; % ROIs size for ROIsMask (pixels)
ROI_Thresh = 0.25; % lower threshold for ROIsMask (A.U.)
ROI_MinSize = 100; % min size for ROIsMask's objects (pixels)
ROI_MaxSize = 50000; % max size for ROIsMask's objects (pixels)

% Tophat filtering
Tophat_Size = 6; % disk size for tophat filtering (pixels)

% Steerable Filter
Steerable_Order = 2; Steerable_Sigma = 2;

% Padding
Pad_Angle = 10; % Window size for angle measurement (nROIs)
Pad_Corr = 5; % Window size for 2D corr. measurement (nROIs)
ROI_Pad_Corr = (ROI_Size+(ROI_Size*(Pad_Corr-1))*2);

% Valid peaks
min_Prom = 0.10; % min prominence for pattern recognition
min_Loc = 0; max_Loc = 100; % min/max size for pattern recognition (µm)
    
%% Initialize

% Set paths
RawPath = strcat(RootPath,'\',RawName,'.tif');

% Open Data
info = imfinfo(RawPath);
% Raw = uint16(imread(RawPath,nChannel,'Info',info));
Raw = imread(RawPath,nChannel,'Info',info);
clear info

% Get variables
nY = size(Raw,1);
nX = size(Raw,2);  
qLow = quantile(Raw,0.001,'all');
qHigh = quantile(Raw,0.999,'all');

%% Mask & ROIsMask

choice = 1;
while choice > 0   

% Process .................................................................

    % Crop Raw
    nGridY = floor(nY/ROI_Size);
    nGridX = floor(nX/ROI_Size);
    Raw(nGridY*ROI_Size+1:end,:) = [];
    Raw(:,nGridX*ROI_Size+1:end) = [];
    nYCrop = size(Raw,1);
    nXCrop = size(Raw,2);
       
    % Substract background
    Raw_BGSub = Raw-imgaussfilt(Raw,20); % Gaussian blur #1 
    Raw_BGSub = imgaussfilt(Raw_BGSub,1); % Gaussian blur #2

    % Create Mask
    Raw_Mask = Raw_BGSub;
    Raw_Mask(Raw_Mask<Mask_Thresh) = 0;
    Raw_Mask(Raw_Mask>=Mask_Thresh) = 1;

    % Create ROIsMask
    Raw_ROIsMask = zeros(nGridY,nGridX);
    for i=1:nGridY
        for j=1:nGridX
            temp = mean(Raw_Mask(ROI_Size*i-(ROI_Size-1):ROI_Size*i,ROI_Size*j-(ROI_Size-1):ROI_Size*j));
            if mean(temp(:)) > ROI_Thresh
                if i >= Pad_Angle && i <= nGridY-(Pad_Angle-1) && j >= Pad_Angle && j <= nGridX-(Pad_Angle-1)
                    Raw_ROIsMask(i,j) = 1;
                end
            end
        end
    end

    % Apply size filter (ROIsMask)
    Raw_ROIsMask = logical(Raw_ROIsMask);
    Raw_ROIsMask = bwareafilt(Raw_ROIsMask,[ROI_MinSize ROI_MaxSize]);
    
    % Update ROI_Pad_Corr
    ROI_Pad_Corr = (ROI_Size+(ROI_Size*(Pad_Corr-1))*2);
    
% Display .................................................................

    % BGSub
    subplot(3,1,1) 
    imshow(Raw,[qLow qHigh])
    title('Raw')

    % Mask
    subplot(3,1,2) 
    imshow(Raw_Mask,[0 1])
    title(strcat(...
        'Mask (thresh. =',{' '},num2str(Mask_Thresh),')'));
    
    % ROIsMask
    subplot(3,1,3) 
    imshow(Raw_ROIsMask,[0 1])
    title(strcat(...
        'ROIsMask (ROIs size = ',{' '},num2str(ROI_Size),{' '},'pix. ;',...
        {' '},'ROIs thresh. =',{' '},num2str(ROI_Thresh),{' '},';',...         
        {' '},'min size =',{' '},num2str(ROI_MinSize),{' '},'nROIs ;',...
        {' '},'max size =',{' '},num2str(ROI_MaxSize),{' '},'nROIs)'));
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
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
            'Mask_Thresh :',...
            'ROI_Size :',...
            'ROI_Thresh :',...
            'ROI_MinSize :',...
            'ROI_MaxSize :'
            };

        definput = {
            num2str(Mask_Thresh),...
            num2str(ROI_Size),....
            num2str(ROI_Thresh),...
            num2str(ROI_MinSize),...
            num2str(ROI_MaxSize)
            };
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        Mask_Thresh = answer(1,1); 
        ROI_Size = answer(2,1); 
        ROI_Thresh = answer(3,1); 
        ROI_MinSize = answer(4,1); 
        ROI_MaxSize = answer(5,1);
        
        close
    end

    if choice == 0
        
        close
    end
    
end

%% Tophat

choice = 1;
while choice > 0   

% Process .................................................................    
    
    % Tophat filtering
    Raw_Tophat = imtophat(Raw,strel('disk',Tophat_Size));
    
% Display .................................................................

    % Raw
    subplot(2,1,1) 
    imshow(Raw,[qLow qHigh])
    title('Raw')

    % Tophat
    subplot(2,1,2) 
    imshow(Raw_Tophat,[qLow qHigh])
    title(strcat(...
        'Tophat (size =',{' '},num2str(Tophat_Size),{' '},'pix.)'));
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
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
        
        prompt = {'Tophat_Size :'};
        
        definput = {num2str(Tophat_Size)};
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        Tophat_Size = answer(1,1);
        
        close
    end

    if choice == 0
        
        close
    end
    
end

%% AngleMap

choice = 1;
while choice > 0   

% Process .................................................................

    % Steerable Filter
    RawRSize = double(imresize(Raw,[nGridY nGridX],'nearest'));
    MaskRSize = double(imresize(Raw_Mask,[nGridY nGridX],'nearest'));
    [res,~,~,rot] = steerableDetector(RawRSize,Steerable_Order,Steerable_Sigma,180);
    for i=1:size(rot,3)
        temp = rot(:,:,i);
        temp(MaskRSize==0) = NaN;
        rot(:,:,i) = temp;
    end

    % Make AngleMap
    AngleMap = zeros(nGridY,nGridX);
    parfor i=1:nGridY
        for j=1:nGridX
            if Raw_ROIsMask(i,j) == 1            
                Crop = rot(i-(Pad_Angle-1):i+(Pad_Angle-1),j-(Pad_Angle-1):j+(Pad_Angle-1),:);
                idxMax = NaN(size(Crop,3),1);
                for k=1:size(Crop,3)
                    temp = Crop(:,:,k);
                    idxMax(k,1) = nanmean(temp(:));
                end
                [M,I] = max(idxMax);
                AngleMap(i,j) = I;            
            end
        end
    end
    
% Display .................................................................

    % Raw
    subplot(3,1,1)
    imshow(Raw,[qLow qHigh])
    title('Raw')

    % Steerable filter
    subplot(3,1,2)
    imshow(res,[min(res(:)) max(res(:))])
    title(strcat(...
        'Steerable filter (sigma =',{' '},num2str(Steerable_Sigma),{' '},'pix.)'));

    % AngleMap
    subplot(3,1,3)
    imshow(AngleMap,[0 180])
    title('AngleMap')
    colormap(gca, jet(256));
    colorbar(gca);

    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
   
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
            'Steerable_Sigma :'
            };

        definput = {
            num2str(Steerable_Sigma)
            };
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        Steerable_Sigma = answer(1,1); 
        
        close
    end

    if choice == 0
        
        close
    end    
    
end

%% 2D cross-correlation & findpeaks

MergedData = cell(nGridY,nGridX);
parfor i=1:nGridY
    for j=1:nGridX
        
        if Raw_ROIsMask(i,j) == 1  

            % Crop Data
            Angle = AngleMap(i,j);
            Crop = Raw_Tophat(...
                ROI_Size*i-((ROI_Size*Pad_Corr)-1):...
                ROI_Size*i+(ROI_Size*(Pad_Corr-1)),...
                ROI_Size*j-((ROI_Size*Pad_Corr)-1):...
                ROI_Size*j+(ROI_Size*(Pad_Corr-1)));

            % Compute Corr2D
            Corr2D = normxcorr2(Crop,Crop); 

            % Rotate Corr2D according to "Angle" 
            Corr2D_Rot = double(imrotate(Corr2D,(90-Angle)*-1,'bilinear')); 

            % Average value on x axis 
            Corr2D_Rot_AvgX = nanmean(Corr2D_Rot(...
                round(size(Corr2D_Rot,1)/2)-10:round(size(Corr2D_Rot,1)/2)+10,:),1)'; 

            % Find max correlation 
            Max = max(Corr2D_Rot_AvgX(...
                round(size(Corr2D_Rot,1)/2)-10:round(size(Corr2D_Rot,1)/2)+10));
            idxMax = find(Corr2D_Rot_AvgX==Max);
            
            % Max normalization
            Corr2D_Rot_AvgX = Corr2D_Rot_AvgX/Max;

            % Get "one-sided" correlation profile
            rProfile = Corr2D_Rot_AvgX(idxMax:idxMax+ROI_Pad_Corr-1); % right corr. profile
            lProfile = flip(Corr2D_Rot_AvgX(idxMax-ROI_Pad_Corr+1:idxMax)); % left corr. profile              
            mProfile = horzcat(rProfile,lProfile); % merged corr. profile  
            CorrProfile = mean(mProfile,2); 
            
            % Find Peaks
            [pks,locs,width,prom] = findpeaks(CorrProfile,'MinPeakDistance',6);
            if ~isempty(pks)
                [~,idx] = max(prom); % Get idx for the most prominent peak
                PeaksInfo = horzcat(i,j,pks(idx,1),locs(idx,1)*pixSize,width(idx,1),prom(idx,1));
            else
                PeaksInfo = NaN(1,6); % NaN if no peak found
            end

            % Append MergedData
            MergedData{i,j}{1,1} = Angle;
            MergedData{i,j}{2,1} = Crop;
            MergedData{i,j}{3,1} = CorrProfile; 
            MergedData{i,j}{4,1} = PeaksInfo; 
            
        end
    end
end

%% Get results

choice = 2;
while choice == 2   

% Process .................................................................

    % Merge "All" data
    PeaksInfo_All = NaN(0,0);
    CorrProfile_All = NaN(0,0);
    Map_Prom = zeros(nGridY,nGridX);
    for i=1:nGridY
        for j=1:nGridX
            if Raw_ROIsMask(i,j) == 1
                PeaksInfo_All = vertcat(PeaksInfo_All,MergedData{i,j}{4,1});
                CorrProfile_All(:,end+1) = MergedData{i,j}{3,1};
                Map_Prom(i,j) = MergedData{i,j}{4,1}(1,6);
            end
        end
    end

    % Get "Valid" data          
    PeaksInfo_Valid = NaN(0,0);
    CorrProfile_Valid = NaN(0,0);
    Map_Valid = zeros(nGridY,nGridX);
    for i=1:size(PeaksInfo_All,1)
        tempY = PeaksInfo_All(i,1);
        tempX = PeaksInfo_All(i,2);
        tempProm = PeaksInfo_All(i,6);
        tempLoc = PeaksInfo_All(i,4);
        if tempProm >= min_Prom
            if tempLoc >= min_Loc && tempLoc <= max_Loc
                PeaksInfo_Valid(end+1,:) = PeaksInfo_All(i,:);
                CorrProfile_Valid(:,end+1) = CorrProfile_All(:,i);
                Map_Valid(tempY,tempX) = 1;
            end
        end
    end

    % Get Avg CorrProfile
    CorrProfile_Avg(:,1) = (0:ROI_Pad_Corr-1)*pixSize; % distance (µm)
    CorrProfile_Avg(:,2) = mean(CorrProfile_All,2); % All 
    CorrProfile_Avg(:,3) = std(CorrProfile_All,0,2); % All S.D.
    if ~isempty(CorrProfile_Valid)
        CorrProfile_Avg(:,4) = mean(CorrProfile_Valid,2); % Valid
        CorrProfile_Avg(:,5) = std(CorrProfile_Valid,0,2); % Valid S.D.
    else
        CorrProfile_Avg(:,4) = NaN;
        CorrProfile_Avg(:,5) = NaN;
    end
    
% Display .................................................................

    % Raw
    subplot(3,3,1:3)
    tempRaw = Raw;
    tempValidMap = imresize(Map_Valid,[nYCrop nXCrop],'nearest');
    tempValidMap = bwmorph(logical(tempValidMap),'remove');
    tempValidMap = bwmorph(tempValidMap,'dilate');
    tempRaw(tempValidMap==1) = max(Raw(:));
    imshow(tempRaw,[qLow qHigh])
	title(strcat(...
        'MapValid (min. prominence =',{' '},num2str(min_Prom),{' '},';',...
        {' '},'min. localization =',{' '},num2str(min_Loc),{' '},';',...
        {' '},'max. localization =',{' '},num2str(max_Loc),{' '},')'));

    % PromMap
    subplot(3,3,4:6)
    imshow(imresize(Map_Prom,[nYCrop nXCrop],'nearest'),[min(Map_Prom(:)) max(Map_Prom(:))])
    title('MapProm')
    colormap(gca, jet(256));
    colorbar(gca);

    % CorrProfile_AvgX (All & Valid)
    subplot(3,3,7)
    plot(CorrProfile_Avg(:,1),CorrProfile_Avg(:,2))
    hold on
    plot(CorrProfile_Avg(:,1),CorrProfile_Avg(:,4))
    xlabel('Distance (µm)') 
    ylabel('Correlation')
    legend('All','Valid')
    xlim([0 ceil(max(CorrProfile_Avg(:,1)))])
    pbaspect([1 1 1])
    title('Avg. Corr2D')

    % Locs distribution ('All') 
    subplot(3,3,8)
    histogram(PeaksInfo_All(:,4),50)
    xlabel('Distance (µm)') 
    ylabel('Number of peaks')
    xlim([0 ceil(max(CorrProfile_Avg(:,1)))])
    pbaspect([1 1 1])
    title('Main peak localization (All)')

    % Locs distribution ('Valid')  
    subplot(3,3,9)
    histogram(PeaksInfo_Valid(:,4),50)
    xlabel('Distance (µm)') 
    ylabel('Number of peaks')
    xlim([0 ceil(max(CorrProfile_Avg(:,1)))])
    pbaspect([1 1 1])
    title('Main peak localization (Valid)')

    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
% Dialog box ..............................................................
    
    choice = questdlg('What next?', ...
        'Menu', ...
        'Modify Parameters','Terminate','Terminate & Save','Terminate & Save');
    switch choice
        case 'Modify Parameters'
            choice = 2;
        case 'Terminate'
            choice = 1;
        case 'Terminate & Save'
            choice = 0;
    end

    if choice == 2

        prompt = {
            'min_Prom :',...
            'min_Loc :',...
            'max_Loc :'
            };

        definput = {
            num2str(min_Prom),...
            num2str(min_Loc),....
            num2str(max_Loc)
            };
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        min_Prom = answer(1,1); 
        min_Loc = answer(2,1); 
        max_Loc = answer(3,1); 
        
        close
    end

    if choice < 2
        
        close
    end
    
end

%% Format and save

% Format outputs ..........................................................

% Parameters
Parameters = vertcat(...
    Mask_Thresh,...
    ROI_Size,ROI_Thresh,ROI_MinSize,ROI_MaxSize,...
    Tophat_Size,...
    Steerable_Order,Steerable_Sigma,...
    Pad_Angle, Pad_Corr, ROI_Pad_Corr,...
    min_Prom, min_Loc, max_Loc);

Parameters = array2table(Parameters,'RowNames',...
    {'Mask_Thresh',...
    'ROI_Size','ROI_Thresh','ROI_MinSize','ROI_MaxSize',...
    'Tophat_Size',...
    'Steerable_Order','Steerable_Sigma',...
    'Pad_Angle', 'Pad_Corr', 'ROI_Pad_Corr',...
    'min_Prom', 'min_Loc', 'max_Loc'});

PeaksInfo_All = array2table(PeaksInfo_All,'VariableNames',...
    {'GridY','GridX','pks','loc','width','prom'});
PeaksInfo_Valid = array2table(PeaksInfo_Valid,'VariableNames',...
    {'GridY','GridX','pks','loc','width','prom'});
CorrProfile_Avg = array2table(CorrProfile_Avg,'VariableNames',...
    {'Distance','All','Valid','AllSD','ValidSD'});

% Clearvars ...............................................................

clearvars -except...
    RootPath RawName choice...
    Raw Raw_BGSub Raw_Mask Raw_ROIsMask Raw_Tophat...
    MergedData... 
    Parameters...
    CorrProfile_All PeaksInfo_All...
    CorrProfile_Valid PeaksInfo_Valid...
    CorrProfile_Avg...
    Map_Prom Map_Valid

% Save ....................................................................
if choice == 0
    save(strcat(RootPath,'/temp/',RawName,'.mat'))
end
clear choice
