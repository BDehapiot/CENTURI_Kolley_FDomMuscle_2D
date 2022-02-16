clearvars
%%% PatScanMuscle_2D
%%% Benoit Dehapiot, PhD
%%% benoit.dehapiot@univ-amu.fr
%%% CENTURI (Turing Center for Living Systems)
%%% Multi-Engineering Platform
%%% Aix Marseille Université

%% Inputs

RawName = 'MAX_Trial8_D12_C1(488-TTNrb)_C2(633-MHCall)_C3(DAPI)_C4(568-Phalloidin)_100x_01_stiched_med';
RootPath = 'C:\Datas\3-GitHub_BDehapiot\PatScanMuscle_2D\data';
pixSize = 0.0830266; % pixel size (µm)
nChannel = 1; % select channel to process

%% Parameters 

% Tophat filtering
TophatSize = 6; % disk size for tophat filtering (pixels)

% Mask  (obtained from BGSub)
mThresh = 125; % lower threshold for Mask (A.U.)

% ROIsMask (obtained from Mask)
rSize = 10; % ROIs size for ROIsMask (pixels)
rThresh = 0.25; % lower threshold for ROIsMask (A.U.)
rMinSize = 100; % min size for ROIsMask's objects (pixels)
rMaxSize = 50000; % max size for ROIsMask's objects (pixels)

Pad1 = 10; % Orientations
rPad1 = rSize+(rSize*(Pad1-1))*2;
Pad2 = 5; % Correlations (Pad1 must be > to Pad2)
rPad2 = rSize+(rSize*(Pad2-1))*2;

% Steerable Filter
sfOrder = 2; sfSigma = 3;

% Valid peaks
minProm = 0.10; % min prominence for pattern recognition
minLoc = 0; maxLoc = 100; % min/max size for pattern recognition (µm)
    
%% Initialize

% Set paths
RawPath = strcat(RootPath,'\',RawName,'.tif');

% Open Data
info = imfinfo(RawPath);
Raw = uint16(imread(RawPath,nChannel,'Info',info));
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
    nGridY = floor(nY/rSize);
    nGridX = floor(nX/rSize);
    Raw(nGridY*rSize+1:end,:) = [];
    Raw(:,nGridX*rSize+1:end) = [];
    nYCrop = size(Raw,1);
    nXCrop = size(Raw,2);
    
    % Substract background
    BGSub = Raw-imgaussfilt(Raw,20); % Gaussian blur #1 
    BGSub = imgaussfilt(BGSub,1); % Gaussian blur #2

    % Create Mask
    Mask = BGSub;
    Mask(Mask<mThresh) = 0;
    Mask(Mask>=mThresh) = 1;

    % Create ROIsMask
    ROIsMask = zeros(nGridY,nGridX);
    for i=1:nGridY
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
    subplot(3,1,1) 
    imshow(BGSub,[qLow qHigh])
    title('BGSub')

    % Mask
    subplot(3,1,2) 
    imshow(Mask,[0 1])
    title(strcat(...
        'Mask (thresh. =',{' '},num2str(mThresh),')'));
    
    % ROIsMask
    subplot(3,1,3) 
    imshow(ROIsMask,[0 1])
    title(strcat(...
        'ROIsMask (size = ',{' '},num2str(rSize),{' '},'pix. ;',...
        {' '},'thresh =',{' '},num2str(rThresh),{' '},';',...
        {' '},'min size =',{' '},num2str(rMinSize),{' '},'pix. ;',...
        {' '},'max size =',{' '},num2str(rMaxSize),{' '},'pix.)'));
    
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

%% Tophat

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

    % Tophat
    subplot(2,1,2) 
    imshow(Tophat,[qLow qHigh])
    title(strcat(...
        'Tophat (size =',{' '},num2str(TophatSize),{' '},'pix.)'));
    
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

%% AngleMap

choice = 1;
while choice > 0   

% Process .................................................................

    % Steerable Filter
    RawRSize = double(imresize(Raw,[nGridY nGridX],'nearest'));
    MaskRSize = double(imresize(Mask,[nGridY nGridX],'nearest'));
    [res,~,~,rot] = steerableDetector(RawRSize,sfOrder,sfSigma,180);
    for i=1:size(rot,3)
        temp = rot(:,:,i);
        temp(MaskRSize==0) = NaN;
        rot(:,:,i) = temp;
    end

    % Make AngleMap
    AngleMap = zeros(nGridY,nGridX);
    parfor i=1:nGridY
        for j=1:nGridX
            if ROIsMask(i,j) == 1            
                Crop = rot(i-(Pad1-1):i+(Pad1-1),j-(Pad1-1):j+(Pad1-1),:);
                tempMax = NaN(size(Crop,3),1);
                for k=1:size(Crop,3)
                    temp = Crop(:,:,k);
                    tempMax(k,1) = nanmean(temp(:));
                end
                [M,I] = max(tempMax);
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
        'Steerable filter (sigma =',{' '},num2str(sfSigma),{' '},'pix.)'));

    % AngleMap
    subplot(3,1,3)
    imshow(AngleMap,[0 180])
    title(strcat(...
        'AngleMap (Pad1 =',{' '},num2str(Pad1),{' '},'nROIs.)'));
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
            'Sigma :'
            };

        definput = {
            num2str(sfSigma)
            };
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        sfSigma = answer(1,1); 
        
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
        
        if ROIsMask(i,j) == 1  

            % Crop Data
            Angle = AngleMap(i,j);
            Crop = Tophat(rSize*i-((rSize*Pad2)-1):rSize*i+(rSize*(Pad2-1)),rSize*j-((rSize*Pad2)-1):rSize*j+(rSize*(Pad2-1)));

            % Determine 2D Correlation
            Corr2D = normxcorr2(Crop,Crop);
            Corr2D_Rot = double(imrotate(Corr2D,(90-Angle)*-1)); % rotate image acc. to angle 
            temp = max(Corr2D_Rot,[],2); [~,tempMid] = max(temp); % find max corr #1
            temp_AvgX = (nanmean(Corr2D_Rot(tempMid-2:tempMid+2,:),1))'; % average value on x axis
            temp = max(temp_AvgX,[],2); [~,tempMid] = max(temp); % find max corr #2 
            temp_AvgX = temp_AvgX/max(temp_AvgX); % max normalization
            temp_AvgX1 = temp_AvgX(tempMid:tempMid+rPad2-1);
            temp_AvgX2 = flip(temp_AvgX(tempMid-rPad2+1:tempMid));                
            Corr2D_AvgX = horzcat(temp_AvgX1,temp_AvgX2);
            Corr2D_AvgX = mean(Corr2D_AvgX,2);

            % Find Peaks
            [pks,locs,width,prom] = findpeaks(Corr2D_AvgX,'MinPeakDistance',6);
            if ~isempty(pks)
                [~,idx] = max(prom); % Get idx for the most prominent peak
                PeaksInfo = horzcat(i,j,pks(idx,1),locs(idx,1)*pixSize,width(idx,1),prom(idx,1));
            else
                PeaksInfo = NaN(1,6); % NaN if no peak found
            end

            % Append MergedData
            MergedData{i,j}{1,1} = Angle;
            MergedData{i,j}{2,1} = Crop;
            MergedData{i,j}{3,1} = Corr2D_AvgX; 
            MergedData{i,j}{4,1} = PeaksInfo; 
            
        end
    end
end

%% Get results

choice = 1;
while choice > 0   

% Process .................................................................

    % Merge "All" data
    All_PeaksInfo = NaN(0,0);
    All_Corr2D_AvgX = NaN(0,0);
    PromMap = zeros(nGridY,nGridX);
    for i=1:nGridY
        for j=1:nGridX
            if ROIsMask(i,j) == 1
                All_PeaksInfo = vertcat(All_PeaksInfo,MergedData{i,j}{4,1});
                All_Corr2D_AvgX(:,end+1) = MergedData{i,j}{3,1};
                PromMap(i,j) = MergedData{i,j}{4,1}(1,6);
            end
        end
    end

    % Get "Valid" data          
    Valid_PeaksInfo = NaN(0,0);
    Valid_Corr2D_AvgX = NaN(0,0);
    ValidMap = zeros(nGridY,nGridX);
    for i=1:size(All_PeaksInfo,1)
        tempY = All_PeaksInfo(i,1);
        tempX = All_PeaksInfo(i,2);
        tempProm = All_PeaksInfo(i,6);
        tempLoc = All_PeaksInfo(i,4);
        if tempProm >= minProm
            if tempLoc >= minLoc && tempLoc <= maxLoc
                Valid_PeaksInfo(end+1,:) = All_PeaksInfo(i,:);
                Valid_Corr2D_AvgX(:,end+1) = All_Corr2D_AvgX(:,i);
                ValidMap(tempY,tempX) = 1;
            end
        end
    end

    % Get Avg Corr2D_AvgX
    Avg_Corr2D_AvgX(:,1) = (0:rPad2-1)*pixSize; % distance (µm)
    Avg_Corr2D_AvgX(:,2) = mean(All_Corr2D_AvgX,2); % All 
    Avg_Corr2D_AvgX(:,3) = std(All_Corr2D_AvgX,0,2); % All S.D.
    if ~isempty(Valid_Corr2D_AvgX)
        Avg_Corr2D_AvgX(:,4) = mean(Valid_Corr2D_AvgX,2); % Valid
        Avg_Corr2D_AvgX(:,5) = std(Valid_Corr2D_AvgX,0,2); % Valid S.D.
    end
    
% Display .................................................................

    % Raw
    subplot(3,3,1:3)
    tempRaw = Raw;
    tempValidMap = imresize(ValidMap,[nYCrop nXCrop],'nearest');
    tempValidMap = bwmorph(logical(tempValidMap),'remove');
    tempValidMap = bwmorph(tempValidMap,'dilate');
    tempRaw(tempValidMap==1) = max(Raw(:));
    imshow(tempRaw,[qLow qHigh])
	title(strcat(...
        'ValidMap (minProm =',{' '},num2str(minProm),{' '},';',...
        {' '},'minLoc =',{' '},num2str(minLoc),{' '},';',...
        {' '},'maxLoc =',{' '},num2str(maxLoc),{' '},')'));

    % PromMap
    subplot(3,3,4:6)
    imshow(imresize(PromMap,[nYCrop nXCrop],'nearest'),[min(PromMap(:)) max(PromMap(:))])
    title('PromMap')
    colormap(gca, jet(256));
    colorbar(gca);

    % Avg_Corr2D_AvgX (All & Valid)
    subplot(3,3,7)
    plot(Avg_Corr2D_AvgX(:,1), Avg_Corr2D_AvgX(:,2))
    hold on
    plot(Avg_Corr2D_AvgX(:,1), Avg_Corr2D_AvgX(:,4))
    xlabel('Distance (µm)') 
    ylabel('Correlation')
    legend('All','Valid')
    xlim([0 ceil(max(Avg_Corr2D_AvgX(:,1)))])
    pbaspect([1 1 1])
    title('Avg. Corr2D')

    % Locs distribution ('All') 
    subplot(3,3,8)
    histogram(All_PeaksInfo(:,4),50)
    xlabel('Distance (µm)') 
    ylabel('Number of peaks')
    xlim([0 ceil(max(Avg_Corr2D_AvgX(:,1)))])
    pbaspect([1 1 1])
    title('Main peak localization (All)')

    % Locs distribution ('Valid')  
    subplot(3,3,9)
    histogram(Valid_PeaksInfo(:,4),50)
    xlabel('Distance (µm)') 
    ylabel('Number of peaks')
    xlim([0 ceil(max(Avg_Corr2D_AvgX(:,1)))])
    pbaspect([1 1 1])
    title('Main peak localization (Valid)')

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
            'minProm :',...
            'minLoc :',...
            'maxLoc :'
            };

        definput = {
            num2str(minProm),...
            num2str(minLoc),....
            num2str(maxLoc)
            };
        
        dlgtitle = 'Input'; dims = 1;
        answer = str2double(inputdlg(prompt,dlgtitle,dims,definput));
        
        minProm = answer(1,1); 
        minLoc = answer(2,1); 
        maxLoc = answer(3,1); 
        
        close
    end

    if choice == 0
        
        close
    end
    
end

%% Format and save

% Clear vars ..............................................................

% Format outputs ..........................................................

All_PeaksInfo = array2table(All_PeaksInfo,'VariableNames',...
    {'GridY','GridX','pks','loc','width','prom','valid'});
Valid_PeaksInfo = array2table(Valid_PeaksInfo,'VariableNames',...
    {'GridY','GridX','pks','loc','width','prom','valid'});
Avg_Corr2D_AvgX = array2table(Avg_Corr2D_AvgX,'VariableNames',...
    {'µm','All','Valid','AllSD','ValidSD'});
