
%% fullbkg
function BackgroundAndFlatfieldCorrectionOfTimeLapseImages(~,datename,channelstoinput,BACKGROUND)
% datename = '2015_12_15 smad3g smFISH';
C = datename;
E={datename}; %beginning of file name for saving
channelinputs =channelregexpmaker(channelstoinput);


 mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
    parentdir = mdir(1:b(end-1));
cd(strcat(parentdir,datename))
folderz = dir('*exp*');
folderzname = {folderz.name};
dirlog = [folderz.isdir];
explist = folderzname(dirlog);

for expdircell = explist
    expdirname = char(expdircell);
    DIRone = strcat(parentdir,datename,'\',expdirname);
    cd(DIRone);
    dirlist = dir('*.tif');
    [~,~,~,ScenesListed] = regexp([dirlist.name],'s[0-9]+');
    sceneList = unique(ScenesListed);
    [~,~,~,timeFrameList] = regexp([dirlist.name],'t[0-9]+');
    timeList = unique(timeFrameList);
    [~,~,~,channelsListed] = regexp([dirlist.name],channelinputs);
    channelList = unique(channelsListed);

    cd .. 
    mkdir(strcat(parentdir,datename,'\flatfield_corrected'));
    cd (strcat(parentdir,datename,'\flatfield_corrected'));
    
    
    
    
    %create the median stack of media only images for all timepoints for
    %each channel -- these will be used for flatfield correction in next
    %step
    for DIRflats = {DIRone} %choose directory (this is a remant of splitz)
        DIRflat = char(DIRflats);
        cd (DIRflat)
        
        for channel = channelList %cycle through one channel at a time
            chan = char(channel);
            
            normBkgimgMat = zeros(512,512,length(timeList));
            parfor i = 1:length(timeList)
                tpoints = timeList{i};
                bkimg = medianBKG(BACKGROUND,tpoints,chan); %load background images and compile into one median image
                bkShape = reshape(bkimg,[1 size(bkimg,2).*size(bkimg,1)]);
                bkSort = sort(bkShape(~isnan(bkShape)));
                normalizedBkgimg = (double(bkimg)./mean([bkSort(round(length(bkSort).*0.9999)) bkSort(round(length(bkSort).*0.99999))]));
                normBkgimgMat(:,:,i) = normalizedBkgimg;
            end
            [a,~] = regexp(chan,'_');
            chanstruct = chan;
            chanstruct(a) = [];
            bkgimgstruct.(chanstruct)=normBkgimgMat;
        end
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %make folders for each scene
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for scenedir = sceneList
        scene = scenedir{1};
    %     sceneN = str2num(scene(2:end));
    %     scenenum = scenestring(sceneN); 
        scenefile{1} = strcat(E{1},'_scene_',scene);
        mkdir(scenefile{1});
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for DIRflats = {DIRone} %choose directory (this is a remant of splitz)
        DIRflat = char(DIRflats);

        cd (DIRflat)

        for channel = channelList %cycle through one channel at a time
        filelist = dir(strcat('*',char(channel),'*.tif'));
        cfile = {filelist.name};


    %     for j = 1:length(cfile) %load image and determine scene and timepoint and channel
            parfor j = 1:length(cfile) %load image and determine scene and timepoint and channel
                filename = char(cfile(j));
                img = imread(filename);

                [a,b] = regexp(filename,'s[0-9]+'); %determine scene
                scene = filename(a:b);

                [c,d] = regexp(filename,'t[0-9]+'); %determine timepoint
                if ~isempty(c)
                    tpoint = filename(c:d);
                    tpoints = tpoint;
                else
                    tpoint='t00';
                    tpoints = '*';
                end


                [e,f] = regexp(filename,channelinputs); %determine channel
                chan = filename(e:f);
                    [a,~] = regexp(chan,'_');
                    chanstruct = chan;
                    chanstruct(a) = [];
                    tlog = strcmp(timeList,tpoints);
                    i = find(tlog==1);
        %         bkimg = medianBKG(BACKGROUND,tpoints,chan); %load background images and compile into one median image
        %         bkShape = reshape(bkimg,[1 size(bkimg,2).*size(bkimg,1)]);
        %         bkSort = sort(bkShape(~isnan(bkShape)));
        %         normalizedBkgimg = (double(bkimg)./mean([bkSort(round(length(bkSort).*0.9999)) bkSort(round(length(bkSort).*0.99999))]));
                normBkgimgMat = bkgimgstruct.(chanstruct);
                normalizedBkgimg = normBkgimgMat(:,:,i);

                flat = uint16(double(img)./normalizedBkgimg); %flatten the experimental image with the background image
                savename = strcat(E,'_flat_',scene,'_',tpoint,'_',chan,'.tif');
                SAVdir = strcat(parentdir,datename,'\flatfield_corrected\',E{1},'_scene_',scene);
                savethatimage(savename,SAVdir,flat,DIRflat,j);
%                 stophere=1;
            end
        end
    end
end

    % lineageyo(A,B);
    % cd('D:\Users\zeiss\Documents\MATLAB')
    % % bleachcorrection(A,B)
    % bleachcorrectionNew(A,B)
    % cd('D:\Users\zeiss\Documents\MATLAB')
    % toughsegmentationforstacks(A,B)
end


function bkimg = medianBKG(BACKGROUND,tpoints,chan)
    bkimgs = zeros(512,512,length(BACKGROUND));
    i = 1;
    for scenenumber = BACKGROUND
    scene = scenestring(scenenumber);
    % files = strcat('*',scene,'*',tpoints,'*',chan,'*');
    files = strcat('*',chan,'*',scene,'*',tpoints,'*');
    filelist = dir(files); cfile = {filelist.name}; filename = char(cfile);
    bkk = imread(filename);

    % %  if ~strcmp(chan,'c3')
    % if ~strcmp(chan,'DIC')
    %      bkklog = segmentationBKGsecond(bkk,[],scenenumber,[],[]);
    %      bkkg = regionfill(bkk,bkklog);
    %  else
    %      bkkg=bkk;
    % end
    %  
    bkkg=bkk;

     bkimgs(:,:,i) = bkkg;
    i=i+1;
    end
    bkimg = median(bkimgs,3);
end


function scene = scenestring(scenenumber)
    scene = char('s00');
        reg = num2str(scenenumber);
        c = size(reg,2);
        for b = 1:c
            scene(3+1-b) = reg(c+1-b);
        end
end

function tpoint = timestring(time)
    tpoint = char('t00');
        reg = num2str(time);
        c = size(reg,2);
        for b = 1:c
            tpoint(3+1-b) = reg(c+1-b);
        end
end
 
function savethatimage(savename,SAVdir,flat,DIRone,j)
    disp(strcat(savename,'...',num2str(j)));
    cd (SAVdir);
    imwrite(flat,char(savename),'tiff');
    cd (DIRone);
end



function channelinputs =channelregexpmaker(channelstoinput)
    channelinputs = '(';
    for i=1:length(channelstoinput) % creates a string of from '(c1|c2|c3|c4)' for regexp functions
        if i ==1
        channelinputs = strcat(channelinputs,channelstoinput{i});
        elseif i < length(channelstoinput)
            channelinputs = strcat(channelinputs,'|',channelstoinput{i});
        else
            channelinputs = strcat(channelinputs,'|',channelstoinput{i},')');
        end
    end
end


