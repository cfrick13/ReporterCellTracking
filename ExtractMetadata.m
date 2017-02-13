% reader = bfGetReader('/Volumes/Seagate Backup Plus Drive/FrickData/2017_01_27 plate/exp1 timelapse of all snail clones in NGSmad3 cells.czi');

%set directory to location of code being used (generally external harddrive
mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
    parentdir = mdir(1:b(end-1));
cd(parentdir)
readerpath = strcat(parentdir,'Tracking/bfmatlab/');
% denoisepath = strcat(parentdir,'Tracking/coherencefilter_version5b/');
addpath(readerpath);
exportdir = strcat(parentdir,'Tracking/Export/');

cd(parentdir)
[FileName,PathName,~] = uigetfile('*.*');
cd(PathName)
reader = bfGetReader(FileName);
omeMeta = reader.getMetadataStore();
% reader = bfGetReader('/Volumes/Seagate Backup Plus Drive/FrickData/2017_01_25 plate/exp1 timelapse of caga and pctgf in cmvNGSmad3ex1 cells.czi');
% omeMeta = reader.getMetadataStore();

% omeMeta = data{1,4};

%determine dimensions of images
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
dimensions = [stackSizeX stackSizeY stackSizeZ]; %[x y z];

timeCount = omeMeta.getPixelsSizeT(0).getValue();
planeCount = omeMeta.getPlaneCount(0); % a plane is number of images for a certain scene (each channel at each timepoint)
channelCount = omeMeta.getChannelCount(0); %number of channels
imageCount = omeMeta.getImageCount(); %number of scenes in experiment
% sceneCount = size(data,1);

t=[];
thet=[];
thec=[];
timeVec=[];

planeCountMat = 0:planeCount-1;
for i = 0:1:imageCount-1
    for pp = 1:timeCount
        p = planeCountMat(pp*3);
    
        deltaT_Hash = omeMeta.getPlaneDeltaT(i,p); %Time getPlaneDeltaT(int imageIndex,int planeIndex)
        deltaT = double(deltaT_Hash.value);
        theT_Hash = omeMeta.getPlaneTheT(i,p); %getPlaneTheT(int imageIndex, int planeIndex)
        theT = double(theT_Hash.getValue);
        theC_Hash = omeMeta.getPlaneTheC(i,p); %getPlaneTheT(int imageIndex, int planeIndex)
        theC = double(theC_Hash.getValue);
        
        deltaTarray(pp) = deltaT;
        thet(pp)=theT+1;
        thec(pp) = theC+1;
        if theT==0
            t(theT+1) = deltaT;
        else
%             t(theT+1) = t(theT)+deltaT;
            t(theT+1) = deltaT;
        end
        %extract the times

    end
    timeVec(i+1,:) = t./60; %save in units of minutes
    clear t
end




for i = 1:channelCount
    cN = omeMeta.getChannelName(0,i-1);
    cname = char(cN);
    a = (cname==' ');
    if sum(a)==1
        g = find(a==1);
        cname = cname((g+1):end);
    end
    channelNames{i} = cname;
end


datastruct  = struct();

datastruct.imageDimensions = dimensions;
datastruct.timeCount= timeCount;
datastruct.channelCount = channelCount;
datastruct.sceneCount = imageCount; 
datastruct.channelNames = channelNames;


[~,b] = regexp(FileName,'exp[0-9]');
savename = FileName(1:b);
cd(exportdir)
savename = strcat(savename,'-metaData.mat');
clear omeMeta reader deltaT_Hash theC_Hash theT_Hash FileName
save(savename)

makeDoseStructFromXLS(savename)