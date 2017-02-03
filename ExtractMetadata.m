% reader = bfGetReader('/Volumes/Seagate Backup Plus Drive/FrickData/2017_01_27 plate/exp1 timelapse of all snail clones in NGSmad3 cells.czi');
reader = bfGetReader('/Volumes/Seagate Backup Plus Drive/FrickData/2017_01_25 plate/exp1 timelapse of caga and pctgf in cmvNGSmad3ex1 cells.czi');
omeMeta = reader.getMetadataStore();

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
% 
% disp(thec)
% disp(thet)
% disp(timeVec./60)


datastruct  = struct();

datastruct.imageDimensions = dimensions;
datastruct.timeCount= timeCount;
datastruct.channelCount = channelCount;
datastruct.sceneCount = imageCount; 
% datastruct.channelNames = 

% t(theT+1) = y.value