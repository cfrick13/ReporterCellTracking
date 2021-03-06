
%% fullbkg
function preprocessBKGandFLATnew(A,B)
% A = 'D:\Frick\';
% B = '2015_12_15 smad3g smFISH';
C = B;

E={B}; %beginning of file name for saving
F = {'_back_c1_t','_back_c2_t','_back_c3_t'};
FF = {'*c1_*.tif','*c2_*.tif','*c3_*.tif'};
% BACKGROUND = [13 14 15 16 17];
% BACKGROUND = [28 29 30 31];
% BACKGROUND = [7:11];
BACKGROUND = [25:29];
BACKGROUND = [1:12];


BACKdir = strcat(A,B,'\back');mkdir (BACKdir);
DIRone = strcat(A,B,'\',C);

cd(DIRone);
dirlist = dir('*.tif');
[~,~,~,ScenesListed] = regexp([dirlist.name],'s[0-9]+');
sceneList = unique(ScenesListed);
[~,~,~,timeFrameList] = regexp([dirlist.name],'t[0-9]+');
timeList = unique(timeFrameList);
[~,~,~,channelsListed] = regexp([dirlist.name],'c[0-9]+');
channelList = unique(channelsListed);

cd .. 
mkdir(strcat(A,B,'\flatfield_corrected'));
cd (strcat(A,B,'\flatfield_corrected'));
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
    DIRflatcomp = DIRflat(end-5:end);
    DIRonecomp = DIRone(end-5:end);
  
    cd (DIRflat)

for channel = channelList %cycle through one channel at a time
filelist = dir(strcat('*',char(channel),'*.tif'));
cfile = {filelist.name};

% for j = 1:length(cfile) %load image and determine scene and timepoint and channel

for j = 1:length(cfile) %load image and determine scene and timepoint and channel
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
    
    [e,f] = regexp(filename,'c[0-9]+'); %determine channel
    chan = filename(e:f);
    
    bkimg = medianBKG(BACKGROUND,tpoints,chan); %load background images and compile into one median image
    bkShape = reshape(bkimg,[1 size(bkimg,2).*size(bkimg,1)]);
    bkSort = sort(bkShape(~isnan(bkShape)));
    normalizedBkgimg = (double(bkimg)./mean([bkSort(round(length(bkSort).*0.9999)) bkSort(round(length(bkSort).*0.99999))]));
    
    flat = uint16(double(img)./normalizedBkgimg); %flatten the experimental image with the background image
    savename = strcat(E,'_flat_',scene,'_',tpoint,'_',chan,'.tif');
    SAVdir = strcat(A,B,'\flatfield_corrected\',E{1},'_scene_',scene);
    savethatimage(savename,SAVdir,flat,DIRflat,j);
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
files = strcat('*',scene,'*',tpoints,'*',chan,'*');
filelist = dir(files); cfile = {filelist.name}; filename = char(cfile);
bkimgs(:,:,i) = imread(filename);
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

function lineageyo(A,B)
%%
F = {'_c1_ORG.tif','_c2_ORG.tif','_c3_ORG.tif'};
FF = {'*c1_*.tif','*c2_*.tif','*c3_*.tif'};
cd (strcat(A,B));

cd ('flatfield_corrected')

%%
folderlist = dir('*_s*');
folderlist.name;

for folder = {folderlist.name}
cd(folder{1})
mkdir('c1_flat')
mkdir('c2_flat')
mkdir('c3_flat')
mkdir('tiffs');

filelist = dir('*c1*.tif');  
for cfile={filelist.name}
       filepath = which(cfile{1});
       disp(filepath);
       movefile(filepath,'c1_flat')
end

filelist = dir('*c2*.tif');
for cfile={filelist.name}
       filepath = which(cfile{1});
       disp(filepath);
       movefile(filepath,'c2_flat')
end

filelist = dir('*c3*.tif');
for cfile={filelist.name}
       filepath = which(cfile{1});
       disp(filepath);
       movefile(filepath,'c3_flat')
end
 cd .. 

end
end



