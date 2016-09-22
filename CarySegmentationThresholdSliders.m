function  CarySegmentationThresholdSliders
%Loosely based on the paper by Arce et al, 2013.
%http://www.nature.com/articles/srep02266  
%Fast and accurate automated cell boundary determination for fluorescence microscopy

global sldPM sldNuc sldCell threshinput PMSegColor NucSegColor CytoSegColor PMerode Imagez ChannelHandle selectedchannel FluorChannelOne FluorChannelTwo NucleusColor CellColor imageaxis f SegmentationFileDirectory quantaxis DateDirectory ImagePositions PositionHandle
close all
% FluorChannelOne = 'H2A';
% FluorChannelTwo = 'VP40';
FluorChannelOne = 'green';
FluorChannelTwo = 'cyan';
% DateOfExperiment = '2016_01_01';
NucleusColor = 3; %1=red, 2=green, 3=blue
CellColor = 2; %1=red, 2=green, 3=blue
threshinput = struct(); %DO NOT ALTER.... keep as struct() initially then adjust manually using the userinterface.
    threshinput.Nuc = []; threshinput.Cell = []; threshinput.PM = [];
    
NucSegColor = 1; %1=red, 2=green, 3=blue, [1 2] =yellow, [1 3]=purple, [2 3]=cyan 
CytoSegColor = [1 3]; %1=red, 2=green, 3=blue, [1 2] =yellow, [1 3]=purple, [2 3]=cyan 
PMSegColor = [1 2]; %1=red, 2=green, 3=blue, [1 2] =yellow, [1 3]=purple, [2 3]=cyan 

PMerode = 16; %define the number of pixels by which to erode the cytoplasm with for creating the pm segmentation 
imagescalefactor = 1;
Imagez = struct();

%Find the directory to analyze
MainDirectory = uigetdir([],'Choose the Experiment Folder Containing the Images you wish to analyze');
% MainDirectory = strcat('/Users/frick/Documents/Goentoro_Lab/Collaborations/Cary/');
DateDirectory = strcat(MainDirectory);
SegmentationFileDirectory = strcat(MainDirectory,'/Segmentation');
mkdir(SegmentationFileDirectory)



ChannelOptions = {'Nucleus','Cell','PM','Nuc+Cell overlay'};
selectedchannel = ChannelOptions{end};
ImagePositions = updatePositionList;

%set the size of the user interface figure
f = figure(1);
f.Units = 'normalized';
f.Position =[0.02,0.02,0.8,0.8];
% f.Position =[10,10,800,1100];

%define the positions of the buttons
xspotone = 0.7;
xspottwo = 0.83;
yspot = sort(0.10:0.03:0.95,'descend');
buttonwidth=0.12;
buttonheight=0.025;



i=1;
%define the user interface buttons on right hand side
uicontrol('Style','text','String','Choose channel to observe | or Push 1,2,3 or 4 to select channel to observe',...
          'Units','normalized','FontSize',8,...
          'Position',[xspottwo,yspot(i),buttonwidth,buttonheight.*2]); i=i+1;
ChannelHandle = uicontrol('Style','popupmenu',...
          'String',ChannelOptions',...
          'Units','normalized','FontSize',10,...
          'Position',[xspottwo,yspot(i),buttonwidth.*1.2,buttonheight],...
          'Callback',@channelpopup_menu_Callback); i=i+3;
      ChannelHandle.Value = length(ChannelHandle.String);

%Give text to explain hotkeys
uicontrol('Style','text','String',...
                '*Add area hotkey is "a" (works for selected channel when one channel is selected)',...
          'Units','normalized','FontSize',8,...
          'Position',[xspottwo,yspot(i),buttonwidth,buttonheight.*2]); i=i+1;
      
uicontrol('Style','text','String',...
                '*Remove area hotkey is "r" (works for selected channel when one channel is selected)',...
          'Units','normalized','FontSize',8,...
          'Position',[xspottwo,yspot(i),buttonwidth,buttonheight.*2]); i=i+1;
uicontrol('Style','text','String',...
                '*Delete area hotkey is "d" (works for selected channel when one channel is selected)',...
          'Units','normalized','FontSize',8,...
          'Position',[xspottwo,yspot(i),buttonwidth,buttonheight.*2]); i=i+2;
      
      
uicontrol('Style','text','String',...
                '*Erode area hotkey is "e", dilate is "l" (the width erosion/dilation is the number of clicks)',...
          'Units','normalized','FontSize',8,...
          'Position',[xspottwo,yspot(i),buttonwidth,buttonheight.*2]); i=i+5;
      
      
uicontrol('Style','text','String',...
                'Save segmentation hotkey is "s" (works at anytime)',...
          'Units','normalized','FontSize',8,...
          'Position',[xspottwo,yspot(i),buttonwidth,buttonheight.*2]); i=i+2;
      
uicontrol('Style','text','String',...
                'Quantify hotkey hotkey is "q" (works at anytime)',...
          'Units','normalized','FontSize',8,...
          'Position',[xspottwo,yspot(i),buttonwidth,buttonheight.*2]); i=i+3;
             
 
%define the user interface buttons on left hand side
i=1;
uicontrol('Style','text','String','Choose Position to analyze | or push <- or -> to scroll through experiment positions',...
          'Units','normalized','FontSize',8,...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight.*2]); i=i+1;
PositionHandle = uicontrol('Style','popupmenu',...
          'String',ImagePositions','FontSize',10,...
          'Units','normalized',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@positionpopup_menu_Callback); i=i+2;
 
    % ADD AREA
uicontrol('Style','pushbutton','String','AddArea)',...
          'Units','normalized',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@addareabuttonNUC_Callback); i=i+1;
      
      %REMOVE AREA
uicontrol('Style','pushbutton','String','RemoveArea',...
          'Units','normalized',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@removeareabuttonNUC_Callback); i=i+1;
      
      %DELETE AREA
uicontrol('Style','pushbutton','String','Delete (Nucleus)',...
          'Units','normalized',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@deletebuttonNUC_Callback); i=i+1;  
      
      %ERODE
uicontrol('Style','pushbutton','String','Erode (e)',...
          'Units','normalized',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@erode_Callback); i=i+1;  
      %DILATE
uicontrol('Style','pushbutton','String','Dilate (l)',...
          'Units','normalized',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@dilate_Callback); i=i+2;   
      
      %SEGMENT
uicontrol('Style','pushbutton','String','Re-segment',...
          'Units','normalized',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@segmentbutton_Callback); 
      %RE-SEGMENT only Plamsa membrane
uicontrol('Style','pushbutton','String','Re-segment PM only',...
          'Units','normalized',...
          'Position',[xspottwo,yspot(i),buttonwidth,buttonheight],...
          'Callback',@updatePMsegCallBack);i=i+3;

      %save
uicontrol('Style','pushbutton',...
          'Units','normalized',...
          'String','Save Segmentation',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@savebutton_Callback); i=i+2; 
      
      %Quantify
uicontrol('Style','pushbutton',...
          'Units','normalized',...
          'String','QUANTIFY',...
          'Position',[xspotone,yspot(i),buttonwidth,buttonheight],...
          'Callback',@quantify); i=i+3;
      
 %NUCLEUS     
uicontrol('Style','text','String',...
                'Manually Threshold for Nucleus segmentation',...
          'Units','normalized','FontSize',8,...
          'Position',[xspotone,yspot(i),buttonwidth.*2,buttonheight.*1.3]);
      
sldNuc = uicontrol('Style', 'slider',...
        'Units','normalized',...
        'String','Nuc',...
        'Min',1,'Max',255,'Value',200,...
        'Position', [xspotone yspot(i) buttonwidth.*2 buttonheight./2],...
        'Callback', @PMthreshslider); i=i+2;
 
%CELL    
uicontrol('Style','text','String',...
                'Manually Threshold for Cyto segmentation',...
          'Units','normalized','FontSize',8,...
          'Position',[xspotone,yspot(i),buttonwidth.*2,buttonheight.*1.3]);
      
sldCell = uicontrol('Style', 'slider',...
        'Units','normalized',...
        'String','Cell',...
        'Min',1,'Max',255,'Value',200,...
        'Position', [xspotone yspot(i) buttonwidth.*2 buttonheight./2],...
        'Callback', @PMthreshslider); i=i+2;
 
%PM
uicontrol('Style','text','String',...
                'Manually Threshold PM high intesnity projection segmentation',...
          'Units','normalized','FontSize',8,...
          'Position',[xspotone,yspot(i),buttonwidth.*2,buttonheight.*1.3]);
      
sldPM = uicontrol('Style', 'slider',...
        'Units','normalized',...
        'String','PM',...
        'Min',1,'Max',255,'Value',200,...
        'Position', [xspotone yspot(i) buttonwidth.*2 buttonheight./2],...
        'Callback', @PMthreshslider); 
    

      
f.Visible = 'on';
f.Units = 'normalized';

%set the size of the image axis
imageaxis = axes();
imageaxis.Position = [0.1 0.1 0.4 0.4];
imageaxis.Units  = 'pixels';
Pos = imageaxis.Position;
Pos(3:4) = 512.*imagescalefactor;
imageaxis.Position = Pos;
imageaxis.Units = 'normalized';


quantaxis = axes();
quantaxis.Position = [0.7 0.05 0.2 0.2];
% quantaxis.Units  = 'pixels';
% Pos = quantaxis.Position;
% Pos(3:4) = 512.*imagescalefactor;
% quantaxis.Position = Pos;
quantaxis.Units = 'normalized';


set(f,'KeyPressFcn',@keypress);
end


function PMthreshslider(source,callbackdata)
    global threshinput Imagez
    str = source.String;
     threshinput.(str) =source.Value;
     
    if strcmp(source.String,'PM')
    elseif strcmp(source.String,'Nuc')
        Imagez.NucSeg = segmentNucleus(Imagez.NucImage); %semgnet the image
    elseif strcmp(source.String,'Cell')
        Imagez.CellSeg = segmentCell(Imagez.CellImage);
    end
updatePMseg
updateImage
end
    
%allows for hotkeys
function keypress(fig_obj,~)
global ChannelHandle selectedchannel
key = get(fig_obj,'CurrentKey');




switch key
    case '1'
            ChannelHandle.Value = 1;
            s = ChannelHandle.String;
            selectedchannel = s{ChannelHandle.Value};
            updateImage
    case '2'
            ChannelHandle.Value = 2;
            s = ChannelHandle.String;
            selectedchannel = s{ChannelHandle.Value};
            updateImage
    case '3'
            ChannelHandle.Value = 3;
            s = ChannelHandle.String;
            selectedchannel = s{ChannelHandle.Value};
            updateImage
    case '4'
            ChannelHandle.Value = 4;
            s = ChannelHandle.String;
            selectedchannel = s{ChannelHandle.Value};
            updateImage
            

    case 'a'
        if strcmp(selectedchannel,'Nucleus')
            addareabuttonNUC_Callback([],[])
        elseif strcmp(selectedchannel,'Cell')
            addareabuttonNUC_Callback([],[]) 
        elseif strcmp(selectedchannel,'PM')
            addareabuttonNUC_Callback([],[]) 
        end
    case 'r' 
        if strcmp(selectedchannel,'Nucleus')
            removeareabuttonNUC_Callback([],[])
        elseif strcmp(selectedchannel,'Cell')
            removeareabuttonNUC_Callback([],[])
        elseif strcmp(selectedchannel,'PM')
            removeareabuttonNUC_Callback([],[]) 
        end
    case 'd'     
        if strcmp(selectedchannel,'Nucleus')
            deletebuttonNUC_Callback([],[])
        elseif strcmp(selectedchannel,'Cell')
            deletebuttonNUC_Callback([],[])
        elseif strcmp(selectedchannel,'PM')
            deletebuttonNUC_Callback([],[])
        end   
      
    case 'e'
        if strcmp(selectedchannel,'Nucleus')
            erode_Callback
        elseif strcmp(selectedchannel,'Cell')
            erode_Callback
        elseif strcmp(selectedchannel,'PM')
            erode_Callback
        end    
        
   case 'l'
        if strcmp(selectedchannel,'Nucleus')
            dilate_Callback
        elseif strcmp(selectedchannel,'Cell')
            dilate_Callback
        elseif strcmp(selectedchannel,'PM')
            dilate_Callback
        end    
        
    case 'rightarrow'
        nextpos([],[])
    case 'leftarrow'
        prevpos([],[])
    case 's'
        savebutton_Callback([],[])
    case 'q'
        quantify
end

end

%identify the positions present
function ImagePositions = updatePositionList
global DateDirectory
cd (DateDirectory); %change to the folder containing images
Imagelist = dir('*.tif');  %look for all .tif files in folder
ImageFilenames = {Imagelist.name}; %convert file structure into array filenames
% [~,~,~,d] = regexp(ImageFilenames,'_[0-9]+.tif');
[~,~,~,d] = regexp(ImageFilenames,'c[0-9]+'); %search all strings in filename array for position numbers
arrayofpositions = cellfun(@(x) x{1},d,'UniformOutput',0); %convert cell array of cells into cell array of strings
ImagePositions = unique(arrayofpositions); %eliminate duplicates to determine how many unique, individual positions are present
end
%chooses the position
function positionpopup_menu_Callback(source,~) 
global Position 
% Determine the selected data set.
 str = source.String;
 val = source.Value;
 Position = char(str{val});


setSceneAndUpdate
end

function channelpopup_menu_Callback(source,~) 
global selectedchannel 
% Determine the selected data set.
 str = source.String;
 val = source.Value;
 selectedchannel = char(str{val});


setSceneAndUpdate
end


%chooses the files and load and display them
function setSceneAndUpdate
global FluorChannelOne FluorChannelTwo Position Imagez SegmentationFileDirectory DateDirectory PMerode

%load the fluorescence image files
cd(DateDirectory)
   NucleusImageFilelist = dir(strcat('*',FluorChannelOne,'*',Position,'*.tif'));  %the purple text in quotes is the fluorescent channel used
       NucleusImageFilename = char(NucleusImageFilelist.name); %convert the directory search structure to the filename
       Imagez.NucImage = imread(NucleusImageFilename); %load image of nuclues (e.g. H2A)

   Imagez.CellImageFilelist = dir(strcat('*',FluorChannelTwo,'*',Position,'*.tif')); %the purple text in quotes is the fluorescent channel used
       Imagez.CellImageFilename = char(Imagez.CellImageFilelist.name); %convert the directory search structure to the filename
       Imagez.CellImage = imread(Imagez.CellImageFilename); %load image of cell (e.g. VP40) 

       
%load the segmentation files or run segmentation if no saved files are present
cd(SegmentationFileDirectory)
   NucleusSegFilelist = dir(strcat('Segment*',FluorChannelOne,'*',Position,'*.mat')); 
   if ~isempty(NucleusSegFilelist)
       NucleusSegFilename = char(NucleusSegFilelist.name); %convert the directory search structure to the filename
       ld = load(NucleusSegFilename); %load image of nuclues (e.g. H2A)
       Imagez.NucSeg = ld.NucSeg;
   else
       Imagez.NucSeg =[];
   end
       
   Imagez.CellSegFilelist = dir(strcat('Segment*',FluorChannelTwo,'*',Position,'*.mat')); 
   if ~isempty(Imagez.CellSegFilelist) %only attempt to load segmentation files if they are present (empty if not)
       Imagez.CellSegFilename = char(Imagez.CellSegFilelist.name); %convert the directory search structure to the filename
       ld = load(Imagez.CellSegFilename); %load image of nuclues (e.g. H2A)
       Imagez.CellSeg = ld.CellSeg;
   else
       Imagez.CellSeg =[];
   end
   
   %if there aren't previously saved segmentation files, then run segmentation
   if isempty(Imagez.CellSeg) || isempty(Imagez.NucSeg) 
       segmentbutton_Callback([],[])
   end
   
   
      Imagez.PMSegFilelist = dir(strcat('Segment*','PM','*',Position,'*.mat')); 
   if ~isempty(Imagez.PMSegFilelist) %only attempt to load segmentation files if they are present (empty if not)
       Imagez.PMSegFilename = char(Imagez.PMSegFilelist.name); %convert the directory search structure to the filename
       ld = load(Imagez.PMSegFilename); %load image of nuclues (e.g. H2A)
       Imagez.PMSeg = ld.PMSeg;
   else
       Imagez.PMSeg =[];
   end
   
   if isempty(Imagez.PMSeg);
   Cyto = Imagez.CellSeg;
   se = strel('disk',PMerode);
   CytoE = imerode(Cyto,se);
   Imagez.CellSeg = CytoE;
   updatePMseg([],[]) %execute the function
   end
   
  
cd(DateDirectory)
updateImage
end

function updatePMseg(~,~)
global  Imagez PMerode
PM = false(size(Imagez.CellSeg));

se = strel('disk',PMerode);
CytoE = logical(Imagez.CellSeg); % select the whole cell segmentation;
CytoED = imdilate(CytoE,se);
PM(CytoED) = 1; 
PM(CytoE) = 0;

highintensityprotrusions = segmentPM(Imagez.CellImage); %
PM(logical(highintensityprotrusions)) = 1; 
% CytoEperimD = imdilate(CytoEperim,se); % widen the perimet by dilating
Imagez.PMSeg = PM; %save it to the Imagez structure

CytoE(PM)=0;
% CytoE(logical(Imagez.NucSeg)) = 0;
Imagez.CellSeg = CytoE; % update the cytoplasm segmentation to exclude the PM segmentation
end

function updatePMsegCallBack(~,~)
updatePMseg([],[])
updateImage
end
    
    
function updateImage
global PMSegColor NucSegColor CytoSegColor Imagez NucleusColor CellColor f imageaxis selectedchannel

%determine which channel to display
if strcmp(selectedchannel,'Nuc+Cell overlay')
    C = uint8(zeros([size(Imagez.NucImage,1) size(Imagez.NucImage,2) 3]));
    C(:,:,NucleusColor) = Imagez.NucImage;
    C(:,:,CellColor) = Imagez.CellImage; 
    
    % set the nuclear segmentation to be the color
    for c = NucSegColor
    colorIM = C(:,:,c);
        imPerimeter = imdilate(bwperim(Imagez.NucSeg),strel('disk',1));
        colorIM(logical(imPerimeter)) = max(max(max(C)));
        C(:,:,c) = colorIM;
    end
    
    % set the cytoplasmic segmentation to be the color 
    for c = CytoSegColor
    colorIM = C(:,:,c);
        imPerimeter = imdilate(bwperim(Imagez.CellSeg),strel('disk',1));
        colorIM(logical(imPerimeter)) = max(max(max(C)));
        C(:,:,c) = colorIM;
    end
    
    % set the plasma membrane (PM) segmentation to be the color
    for c = PMSegColor
    colorIM = C(:,:,c);
        imPerimeter = bwperim(Imagez.PMSeg);
%         imPerimeter = imdilate(bwperim(Imagez.PMSeg),strel('disk',2));
        colorIM(logical(imPerimeter)) = max(max(max(C)));
        C(:,:,c) = colorIM;
    end
        
elseif strcmp(selectedchannel,'Nucleus')
    C = uint8(zeros([size(Imagez.NucImage,1) size(Imagez.NucImage,2) 3]));
    C(:,:,NucleusColor) = Imagez.NucImage;
   
   % set the nuclear segmentation to be the color red (1 is red for RGB image)
   for c = NucSegColor
    colorIM = C(:,:,c);
        imPerimeter = imdilate(bwperim(Imagez.NucSeg),strel('disk',1));
        colorIM(logical(imPerimeter)) = max(max(max(C)));
        C(:,:,c) = colorIM;
   end
    
elseif strcmp(selectedchannel,'Cell')
    C = uint8(zeros([size(Imagez.NucImage,1) size(Imagez.NucImage,2) 3]));
    C(:,:,CellColor) = Imagez.CellImage; 
    
    % set the cytoplasmic segmentation to color 
   for c = CytoSegColor
    colorIM = C(:,:,c);
        imPerimeter = imdilate(bwperim(Imagez.CellSeg),strel('disk',1));
        colorIM(logical(imPerimeter)) = max(max(max(C)));
        C(:,:,c) = colorIM;
   end
   
elseif strcmp(selectedchannel,'PM')
    C = uint8(zeros([size(Imagez.NucImage,1) size(Imagez.NucImage,2) 3]));
    C(:,:,CellColor) = Imagez.CellImage; 
    
    % set the cytoplasmic segmentation to color 
   for c = PMSegColor
    colorIM = C(:,:,c);
        imPerimeter = bwperim(Imagez.PMSeg);
        colorIM(logical(imPerimeter)) = max(max(max(C)));
        C(:,:,c) = colorIM;
   end
   
end


bkg = true(size(Imagez.CellSeg));
celldilate = (logical(Imagez.CellSeg) | logical(Imagez.PMSeg));
celldilate = imdilate(celldilate,strel('disk',50));
bkg(celldilate) = 0;
Imagez.bkg = bkg;

%make a figure with the image
figure(f)
axes(imageaxis)
imagesc(C);
titlestr = Imagez.CellImageFilename;
[a] = regexp(titlestr,'_');
titlestr(a) = '-';
title(titlestr,'FontSize',12);
end

function quantify(~,~)
global quantaxis f Imagez FluorChannelTwo
figure(f)
axes(quantaxis)

cl = Imagez.CellImage;
background = mean(cl(logical(Imagez.bkg)));
CytoSeg = false(size(Imagez.NucSeg));
CytoSeg(logical(Imagez.CellSeg))=1;
CytoSeg(logical(Imagez.NucSeg))=0;

CellQuant = Imagez.CellImage(logical(Imagez.CellSeg));
NucQuant = Imagez.CellImage(logical(Imagez.NucSeg));
CytoQuant = Imagez.CellImage(logical(CytoSeg));
CytoNucRatio = median(CytoQuant)./median(NucQuant);
PMQuant = Imagez.CellImage(logical(Imagez.PMSeg));

ValueVector = [median(CellQuant) median(NucQuant) median(CytoQuant) median(PMQuant)]-background;
RatioVector = ValueVector(3)./ValueVector(2);
b = bar(horzcat(ValueVector,RatioVector),'FaceColor',[0 1 1]);
quantaxis.XTickLabels = {'cell','nucleus','cytoplasm','PM','Cyto/NucRatio'};
quantaxis.XTickLabelRotation = 45;
title(strcat('quantified ',FluorChannelTwo,' fluorescence of cellular compartments'));
end

%segmentation code
function If = segmentNucleus(img)
global NucimgSmoothingFactor threshinput sldNuc
imgsize = 512;
ratio = imgsize./2048;

% if isa(img,'uint8')
%     maxint=255;
% else
%     maxint = 65535;
% end


NucimgSmoothingFactor = 30;
sigmafirst = NucimgSmoothingFactor.*3.*ratio;
kernelgsizefirst = NucimgSmoothingFactor.*6.*ratio;

% start
imgorig = img;

img = wiener2(img,[5 5]);
se =strel('disk',2);
    Ie = imerode(imgorig,se);
    Iobr = imreconstruct(Ie,img);
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    IobrcbrF = imcomplement(Iobrcbr);
    refinedImage = double(IobrcbrF);

se =strel('disk',NucimgSmoothingFactor);
    Ie = imerode(refinedImage,se);
    Iobr = imreconstruct(Ie,refinedImage);
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    refinedImage = double(Iobrcbr);


sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
    gauss = gaussianBlurz(double(refinedImage),sigma,kernelgsize);

%See paper
%http://www.nature.com/articles/srep02266  
%Fast and accurate automated cell boundary determination for fluorescence microscopy
sub = double(refinedImage) -double(gauss);%%%%%%% key step!
b = find(sub == min(min(sub)));
rattio = gauss(b)./refinedImage(b);
gaustwocorr = gauss./rattio;
sub_scale_corr = double(refinedImage) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
% [numbers,bincenters] = hist(double(vec),0:0.5:maxint);
[numbers,bincenters] = hist(double(vec),min(vec):0.5:max(vec));


numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);

%%%%%%%%%%%%%%%%%%%%
left=0.3*mf;
slopedown=0.2*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
insideslopeup = find(fraction(leftedge+insideslopedown:end) > left,1,'first');
% threshlocation = bincenters(leftedge+insideslopedown);
threshlocation = bincenters(leftedge+insideslopedown+insideslopeup);

if ~isempty(threshinput.Nuc)
    threshlocation = threshinput.Nuc;
end
% %%%
% figure(191)
% bar(bincenters,fraction),xlim([min(vec) max(vec)]);hold on
% stem(threshlocation,mf)
% %%%%


subtractionthreshold = threshlocation;
    if size(subtractionthreshold,1)==size(subtractionthreshold,2)
    else
         subtractionthreshold = mean(threshlocation);
    end


subtracted = sub_scale_corr-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
Ie(Ie>0)=50;
If = imfill(Ie,'holes');
sldNuc.Value = threshlocation;
end
function If = segmentCell(img)
global CellimgSmoothingFactor threshinput sldCell

imgsize = size(img);
ratio = imgsize(1)./2048;

% if isa(img,'uint8')
%     maxint=255;
% else
%     maxint = 65535;
% end

CellimgSmoothingFactor = 30;
sigmafirst = CellimgSmoothingFactor.*3.*ratio;
kernelgsizefirst = CellimgSmoothingFactor.*6.*ratio;

% start
imgorig = img;

img = wiener2(img,[5 5]);
se =strel('disk',2);
    Ie = imerode(imgorig,se);
    Iobr = imreconstruct(Ie,img);
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    IobrcbrF = imcomplement(Iobrcbr);
    refinedImage = double(IobrcbrF);

se =strel('disk',CellimgSmoothingFactor);
    Ie = imerode(refinedImage,se);
    Iobr = imreconstruct(Ie,refinedImage);
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    refinedImage = double(Iobrcbr);


sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
    gauss = gaussianBlurz(double(refinedImage),sigma,kernelgsize);

%See paper
%http://www.nature.com/articles/srep02266  
%Fast and accurate automated cell boundary determination for fluorescence microscopy
sub = double(refinedImage) -double(gauss);%%%%%%% key step!
b = find(sub == min(min(sub)));
rattio = gauss(b)./refinedImage(b);
gaustwocorr = gauss./rattio;
sub_scale_corr = double(refinedImage) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
% [numbers,bincenters] = hist(double(vec),0:0.5:maxint);
[numbers,bincenters] = hist(double(vec),min(vec):0.5:max(vec));


numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);

%%%%%%%%%%%%%%%%%%%%
left=0.3*mf;
slopedown=0.2*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
insideslopeup = find(fraction(leftedge+insideslopedown:end) > left,1,'first');
% threshlocation = bincenters(leftedge+insideslopedown);
threshlocation = bincenters(leftedge+insideslopedown+insideslopeup);


if ~isempty(threshinput.Cell)
    threshlocation = threshinput.Cell;
end

% threshlocation = bincenters(leftedge+insideslopedown);

% %%%
% figure(192)
% bar(bincenters,fraction),xlim([min(vec) max(vec)]);hold on
% stem(threshlocation,mf,'Color','g')
% %%%%

subtractionthreshold = threshlocation;
    if size(subtractionthreshold,1)==size(subtractionthreshold,2)
    else
         subtractionthreshold = mean(threshlocation);
    end


subtracted = sub_scale_corr-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
Ie(Ie>0)=50;
If = imfill(Ie,'holes');
sldCell.Value = threshlocation;
end




function If = segmentPM(img)
global threshinput sldPM
imgsize = size(img);
ratio = imgsize(1)./2048;


CellimgSmoothingFactor = 30;
sigmafirst = CellimgSmoothingFactor.*3.*ratio;
kernelgsizefirst = CellimgSmoothingFactor.*6.*ratio;

% start
imgorig = img;
sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
subtractionref = double(imgorig);
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
vec = vec(vec>0);
[numbersz,bincenters] = hist(double(vec),min(vec):1:max(vec));

initialsmooth=2;
    numbers = medfilt1(numbersz, initialsmooth); %smooths curve
    fraction = numbers./sum(numbers);
    [~,locs] = findpeaks(fraction);
    threshlocation = bincenters(max(locs));
while length(locs)>4
    numbers = medfilt1(numbersz, initialsmooth); %smooths curve
    fraction = numbers./sum(numbers);
    [~,locs] = findpeaks(fraction);
    threshlocation = bincenters(locs(2));
    initialsmooth=initialsmooth+1;
end

j = imregionalmax(subtractionref);
img = j&(subtractionref>threshlocation);
img = imerode(img,strel('disk',1));



th = min(subtractionref(logical(img)&(subtractionref>threshlocation)));
threshlocation = bincenters(locs(3));

if ~isempty(threshinput.PM)
    threshlocation = threshinput.PM;
end

subtractionthreshold = threshlocation;

subtracted = subtractionref-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);






Ie = subtractedzero;
Ie(Ie>0)=50;
If = imfill(Ie,'holes');
If = imdilate(If,strel('disk',1));
sldPM.Value = threshlocation;
end
function segmentbutton_Callback(~,~) 
global Imagez threshinput
       threshinput.Nuc = []; threshinput.Cell = []; threshinput.PM = [];
       
       Imagez.NucSeg = segmentNucleus(Imagez.NucImage); %semgnet the image
       Imagez.CellSeg = segmentCell(Imagez.CellImage);
       updatePMseg
       updateImage
end

function bw = gaussianBlurz(im,sigma,kernelgsize,varargin)

filtersize = [kernelgsize kernelgsize];
kernelg = fspecial('gaussian',filtersize,sigma);

% image filtering
gFrame = imfilter(im,kernelg,'repl');

if ~isempty(varargin)
    bw=gFrame.*uint16(varargin{1}>0);
else
    bw=gFrame;
end
end

%go forward to next position
function nextpos(~,~)
global ImagePositions Position PositionHandle


% Determine the selected data set.
isitthere = strcmp(ImagePositions,Position);
idx = find(isitthere ==1);
v = PositionHandle.Value;

    if idx ==1
        idx = idx+1;
        v=v+1;
        
    elseif idx == length(ImagePositions)
    %     idx=idx-1;
    else
        idx=idx+1;
        v=v+1;
    %     idx=idx-1;
    end
    
PositionHandle.Value = v;
Position = ImagePositions{idx};
setSceneAndUpdate
end

%go reverse to previous position
function prevpos(~,~)
global ImagePositions Position PositionHandle
% Determine the selected data set.

% Determine the selected data set.

isitthere = strcmp(ImagePositions,Position);
idx = find(isitthere ==1);
v = PositionHandle.Value;

    if idx ==1
%         idx = idx+1; do nothing
    elseif idx == length(ImagePositions)
        idx=idx-1;
        v=v-1;
    else
%         idx=idx+1;
        idx=idx-1;
        v=v-1;
    end
PositionHandle.Value = v;
Position = ImagePositions{idx};
setSceneAndUpdate
end

%add area by drawing
function addareabuttonNUC_Callback(~,~) 
 global Imagez selectedchannel
 
        if strcmp(selectedchannel,'Nucleus')
            str = 'NucSeg';
        elseif strcmp(selectedchannel,'Cell')
            str = 'CellSeg';
        elseif strcmp(selectedchannel,'PM')
            str = 'PMSeg';
        end    
       
        im  = logical(Imagez.(str));
       [polyx,polyy] = ginput();
       M  = zeros(1,length(polyx)*2);
      M(1:2:end) = polyx;
      M(2:2:end) = polyy;
      zeroImage = zeros(size(im));
      zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
      zerogray = rgb2gray(zeroImage);
      im(zerogray>0) = 1;
      
      Imagez.(str) = im;
      
    if strcmp(str,'PMSeg')
        %update cytoplasmic segmentation to exclude new PM segmentation
        cyto = Imagez.CellSeg;
        cyto(im)=0;
        Imagez.CellSeg = cyto;    
    end
      
    updateImage
end

%remove area by drawing
function removeareabuttonNUC_Callback(~,~) 
global Imagez selectedchannel

       if strcmp(selectedchannel,'Nucleus')
            str = 'NucSeg';
        elseif strcmp(selectedchannel,'Cell')
            str = 'CellSeg';
        elseif strcmp(selectedchannel,'PM')
            str = 'PMSeg';
       end    
       
      [polyx,polyy] = ginput();
      

      M  = zeros(1,length(polyx)*2);
      M(1:2:end) = polyx;
      M(2:2:end) = polyy;
      zeroImage = zeros(size(Imagez.(str)));
      zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
      zerogray = rgb2gray(zeroImage);
      Imagez.(str)(zerogray>0) = 0;
      
    if strcmp(str,'PMSeg')
        %update cytoplasmic segmentation to exclude new PM segmentation
        cyto = Imagez.CellSeg;
        cyto(im)=0;
        Imagez.CellSeg = cyto;    
    end
      
    updateImage
end


%delte whole areas
function deletebuttonNUC_Callback(~,~) 
   global Imagez
   stats = regionprops(Imagez.NucSeg,'PixelIdxList');
   % Display mesh plot of the currently selected data.
     [cellxx,cellyy] = ginput();
 
    for j = 1:length(cellxx)  
    celly = cellyy(j);
    cellx = cellxx(j);
        for i = 1:length(stats)
        m = stats(i).PixelIdxList;
        index = sub2ind(size(Imagez.NucSeg),round(celly),round(cellx));
            if ~isempty(find(m==index,1))
            Imagez.NucSeg(m)=0;
            end
        end
    end
updateImage
end

%erode whole areas
function erode_Callback(~,~) 
   global Imagez selectedchannel
   
   
       if strcmp(selectedchannel,'Nucleus')
            str = 'NucSeg';
        elseif strcmp(selectedchannel,'Cell')
            str = 'CellSeg';
        elseif strcmp(selectedchannel,'PM')
            str = 'PMSeg';
       end    
        
   stats = regionprops(Imagez.(str),'PixelIdxList');
   % Display mesh plot of the currently selected data.
     [cellxx,cellyy] = ginput();
   

    
    for j = 1:length(cellxx)  
    celly = cellyy(j);
    cellx = cellxx(j);
        for i = 1:length(stats)
        m = stats(i).PixelIdxList;
        index = sub2ind(size(Imagez.(str)),round(celly),round(cellx));
            if ~isempty(find(m==index,1))
            im = false(size(Imagez.(str)));
            im(m) = 1;
            ime = imerode(im,strel('disk',1));
                if ~isempty(ime>0)
                stats(i).PixelIdxList = find((ime>0)==1);
                end
            Imagez.(str) = ime;
            end
        end
    end
    
    if strcmp(str,'PMSeg')
        %update cytoplasmic segmentation to exclude new PM segmentation
        cyto = Imagez.CellSeg;
        cyto(im)=0;
        Imagez.CellSeg = cyto;    
    end
    
    
updateImage
end

%dilate whole areas
function dilate_Callback(~,~) 
   global Imagez selectedchannel
   
   
       if strcmp(selectedchannel,'Nucleus')
            str = 'NucSeg';
        elseif strcmp(selectedchannel,'Cell')
            str = 'CellSeg';
        elseif strcmp(selectedchannel,'PM')
            str = 'PMSeg';
       end    
        
   stats = regionprops(Imagez.(str),'PixelIdxList');
   % Display mesh plot of the currently selected data.
     [cellxx,cellyy] = ginput();
   

    
    for j = 1:length(cellxx)  
    celly = cellyy(j);
    cellx = cellxx(j);
        for i = 1:length(stats)
        m = stats(i).PixelIdxList;
        index = sub2ind(size(Imagez.(str)),round(celly),round(cellx));
            if ~isempty(find(m==index,1))
            im = false(size(Imagez.(str)));
            im(m) = 1;
            ime = imdilate(im,strel('disk',1));
                if ~isempty(ime>0)
                stats(i).PixelIdxList = find((ime>0)==1);
                end
            Imagez.(str) = ime;
            end
        end
    end
    
    if strcmp(str,'PMSeg')
        %update cytoplasmic segmentation to exclude new PM segmentation
        cyto = Imagez.CellSeg;
        cyto(im)=0;
        Imagez.CellSeg = cyto;    
    end
    
    
updateImage
end

 %save your corrected semgentation
function savebutton_Callback(~,~)
global Imagez Position FluorChannelOne FluorChannelTwo SegmentationFileDirectory

savethatimage(Imagez.NucSeg,Imagez.CellSeg,Imagez.PMSeg,Position,FluorChannelOne,FluorChannelTwo,SegmentationFileDirectory)
end
function savethatimage(NucSeg,CellSeg,PMSeg,Position,FluorChannelOne,FluorChannelTwo,SegmentationFileDirectory)
cd (SegmentationFileDirectory)

NucleusSegFileName = strcat('Segment-',FluorChannelOne,'-',Position,'.mat'); 
CellSegFileName = strcat('Segment-',FluorChannelTwo,'-',Position,'.mat'); 
PMSegFileName = strcat('Segment-','PM','-',Position,'.mat'); 

if isempty(dir(strcat(NucleusSegFileName(1:end-4),'.mat')))
    save(strcat(NucleusSegFileName(1:end-4),'.mat'),'NucSeg');
    else
    save(strcat(NucleusSegFileName(1:end-4),'.mat'),'NucSeg','-append');
end

if isempty(dir(strcat(CellSegFileName(1:end-4),'.mat')))
    save(strcat(CellSegFileName(1:end-4),'.mat'),'CellSeg');
    else
    save(strcat(CellSegFileName(1:end-4),'.mat'),'CellSeg','-append');
end

if isempty(dir(strcat(PMSegFileName(1:end-4),'.mat')))
    save(strcat(PMSegFileName(1:end-4),'.mat'),'PMSeg');
    else
    save(strcat(PMSegFileName(1:end-4),'.mat'),'PMSeg','-append');
end

disp(strcat('image ... ',NucleusSegFileName,'and',CellSegFileName,'...saved'));  
end

