function  uiSegmentTimeLapseImages

global denoisepath  sliderOne sliderOneTxt pStruct subaxes   exportdir  channelinputs adjuster cmapper tcontrast lcontrast   OGExpDate   cmap  A AA timeFrames framesForDir ImageDetails  SceneList  imgsize ExpDate
adjuster=0;
imgsize = [512 512];
tcontrast = 99;
lcontrast = 1;


channelstoinput = {'_mKate','_EGFP','_CFP','DIC','_Hoechst'};
channelinputs =channelregexpmaker(channelstoinput);

clearvars -global SceneDirectoryPath

ImageDetails = InitializeImageDetails;

%%% set colormap for the images %%%
cmap = colormap(gray(255));
% cmap = colormap(magma(255));
% cmap = colormap(inferno(255));
% cmap = colormap(plasma(255));
cmap(255,:)=[1 0 0];
cmapper = cmap;
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose directory of experiment to track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set directory to location of code being used (generally external harddrive
mdir = mfilename('fullpath');
    [~,b ] = regexp(mdir,'/');
        if isempty(b)
            [~,b] = regexp(mdir,'\');
        end
    parentdir = mdir(1:b(end-1));
cd(parentdir)
denoisepath = strcat(parentdir,'Tracking/coherencefilter_version5b/');
addpath(denoisepath);
exportdir = strcat(parentdir,'Tracking\Export\');


%set parent directory
A = uigetdir;
AA = 'D:\Users\zeiss\Documents\MATLAB';
cd(A)
subdirname = 'flatfield_corrected';

%subdirectories should include
%> [ flatfield_corrected ]
    %> [ ####date## smad3g smFISH_scene_s## ]
        %> [ c#_flat ]     [ tiffs ]
            %need to load up the NucleusBinary_flat images

%first determine how many scenes are present
dirlist = dir(subdirname);
[~,~,~,d] = regexp({dirlist.name},'s[0-9]+');
dlog = ~cellfun(@isempty,d,'UniformOutput',1); 
dcell = d(dlog);
SceneList = cellfun(@(x) x{1},dcell,'UniformOutput',0);

%determine date of experiment
cd (subdirname)
[a,b] = regexp(A,'201[0-9]');
ExpDate = A(a:b+6);OGExpDate = ExpDate; [a,~] = regexp(ExpDate,'_');ExpDate(a) = '-';

%determine the number of time frames per scene
folderlist = dir(strcat('*',SceneList{1},'*'));
foldername = folderlist.name;
cd (char(foldername))
spec_directory = '_mKate_flat';
[timeFrames,framesForDir] = determineTimeFrames(spec_directory);
cd .. 


%determine the number of channels
folderlist = dir(strcat('*','*'));
    [~,~,~,channelsListed] = regexp([folderlist.name],channelinputs);
    channelList = unique(channelsListed);
            for i=1:length(channelList)
                chan = channelsListed{i};
                [a,~] = regexp(chan,'_');
                chan(a) = [];
                channelList{i} = chan;
            end
            
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Set up  user interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure;
% f.Visible ='off';
f.Units = 'pixels';
f.Position =[10,10,1800,1000];


buttonwidth = 80;
buttonheight = 50;

ypositions = sort([100:20:1000],'descend');
xpositions = ones(1,length(ypositions)).*1600;

        mmm=1;
% htexttwo 
uicontrol('Style','text','String','To choose channel push 1, 2, or 3',...
          'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight]);
        mmm=mmm+1;mmm=mmm+1;mmm=mmm+1;

% hNextFrame
uicontrol('Style','pushbutton',...
    'String','NextFrame [f]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@nextbutton_callback);
% hPreviousFrame
uicontrol('Style','pushbutton',...
    'String','Previous frame [a]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@prevbutton_callback);
        mmm=mmm+1;mmm=mmm+1;
        
% hGoToFrame
uicontrol('Style','pushbutton',...
    'String','Go to Frame',...
    'Position',[xpositions(mmm)-120,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@gotobutton_callback);

        
% hFinalFrame
uicontrol('Style','pushbutton',...
    'String','FinalFrame [g]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@finalbutton_callback);
% hFirstFrame
uicontrol('Style','pushbutton',...
    'String','First frame [z]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@firstbutton_callback);
% hFirstFrame
uicontrol('Style','pushbutton',...
    'String','saveSomethingCallback',...
    'Position',[xpositions(mmm)-300,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@saveSomethingCallback);
        mmm=mmm+1; mmm=mmm+1;

% htextone
uicontrol('Style','text','String','Choose Scene',...
    'Position',[xpositions(mmm),ypositions(mmm)-buttonheight./2,buttonwidth,buttonheight]);

% htextone
uicontrol('Style','text','String','Choose Channel',...
    'Position',[xpositions(mmm)-200,ypositions(mmm)-buttonheight./2,buttonwidth,buttonheight]);
        mmm=mmm+1;
% hpopup
uicontrol('Style','popupmenu',...
    'String',SceneList',...
    'Position',[xpositions(mmm)-buttonwidth./2,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@popup_menu_Callback);
        
        
% hpopup
uicontrol('Style','popupmenu',...
    'String',channelList',...
    'Position',[xpositions(mmm)-200-buttonwidth./2,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@popup_menu_Callback_channels);
        mmm=mmm+1;mmm=mmm+1; mmm=mmm+1;
        
        
% hAddArea
uicontrol('Style','pushbutton','String','AddArea [v]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@addareabutton_Callback);
% hLinkCells
uicontrol('Style','pushbutton','String','LinkCells [r]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@linkCells_Callback);
        mmm=mmm+1;mmm=mmm+1;
% hDelete
uicontrol('Style','pushbutton','String','Delete [d]',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@deletebutton_Callback);
% hEliminate
uicontrol('Style','pushbutton','String','Eliminate [e]',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@eliminatebutton_Callback);
        mmm=mmm+1;mmm=mmm+1;
        
        
%hDestroy
uicontrol('Style','pushbutton','String','Destroy',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@destroybutton_Callback);

%hchosenOnes 
uicontrol('Style','pushbutton','String','Chosen Ones',...
    'Position',[xpositions(mmm)+40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@chosenOnes_Callback);
       mmm=mmm+1;mmm=mmm+1;
       
% hRemoveArea
uicontrol('Style', 'pushbutton', 'String', 'Remove area',...
    'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight/1.5],...
    'Callback',@removeArea_Callback);
        mmm= mmm+1;
% hchosenOnesEnd 
uicontrol('Style','pushbutton','String','Chosen Ones EndOnly',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*2,buttonheight./2],...
    'Callback',@chosenOnesEnd_Callback);
       mmm=mmm+1;

        
% hErode
uicontrol('Style','pushbutton','String','Erode Chosen Only',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*2,buttonheight./2],...
    'Callback',@erodeOnes_Callback);
       mmm=mmm+1;


% hDilate
uicontrol('Style','pushbutton','String','Dilate Chosen Only',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*2,buttonheight./2],...
    'Callback',@dilateOnes_Callback);
       mmm=mmm+1; mmm=mmm+1;mmm=mmm+1;

        
% hDisplayTracking
uicontrol('Style','pushbutton',...
    'String','DisplayTracking [m]',...
    'Position',[xpositions(mmm)-(buttonwidth./2),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@displayTrackingButton_Callback);
% hTrack
uicontrol('Style','pushbutton',...
    'String','Run Tracking [t]',...
    'Position',[xpositions(mmm)+(buttonwidth./2),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@trackbutton_Callback);
        mmm=mmm+1; mmm=mmm+1;
        mmm=mmm+1;

%  hContrast
uicontrol('Style','pushbutton',...
    'String','contrast user',...
    'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@contrast_Callback);       
        mmm=mmm+1; mmm=mmm+1;      
        
        
        %%%%
        %%%%
% hSaveTrackingAs
uicontrol('Style','pushbutton',...
    'String','SaveTrackingAs',...
    'Position',[xpositions(mmm)-buttonwidth./2,ypositions(mmm)-buttonheight,buttonwidth.*2,buttonheight.*2],...
    'Callback',@saveTrackingFileAs_callback);
        mmm=mmm+1; mmm=mmm+1; mmm=mmm+1; mmm=mmm+1; mmm=mmm+1;

% hLoadTracking
uicontrol('Style','pushbutton',...
    'String','LoadTracking',...
    'Position',[xpositions(mmm),ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@loadTrackingFile_callback);
        mmm=mmm+1; mmm=mmm+1; mmm=mmm+1;

        
% hPlot
uicontrol('Style','pushbutton',...
    'String','PLOT!',...
    'Position',[xpositions(mmm)-0,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@Plot_callback);

% hPlotCFPnorm
uicontrol('Style','pushbutton',...
    'String','plotCFPnorm?',...
    'Position',[xpositions(mmm)-90,ypositions(mmm)+20,buttonwidth,buttonheight./2],...
    'Callback',@PlotCFPnorm_callback);

% hPlotCFPnotnorm 
uicontrol('Style','pushbutton',...
    'String','plotCFPnotnorm?',...
    'Position',[xpositions(mmm)-90,ypositions(mmm),buttonwidth,buttonheight./2],...
    'Callback',@PlotCFPnotnorm_callback);


% hPlotSpecificCell
uicontrol('Style','pushbutton',...
    'String','Plot Specific Cell!',...
    'Position',[xpositions(mmm)+80,ypositions(mmm)+buttonheight./2,buttonwidth,buttonheight./2],...
    'Callback',@Plot_SpecificCell_callback);
% hPlotSettings 
uicontrol('Style','pushbutton',...
    'String','Plot Settings!',...
    'Position',[xpositions(mmm)+80,ypositions(mmm),buttonwidth,buttonheight./2],...
    'Callback',@PlotSettings_callback);
      mmm=mmm+1;  mmm=mmm+1;
      
% hExportCells
uicontrol('Style','pushbutton',...
    'String','Export Cells',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@exportCells);
        mmm=mmm+1; 
        
% hLabelCells
uicontrol('Style','pushbutton',...
    'String','Label Cells',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*1.5,buttonheight./2],...
    'Callback',@labelCells);
        mmm=mmm+1; 
        
        
% hcomment
uicontrol('Style','pushbutton','String','Comments',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth.*2,buttonheight./2],...
    'Callback',@comment_Callback);
       mmm=mmm+1; mmm=mmm+1; mmm=mmm+1;
        
% hExportLabelsCells
uicontrol('Style','pushbutton',...
    'String','ExportLabels',...
    'Position',[xpositions(mmm)-40,ypositions(mmm),buttonwidth,buttonheight],...
    'Callback',@ exportLabels);
%     'Callback',@ exportFrames);
        
       

        
       
f.Visible = 'on'   ;
f.Units = 'normalized';
for i = 1:length(f.Children)
   hhh = f.Children(i);
   hhh.Units = 'normalized';
end
% 
% MainAxes = axes;
% MainAxes.Units = 'pixels';
% MainAxes.XTick=[];
% MainAxes.YTick = [];
% imgdim = 512.*1.8;
% Position = [25 25 imgdim imgdim];
% % Position = [0.1 0.3 0.65 0.65];
% MainAxes.Position = Position;


channelimglength = 9;
xinit = 0.02;
yinit = 0.025;
w = 0.15;
h = 0.15;
xspacefactor = 0.3;
yspacefactor = 1;
x = [xinit xinit+w+(xspacefactor*w) xinit+(w+(xspacefactor*w)).*2 xinit xinit+w+(xspacefactor*w) xinit+(w+(xspacefactor*w)).*2 xinit xinit+w+(xspacefactor*w) xinit+(w+(xspacefactor*w)).*2];
y = fliplr([yinit yinit yinit yinit+h+(yspacefactor*w) yinit+h+(yspacefactor*w) yinit+h+(yspacefactor*w)   yinit+(w+(yspacefactor*w)).*2 yinit+(w+(yspacefactor*w)).*2 yinit+(w+(yspacefactor*w)).*2]);

for i=1:channelimglength
    ax= axes();
    ax.Position = [x(i) y(i) w h];
    ax.Units = 'inches';
    pos = ax.Position;
    pos(4) = pos(3);
    ax.Position = pos;
    ax.Units = 'normalized';
    ax.XTick = [];
    ax.YTick = [];
    subaxes(i) = ax;
end





%parameter structure
pStruct = struct();
parameterDefaults = [106 106 40;...
                    0.5 0.5 1.3;...
                    20 20 5;...
                    1 1 1];
                    
parameterStrings = {'nucDiameter','threshFactor','sigmaScaledToParticle','noparametercurrently'};

for p = 1:length(parameterStrings)
    pString = char(parameterStrings{p});
    for c = 1:length(channelList)
        cString = char(channelList{c});
        pStruct.(cString).(pString) = parameterDefaults(p,c); 
    end
end
% str = 'nucDiameter';
% channel = 'mKate';
% val.(channel).(str) = 60;
% pStruct.EGFP.(str)=val.(str);
% pStruct.mKate.(str)=val.(str);
% pStruct.CFP.(str)=val.(str);
% 
% 
% str = 'threshFactor';
% pStruct.EGFP.(str)=val.(str);
% pStruct.mKate.(str)=val.(str);
% pStruct.CFP.(str)=val.(str);
% 
% str = 'sigmaScaledToParticle'; 
% pStruct.EGFP.(str)=val.(str);
% pStruct.mKate.(str)=val.(str);
% pStruct.CFP.(str)=val.(str);
% 
% str = 'noparemtercurrently';
% pStruct.EGFP.(str)=val.(str);
% pStruct.mKate.(str)=val.(str);
% pStruct.CFP.(str)=val.(str);
% 
% nucDiameter=122;
% threshFactor = 0.94;
% sigmScaledToParticle = 2;
% noparametercurrently=60;

% ImageDetails.Channel = 'EGFP';

set(f,'KeyPressFcn',@keypress);

end








function updateSliders
global pStruct ImageDetails sliderOne sliderOneText

%slider1
sliderx = 0.72;
% slidery = 0.6;
sliderw = 0.1;
sliderh = 0.05;
slidertextw = 0.05;
sliderspace = 0.1;





channel = ImageDetails.Channel;

%nucDiameter
str = 'nucDiameter';
val.(str) = pStruct.(channel).(str);
slidery = 0.55;
minz = 1;
maxz = 400;
sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./(maxz-minz),'Value',val.(str),'Position', [1 1 1 1],...
        'Callback', @sliderOneAdjust); 
    sliderOne.(str).Units='normalized';
    sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str(uint8(val.(str)))));
    sliderOneTxt.(str).Units= 'Normalized';
    sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];  

    
    
    
% threshFactor
str = 'threshFactor';
val.(str) = pStruct.(channel).(str);
slidery = 0.5;
minz = 0.5;
maxz = 3;
sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./(10*(maxz-minz)),'Value',val.(str),'Position', [1 1 1 1],...
        'Callback', @sliderOneAdjust); 
    sliderOne.(str).Units='normalized';
    sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str((val.(str)))));
    sliderOneTxt.(str).Units= 'Normalized';
    sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];  

    
   
    
    
% sigmaScaledToParticle (the divide diameter by this factor to get sigma for
% gaussian smoothing)
str = 'sigmaScaledToParticle';
val.(str) = pStruct.(channel).(str);
slidery = 0.45;
minz = 1;
maxz = 40;
sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./(maxz-minz),'Value',val.(str),'Position', [1 1 1 1],...
        'Callback', @sliderOneAdjust); 
    sliderOne.(str).Units='normalized';
    sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str(uint8(val.(str)))));
    sliderOneTxt.(str).Units= 'Normalized';
    sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];   
  
    
    
% noparemtercurrently
str = 'noparametercurrently';
val.(str) = pStruct.(channel).(str);
slidery = 0.40;
minz = 1;
maxz = 400;
sliderOne.(str) = uicontrol('Style', 'slider','String',str,'Min',minz,'Max',maxz,'SliderStep',[1 1]./(maxz-minz),'Value',val.(str),'Position', [1 1 1 1],...
        'Callback', @sliderOneAdjust); 
    sliderOne.(str).Units='normalized';
    sliderOne.(str).Position = [sliderx slidery sliderw sliderh];
sliderOneTxt.(str) = uicontrol('Style','text','Units','Normalized','Position',[1 1 1 1],'String',strcat(str,'=',num2str(uint8(val.(str)))));
    sliderOneTxt.(str).Units= 'Normalized';
    sliderOneTxt.(str).Position = [sliderx-sliderspace slidery slidertextw sliderh];   

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






function sliderOneAdjust(source,~)
    global pStruct  ImageDetails sliderOneTxt
    channel = ImageDetails.Channel;
    str = source.String;
%     threshinput.(str) =source.Value;
%     zerostrel = round(source.Value);

    if strcmpi(str,'threshFactor')
        valupdate = source.Value;  
    else
        valupdate = round(source.Value);
    end
    source.Value = valupdate;
    pStruct.(channel).(str) = valupdate; 
    disp(valupdate)
   source.Visible = 'off';
   sliderOneTxt.(str).String = 'waiting...';
   pause(0.001);
    setSceneAndTime
    disp('done')
    source.Visible = 'on';
    
    sliderOneTxt.(str).String = strcat(str,'=',num2str(pStruct.(channel).(str)));
    
end









function FinalImage=loadStack(FileTif)
% [a,b] = uigetfile;
% FileTif = a;
% cd (b)
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   FinalImage(:,:,i)=TifLink.read();
end
TifLink.close();
end




function IfFinal = segmentationNucleus(FinalImage,subdirname,scenename,filename,channel)
global  pStruct 
cd(subdirname)

foldername = '_mKate_flat';
tsn = determineTimeFrame(foldername);

% channel ='mKate';
% parameters
nucDiameter = pStruct.(channel).nucDiameter;
threshFactor = pStruct.(channel).threshFactor;
sigmaScaledToParticle = pStruct.(channel).sigmaScaledToParticle;
finalerode=2;



% prepareCcodeForAnisotropicDiffusionDenoising(denoisepath)

%start
for frames = 1:size(FinalImage,3)
%Smooth Image using Anisotropic Diffusion
% Options.Scheme :  The numerical diffusion scheme used
%                     'R', Rotation Invariant, Standard Discretization 
%                          (implicit) 5x5 kernel (Default)
%                     'O', Optimized Derivative Kernels
%                     'I', Implicit Discretization (only works in 2D)
%                     'S', Standard Discretization
%                     'N', Non-negativity Discretization
%   Options.T  :      The total diffusion time (default 5)
%   Options.dt :      Diffusion time stepsize, in case of scheme H,R or I
%                     defaults to 1, in case of scheme S or N defaults to
%                     0.15. 
%   Options.sigma :   Sigma of gaussian smoothing before calculation of the
%                     image Hessian, default 1.                   
%   Options.rho :     Rho gives the sigma of the Gaussian smoothing of the 
%                     Hessian, default 1.
%   Options.verbose : Show information about the filtering, values :
%                     'none', 'iter' (default) , 'full'
%   Options.eigenmode : There are many different equations to make an diffusion tensor,
%						this value (only 3D) selects one.
%					    0 (default) : Weickerts equation, line like kernel
%						1 : Weickerts equation, plane like kernel
%						2 : Edge enhancing diffusion (EED)
%						3 : Coherence-enhancing diffusion (CED)
%						4 : Hybrid Diffusion With Continuous Switch (HDCS)
    img = FinalImage(:,:,frames); 
    imgRaw = img;
        Options.T = 5;
        Options.dt = 1;
        Options.Scheme = 'R';
        Options.rho = 5;
        Options.sigma = 5;
        Options.verbose = 'none';
%     imgRawDenoised = CoherenceFilter(imgRaw, Options);
imgRawDenoised = imgRaw;

    
    %Based on algorithm of Fast and accurate automated cell boundary determination for fluorescence microscopy by Arce et al (2013)   
    %LOW PASS FILTER THE IMAGE (scale the gaussian filter to diameter of
    %nuclei -- diameter of nuclei is about 50 to 60))
    kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
    sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize
    imgLowPass = gaussianBlurz(double(imgRawDenoised),sigma,kernelgsize);
    rawMinusLP = double(imgRawDenoised) -double(imgLowPass);%%%%%%% key step!
    rawMinusLPvec = reshape(rawMinusLP,size(rawMinusLP,1)^2,1);
    globalMinimaValues = prctile(rawMinusLPvec,0.01);
    globalMinimaIndices = find(rawMinusLP < globalMinimaValues);
    LPscalingFactor = imgRawDenoised(globalMinimaIndices)./imgLowPass(globalMinimaIndices);
    imgLPScaled = imgLowPass.*nanmedian(LPscalingFactor);
    rawMinusLPScaled = double(imgRawDenoised) - double(imgLPScaled);


    %determine the threshold by looking for minima in log-scaled histogram
    %of pixels from rawMinusLPScaled
    rawMinusLPScaledContrasted = imadjust(uint16(rawMinusLPScaled));
    vecOG = double(reshape(rawMinusLPScaledContrasted,size(rawMinusLPScaledContrasted,1)^2,1));
    logvecpre = vecOG; logvecpre(logvecpre==0)=[];
    logvec = log10(logvecpre);
    vec = logvec;
    [numbers,bincenters] = hist(vec,prctile(vec,1):(prctile(vec,99)-prctile(vec,1))/1000:prctile(vec,99));
    numbersone = medfilt1(numbers, 10); %smooths curve
    numberstwo = medfilt1(numbersone, 100); %smooths curve
    fraction = numberstwo./sum(numberstwo);
    mf = max(fraction);
        %%%%%%%%%%%%%%%%%%%% Important parameters for finding minima of
        %%%%%%%%%%%%%%%%%%%% histogram
        left=0.5*mf;
        slopedown=0.1*mf;
        %%%%%%%%%%%%%%%%%%%%%
    leftedge = find(fraction > left,1,'first');
    insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
    threshLocation = bincenters(leftedge+insideslopedown);
    subtractionThreshold = threshLocation;

    if size(subtractionThreshold,1)==size(subtractionThreshold,2)
        else
         subtractionThreshold = mean(threshLocation);
    end


    subtractionThresholdScaled = (10.^subtractionThreshold).*threshFactor;
    subtracted = double(rawMinusLPScaledContrasted)-subtractionThresholdScaled;
    subzero = (subtracted<0);
    Ih = ~subzero;

%     Ihd = imdilate(Ih,strel('disk',1));
%     Ihdc = imclose(Ihd,strel('disk',2));
%     Ihdcf = imfill(Ihdc,'holes');
%     Im = Ihdcf;
    Ihc = imclose(Ih,strel('disk',4));
    Ihcf = imfill(Ihc,'holes');
    Im=Ihcf;



    %%%%% this is the ultimate addition for watershed segmentation!!!
    see = strel('disk',1);
    seo = strel('disk',8);
    Isum = Im;
    Ier = Isum;
%     figure(2)
    for i=1:(nucDiameter/2)
        Ier = imerode(Ier,see);
%                 subplot(1,2,1);
%                 imagesc(Ier)
%         Iero = imopen(Ier,seo);
    %     Isum = Isum+(Iero.*i);
        Isum = Isum+Ier;
%                 subplot(1,2,2);
%                 imagesc(Isum)
%                 drawnow
%                 pause(0.1)
%         Ier=Iero;
%  disp(max(max(Isum)))
    end
    Isum(Isum>nucDiameter) = nucDiameter;

    gausshed = gaussianBlurz(Isum,ceil(sigma./4),ceil(kernelgsize./1));
    imgt = -double(gausshed);
    waterBoundary = imerode(Im,strel('disk',1));
    imgt(~(waterBoundary>0)) = -Inf;
    L=watershed(imgt);

    L(waterBoundary<1) = 0;
    If = L>1;
%     If = imerode(If,strel('disk',finalerode));


    
    time = tsn{frames};
    tim = time(2:end);
    IfFinal(:,:,frames)=If;
end


stophere=1;
end

function IfFinal = segmentationEGFP(FinalImage,subdirname,scenename,filename,channel)
global  pStruct 
cd(subdirname)

foldername = '_mKate_flat';
tsn = determineTimeFrame(foldername);

% parameters
nucDiameter = pStruct.(channel).nucDiameter;
threshFactor = pStruct.(channel).threshFactor;
sigmaScaledToParticle = pStruct.(channel).sigmaScaledToParticle;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize

finalerode=2;



% prepareCcodeForAnisotropicDiffusionDenoising(denoisepath)

%start
for frames = 1:size(FinalImage,3)
%Smooth Image using Anisotropic Diffusion
% Options.Scheme :  The numerical diffusion scheme used
%                     'R', Rotation Invariant, Standard Discretization 
%                          (implicit) 5x5 kernel (Default)
%                     'O', Optimized Derivative Kernels
%                     'I', Implicit Discretization (only works in 2D)
%                     'S', Standard Discretization
%                     'N', Non-negativity Discretization
%   Options.T  :      The total diffusion time (default 5)
%   Options.dt :      Diffusion time stepsize, in case of scheme H,R or I
%                     defaults to 1, in case of scheme S or N defaults to
%                     0.15. 
%   Options.sigma :   Sigma of gaussian smoothing before calculation of the
%                     image Hessian, default 1.                   
%   Options.rho :     Rho gives the sigma of the Gaussian smoothing of the 
%                     Hessian, default 1.
%   Options.verbose : Show information about the filtering, values :
%                     'none', 'iter' (default) , 'full'
%   Options.eigenmode : There are many different equations to make an diffusion tensor,
%						this value (only 3D) selects one.
%					    0 (default) : Weickerts equation, line like kernel
%						1 : Weickerts equation, plane like kernel
%						2 : Edge enhancing diffusion (EED)
%						3 : Coherence-enhancing diffusion (CED)
%						4 : Hybrid Diffusion With Continuous Switch (HDCS)
    img = FinalImage(:,:,frames); 
    imgRaw = gaussianBlurz(double(img),ceil(sigma./10),ceil(kernelgsize./10));

    imgW = wiener2(img,[1 20]);
    imgWW = wiener2(imgW,[20 1]);
    imgWWW = wiener2(imgWW,[5 5]);
    imgRawDenoised = imgWWW;
    denoiseVec = double(reshape(imgRawDenoised,size(imgRawDenoised,1)^2,1));
    highpoints = prctile(denoiseVec,95);
    imgRawDenoised(imgRawDenoised>highpoints) = highpoints;
    
%         Options.T = 5;
%         Options.dt = 1;
%         Options.Scheme = 'R';
%         Options.rho = 20;
%         Options.sigma = 20;
%         Options.verbose = 'none';ii
% %     imgRawDenoised = CoherenceFilter(imgRaw, Options);
% % imgRawDenoised = imgRaw;

    
    %Based on algorithm of Fast and accurate automated cell boundary determination for fluorescence microscopy by Arce et al (2013)   
    %LOW PASS FILTER THE IMAGE (scale the gaussian filter to diameter of
    %nuclei -- diameter of nuclei is about 50 to 60))
    
    imgLowPass = gaussianBlurz(double(imgRawDenoised),sigma,kernelgsize);
    rawMinusLP = double(imgRawDenoised) -double(imgLowPass);%%%%%%% key step!
    rawMinusLPvec = reshape(rawMinusLP,size(rawMinusLP,1)^2,1);
    globalMinimaValues = prctile(rawMinusLPvec,0.01);
    globalMinimaIndices = find(rawMinusLP < globalMinimaValues);
    LPscalingFactor = imgRawDenoised(globalMinimaIndices)./imgLowPass(globalMinimaIndices);
    imgLPScaled = imgLowPass.*nanmedian(LPscalingFactor);
    rawMinusLPScaled = double(imgRawDenoised) - double(imgLPScaled);


    %determine the threshold by looking for minima in log-scaled histogram
    %of pixels from rawMinusLPScaled
    rawMinusLPScaledContrasted = imadjust(uint16(rawMinusLPScaled));
    vecOG = double(reshape(rawMinusLPScaledContrasted,size(rawMinusLPScaledContrasted,1)^2,1));
    logvecpre = vecOG; logvecpre(logvecpre==0)=[];
    logvec = log10(logvecpre);
    vec = logvec;
    [numbers,bincenters] = hist(vec,prctile(vec,1):(prctile(vec,99)-prctile(vec,1))/1000:max(vec));
    numbersone = medfilt1(numbers, 10); %smooths curve
    numberstwo = medfilt1(numbersone, 100); %smooths curve
    fraction = numberstwo./sum(numberstwo);
    mf = max(fraction);
        %%%%%%%%%%%%%%%%%%%% Important parameters for finding minima of
        %%%%%%%%%%%%%%%%%%%% histogram
        left=0.5*mf;
        slopedown=0.4*mf;
        %%%%%%%%%%%%%%%%%%%%%
    leftedge = find(fraction > left,1,'first');
    insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
    threshLocation = bincenters(leftedge+insideslopedown-1);
    subtractionThreshold = threshLocation;

    if size(subtractionThreshold,1)==size(subtractionThreshold,2)
        else
         subtractionThreshold = mean(threshLocation);
    end


    subtractionThresholdScaled = (10.^subtractionThreshold).*threshFactor;
    subtracted = double(rawMinusLPScaledContrasted)-subtractionThresholdScaled;
    subzero = (subtracted<0);
    Ih = ~subzero;

%     Ihd = imdilate(Ih,strel('disk',1));
%     Ihdc = imclose(Ihd,strel('disk',2));
%     Ihdcf = imfill(Ihdc,'holes');
%     Im = Ihdcf;
    Ihc = imclose(Ih,strel('disk',20));
    Ihcf = imfill(Ihc,'holes');
%     Ihcf = Ihc;
    Ihcfd = imdilate(Ihcf,strel('disk',20));
%     Ihcfd = Ihcf;
    Im=Ihcfd;
    If =Im;





    
    time = tsn{frames};
    tim = time(2:end);
    IfFinal(:,:,frames)=If;
end


stophere=1;
end


function IfFinal = segmentationNucleusRambling(FinalImage,subdirname,scenename,filename,channel)
global  pStruct
cd(subdirname)

dirlist = dir('_mKate_flat');
foldername = '_mKate_flat';
tsn = determineTimeFrame(foldername);

channel ='mKate';
% parameters

nucDiameter = pStruct.mKate.nucDiameter;
threshFactor = pStruct.mKate.threshFactor;
sigmaScaledToParticle = pStruct.mKate.sigmaScaledToParticle;
noparemtercurrently = pStruct.mKate.noparemtercurrently;

finalerode=2;




%start
for frames = 1:size(FinalImage,3)
    img = FinalImage(:,:,frames); 
    imgorig = img;

     %find circles works quite well
%     figure(1029)
%     imagesc(imgorig);
%     [centers, radii] = imfindcircles(imgorig,[20 50],'ObjectPolarity','bright');
%     h=viscircles(centers,radii);
%     stophere=1;


% Options.Scheme :  The numerical diffusion scheme used
%                     'R', Rotation Invariant, Standard Discretization 
%                          (implicit) 5x5 kernel (Default)
%                     'O', Optimized Derivative Kernels
%                     'I', Implicit Discretization (only works in 2D)
%                     'S', Standard Discretization
%                     'N', Non-negativity Discretization
%   Options.T  :      The total diffusion time (default 5)
%   Options.dt :      Diffusion time stepsize, in case of scheme H,R or I
%                     defaults to 1, in case of scheme S or N defaults to
%                     0.15. 
%   Options.sigma :   Sigma of gaussian smoothing before calculation of the
%                     image Hessian, default 1.                   
%   Options.rho :     Rho gives the sigma of the Gaussian smoothing of the 
%                     Hessian, default 1.
%   Options.verbose : Show information about the filtering, values :
%                     'none', 'iter' (default) , 'full'
%   Options.eigenmode : There are many different equations to make an diffusion tensor,
%						this value (only 3D) selects one.
%					    0 (default) : Weickerts equation, line like kernel
%						1 : Weickerts equation, plane like kernel
%						2 : Edge enhancing diffusion (EED)
%						3 : Coherence-enhancing diffusion (CED)
%						4 : Hybrid Diffusion With Continuous Switch (HDCS)
    figure(22)
    subplot(1,4,1);
    imagesc(imgorig);

    Options.T = 5;
    Options.dt = 1;
    Options.Scheme = 'R';
    Options.rho = 5;
    Options.sigma = 5;
    Options.verbose = 'none';
    imgdiff = CoherenceFilter(imgorig, Options);
    subplot(1,4,2);
    imagesc(imgdiff);
       
    %diameter of nuclei is about 50 to 60
    kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
    sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize
    imgLowPass = gaussianBlurz(double(imgdiff),sigma,kernelgsize);

    subplot(1,4,3);
    imagesc(imgLowPass)

    sub = double(imgdiff) -double(imgLowPass);%%%%%%% key step!
    subvec = reshape(sub,size(sub,1)^2,1);
    bb = prctile(subvec,0.01);
    b = find(sub < bb);
    imgLowPassMinusOrigDiff = imgdiff(b)./imgLowPass(b);
    gausscaled = imgLowPass.*nanmedian(imgLowPassMinusOrigDiff);
    sub_scale_corr = double(imgdiff) - double(gausscaled);
    
    subplot(1,4,4);
    imagesc(sub_scale_corr);


    subtractionref = imadjust(uint16(sub_scale_corr));
    vecOG = double(reshape(subtractionref,size(subtractionref,1)^2,1));
    logvecpre = vecOG; logvecpre(logvecpre==0)=[];
    logvec = log10(logvecpre);
    vec = logvec;
    [numbers,bincenters] = hist(vec,prctile(vec,1):(prctile(vec,99)-prctile(vec,1))/1000:prctile(vec,99));


    numberstwo = numbers;
    numbersnorm = numberstwo./sum(numbers);
    numbers = medfilt1(numberstwo, 10); %smooths curve
    numbers = medfilt1(numbers, 100); %smooths curve
    fraction = numbers./sum(numbers);
    hold off
    plot(numbersnorm)
    hold on
    plot(fraction)
    

    mf = max(fraction);

    %%%%%%%%%%%%%%%%%%%%
    left=0.5*mf;
    slopedown=0.1*mf;
    %%%%%%%%%%%%%%%%%%%%%

    leftedge = find(fraction > left,1,'first');
    insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
    insideslopeup = find(fraction(leftedge+insideslopedown:end) >0.0012,1,'first');
    trough = min(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup));
    troughindex = find(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup) == trough);


    threshlocation = bincenters(leftedge+insideslopedown);

    figure(242)
        bar(bincenters,fraction);hold on
        xlim([prctile(vec,1) prctile(vec,99)])
        ylim([0 mf.*2])
    stem ([threshlocation threshlocation],[0 1]);hold off
    drawnow

    subtractionthreshold = threshlocation;

    if size(subtractionthreshold,1)==size(subtractionthreshold,2)
        else
         subtractionthreshold = mean(threshlocation);
    end


    subtractionThresholdScaled = (10.^subtractionthreshold).*threshFactor;
    subtracted = double(subtractionref)-subtractionThresholdScaled;
%     subtracted = double(subtractionref)-subtractionthreshold;
    subzero = (subtracted<0);
    subtractedzero = subtracted.*(~subzero);


    Ie = subtractedzero;
    Ie(Ie>0) = 1;

    Ih = logical(Ie);
    Ihd = imdilate(Ih,strel('disk',1));
    Ihdc = imclose(Ihd,strel('disk',2));
    Ihdcf = imfill(Ihdc,'holes');
    Im = Ihdcf;
%     Im=Ih;



    %%%%% this is the ultimate addition for watershed segmentation!!!
    see = strel('disk',1);
    seo = strel('disk',8);
    Isum = Im;
    Ier = Isum;
%     figure(2)
    for i=1:(nucDiameter./1.5)
        Ier = imerode(Ier,see);
%                 subplot(1,2,1);
%                 imagesc(Ier)
%         Iero = imopen(Ier,seo);
    %     Isum = Isum+(Iero.*i);
        Isum = Isum+Ier;
%                 subplot(1,2,2);
%                 imagesc(Isum)
%                 drawnow
%                 pause(0.1)
%         Ier=Iero;
%  disp(max(max(Isum)))
    end
    Isum(Isum>nucDiameter) = nucDiameter;

    gausshed = gaussianBlurz(Isum,sigma./4,kernelgsize./4);
    imgt = -double(gausshed);
    waterBoundary = imerode(Im,strel('disk',1));
    imgt(~(waterBoundary>0)) = -Inf;
    L=watershed(imgt);

    L(waterBoundary<1) = 0;
    If = L>1;
    If = imerode(If,strel('disk',finalerode));


    
    
%  %exclude non-round objects
%     [B,L] = bwboundaries(If,'noholes');
%     Cells = bwconncomp(If);
%     PX = Cells.PixelIdxList;
%     Ifbefore = If;
%     
%     stats = regionprops(L,'Area','Centroid');
%         threshold = 0.6;
%         % loop over the boundaries
%         pxlog = false(length(PX));
%         for k = 1:length(B)
%           % obtain (X,Y) boundary coordinates corresponding to label 'k'
%           boundary = B{k};
%           % compute a simple estimate of the object's perimeter
%           delta_sq = diff(boundary).^2;
%           perimeter = sum(sqrt(sum(delta_sq,2)));
%           % obtain the area calculation corresponding to label 'k'
%           area = stats(k).Area;
%           % compute the roundness metric
%           metric = 4*pi*area/perimeter^2;
%           
%           if metric<threshold
%               If(PX{k})=0;
%           end
%           
%           % display the results
% %           metric_string = sprintf('%2.2f',metric);
%     
%           % mark objects above the threshold with a black circle
% %                 figure(232)
% %           if metric > threshold
% %             centroid = stats(k).Centroid;
% %             plot(centroid(1),centroid(2),'ko');hold on
% %           end
% % 
% %           text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','k',...
% %                'FontSize',14,'FontWeight','bold');
% 
%         end
% % 
% %         title(['Metrics closer to 1 indicate that ',...
% %                'the object is approximately round']);
% 
      
% % Ifbefore = If;
% figure(33)
% i = imagesc(Ifbefore);
% ax = i.Parent;
% ax.YDir = 'normal';
% 
% figure(44)
% i = imagesc(If);
% ax = i.Parent;
% ax.YDir = 'normal';


    time = tsn{frames};
    tim = time(2:end);
%     if frames==10
%         figure(1)
%         imagesc(If)
%         stophere=1;
%     end

    % savethatimage(scenename,time,If.*255,frames,filename,channel)
    IfFinal(:,:,frames)=If;
end


stophere=1;
end


function If = segmentationMNG(FinalImage,subdirname,scenename,filename,channel)
global nucleus_seg
fig=1;
mkdir(strcat(channel));

% parameters
left = 0.004;
slopedown = 0.003;

dimdiff = 2048./size(FinalImage(:,:,1),1);

zerostrel = 2;
firststrel = round(30./dimdiff);
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;
fracsmoothing = 0.5.*dimdiff;

dirlist = dir(nucleus_seg);
if isempty(dirlist)
    dirlist = dir('mKate_flat');
    foldername = 'mKate_flat';
else
    foldername = nucleus_seg;
end

tsn = determineTimeFrame(foldername);


% start
for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 
imgorig = img;

img = wiener2(img,[5 5]);
se =strel('disk',zerostrel);
Ie = imerode(imgorig,se);
Iobr = imreconstruct(Ie,img);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
IobrcbrF = imcomplement(Iobrcbr);
gaus = double(IobrcbrF);

se =strel('disk',firststrel);
Ie = imerode(gaus,se);
Iobr = imreconstruct(Ie,gaus);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
gaus = double(Iobrcbr);


sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
gaustwo = gaussianBlurz(double(gaus),sigma,kernelgsize);

sub = double(gaus) -double(gaustwo);%%%%%%% key step!
b = find(sub == min(min(sub)),1,'first');
rattio = gaustwo(b)./gaus(b);
gaustwocorr = gaustwo./rattio;
sub_scale_corr = double(gaus) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:fracsmoothing:10000);



numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);

%%%%%%%%%%%%%%%%%%%%
left=0.3*mf;
slopedown=0.2*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
insideslopeup = find(fraction(leftedge+insideslopedown:end) >0.0012,1,'first');
trough = min(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup));
troughindex = find(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup) == trough);
troughindexrounded = round(median(troughindex));



threshlocation = bincenters(leftedge+insideslopedown);
% 
% figure(22)
%     bar(bincenters,fraction);hold on
%     xlim([-500 1000])
%     ylim([0 0.1])
% stem ([threshlocation threshlocation],[0 1]);hold off
% drawnow

subtractionthreshold = threshlocation;

if size(subtractionthreshold,1)==size(subtractionthreshold,2)
else
     subtractionthreshold = mean(threshlocation);
end

% subtractionthreshold = graythresh(subtractionref);

subtracted = sub_scale_corr-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;





Ih = Ie>0;
Igclose = imclose(Ih,strel('disk',round(30./dimdiff)));
Igclosemax = imclose(Ih,strel('disk',round(80./dimdiff)));
Igcopenmax = imopen(Igclosemax,strel('disk',round(10./dimdiff)));
Igcopen = imopen(Igclose,strel('disk',2));
Igcofill = imfill(Igcopen,'holes');
Igcfopen = bwareaopen(Igcofill,round(5000./(dimdiff.^2)));
Igcfopendil = imerode(Igcfopen,strel('disk',round(5)));

% finalsigma=20;
% finalkernelgsize=40;
gaus = gaussianBlurz(IobrcbrF,round(sigmafirst./dimdiff),(kernelgsizefirst./dimdiff));

% sigma=40;
% kernelgsize=80;
% gaus = gaussianBlurz(gaus,sigmafirst.*2,kernelgsizefirst.*2);

imgt = -double(gaus);
% imgt(~(Igcfopen>0)) = -Inf;
imgt(~(Igcopenmax>0)) = -Inf;

L=watershed(imgt);

L(Igcfopendil<1) = 0;
% imagesc(L)
% colormap parula
If = L>1;

stophere=1;
% time = settimecharacter(frames);
time = tsn{frames};
tim = time(2:end);
savethatimage(scenename,time,If.*255,frames,filename,channel)
end


end

function If = segmentationBKGsecond(FinalImage,subdirname,scenename,filename,channel)
global nucleus_seg
fig=1;
mkdir(strcat(channel));
% parameters
dimdiff = 2048./size(FinalImage(:,:,1),1);

zerostrel = 2;
% firststrel = round(30./dimdiff);
firststrel = round(20./dimdiff);
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;
fracsmoothing = 0.5.*dimdiff;

%set the time
dirlist = dir(nucleus_seg);
if isempty(dirlist)
    dirlist = dir('mKate_flat');
    foldername = 'mKate_flat';
else
    foldername = nucleus_seg;
end

tsn = determineTimeFrame(foldername);

% start
for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 
imgorig = img;

img = wiener2(img,[5 5]);
se =strel('disk',zerostrel);
Ie = imerode(imgorig,se);
Iobr = imreconstruct(Ie,img);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
IobrcbrF = imcomplement(Iobrcbr);
gaus = double(IobrcbrF);

se =strel('disk',firststrel);
Ie = imerode(gaus,se);
Iobr = imreconstruct(Ie,gaus);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
gaus = double(Iobrcbr);


sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
gaustwo = gaussianBlurz(double(gaus),sigma,kernelgsize);

sub = double(gaus) -double(gaustwo);%%%%%%% key step!
b = find(sub == min(min(sub)),1,'first');
rattio = gaustwo(b)./gaus(b);
gaustwocorr = gaustwo./rattio;
sub_scale_corr = double(gaus) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:fracsmoothing/10:10000);
numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);
%%%%%%%%%%%%%%%%%%%%
left=0.5*mf;
slopedown=0.4*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
threshlocation = bincenters(leftedge+insideslopedown);

% figure(22)
%     bar(bincenters,fraction./mf);hold on
%     xlim([-500 1000])
%     ylim([0 0.8])
% stem ([threshlocation threshlocation],[0 10],'g');hold off
% drawnow

subtractionthreshold = threshlocation;

if size(subtractionthreshold,1)==size(subtractionthreshold,2)
else
     subtractionthreshold = mean(threshlocation);
end

% subtractionthreshold = graythresh(subtractionref);

subtracted = sub_scale_corr-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;





Ih = Ie>0;
Igclose = imclose(Ih,strel('disk',round(30./dimdiff)));
Igclosemax = imclose(Ih,strel('disk',round(80./dimdiff)));
Igcopenmax = imopen(Igclosemax,strel('disk',round(10./dimdiff)));
Igcopen = imopen(Igclose,strel('disk',2));
Igcofill = imfill(Igcopen,'holes');
Igcfopen = bwareaopen(Igcofill,round(5000./(dimdiff.^2)));
Igcfopendil = imerode(Igcfopen,strel('disk',round(5)));

% finalsigma=20;
% finalkernelgsize=40;
% gaus = gaussianBlurz(IobrcbrF,round(sigmafirst./dimdiff),(kernelgsizefirst./dimdiff));
% 
% % sigma=40;
% % kernelgsize=80;
% % gaus = gaussianBlurz(gaus,sigmafirst.*2,kernelgsizefirst.*2);
% 
% imgt = -double(gaus);
% % % imgt(~(Igcfopen>0)) = -Inf;
% % imgt(~(Igcopenmax>0)) = -Inf;
% imgt(~(Igcopenmax>0)) = 0;
% 
% % L=watershed(imgt);
% L=imgt;
L = Igcopenmax;

% L(Igcfopendil<1) = 0;
% imagesc(L)
% colormap parula
If = L>0;
% se = strel('disk',16);
se = strel('disk',12);
ifd = imdilate(If,se);
ifdc = imclose(ifd,se);
    ifdc(1:512,1)=255;
    ifdc(1:512,512)=255;
    ifdc(1,1:512)=255;
    ifdc(512,1:512)=255;
    

If = ifdc;
cd('NucleusBinary_flat')
time = tsn{frames};
tim = time(2:end);


% time = settimecharacter(frames);
% tim = time(1:end);
filelist = dir(strcat('*t',tim,'_*'));
% if isempty(filelist)
%     stophere=1;
%     filelist = dir(strcat('*t',tim(2:end),'_*'));
% end
clog = imread(char(filelist.name));

% se = strel('disk',4);
se = strel('disk',8);
clogd = imdilate(clog,se);
If(logical(clogd))=1;
cd .. 

prim = bwperim(If);
se = strel('disk',2);
primd = imdilate(prim,se);
imgorig(primd) = max(max(imgorig));
% figure(99),imagesc(imgorig)
% drawnow
stophere=1;
savethatimage(scenename,time,If.*255,frames,filename,channel)
end


end


function savethatimage(scenename,time,Ie,frames,filename,channel)
cd(channel)

[a,b] = regexp(filename,'(_mKate|CFP|EGFP|DIC)');
fname = strcat(filename(1:a-1),channel,filename(b+1:end));
fname = filename;

[a,b] = regexp(fname,'_t[0-9]+');
fname(a:b) = strcat('_',time);
% % % if isempty(a)
% % %     [a,b] = regexp(fname,'t[0-9][0-9]');
% % %     tnum = str2double(fname(a+1:b))-1;
% % %     tm = round(tnum+str2double(time));
% % %     time = settimecharacter(tm);
% % %     if (length(time)>length(fname(a+1:b)))
% % %     time(1) = 't';
% % %     else
% % %     time = horzcat('t',time);
% % % %     disp(strcat('time',time));
% % %     end
% % %     fname(a:b) = time;
% % % else
% % %     time = horzcat('t',time);
% % % %     disp(strcat('2time',time));
% % %     fname(a:b) = time;
% % % end
% if length(time)==length(fname(a+1:b))
% time = horzcat('t',time);
% fname(a:b) = time;
% elseif length(fname(a+1:b))<length(time)
% time = horzcat('t',time(2:end));
% fname(a:b) = time;
% end

imwrite(uint8(Ie),fname,'tiff','WriteMode','overwrite');
cd .. 

% cd('tiffs')
% % imwrite(uint8(Ie),strcat(scenename,'_','t',time,'_NucleusBinary_flat.tif'));
% if ~isempty(dir(strcat(channel(1:2),'*'))) && frames == 1
% imwrite(uint8(Ie),strcat(channel,'.tif'),'tiff','WriteMode','overwrite');
% else
% imwrite(uint8(Ie),strcat(channel,'.tif'),'tiff','WriteMode','append');
% end
% cd ..

end

function displaynormalizationfactor(b,c,frames,scenename)
disp(strcat(num2str(b),'_',num2str(c),'_',num2str(frames),'_',scenename))
end

function bw = gaussianBlurz(im,sigma,kernelgsize,varargin)

filtersize = [kernelgsize kernelgsize];
kernelg = fspecial('gaussian',filtersize,sigma);

%% image filtering
gFrame = imfilter(im,kernelg,'repl');

if ~isempty(varargin)
    bw=gFrame.*uint16(varargin{1}>0);
else
    bw=gFrame;
end
end


function bw = logMasked(im,ksize,varargin)
% Discrete Laplacian
kernel = chooseKernel(ksize);
%% image filtering
lapFrame = imfilter(im,kernel,'repl');
if ~isempty(varargin)
    bw=lapFrame.*uint16(varargin{1}>0);
else
    bw=lapFrame;
end
end



function kernel = chooseKernel(ksize)
if ksize ==5
kernel = [-4 -1  0 -1 -4;...
     -1  2  3  2 -1;...
      0  3  4  3  0;...
     -1  2  3  2 -1;...
     -4 -1  0 -1 -4];


% % % -4 -1  0 -1 -4
% % % -1  2  3  2 -1
% % % 0  3  4  3  0
% % % -1  2  3  2 -1
% % % -4 -1  0 -1 -4

elseif ksize == 7
kernel =[-10 -5 -2 -1 -2 -5 -10;... 
    -5  0  3  4  3  0  -5;... 
    -2  3  6  7  6  3  -2;... 
    -1  4  7  8  7  4  -1;... 
    -2  3  6  7  6  3  -2;... 
    -5  0  3  4  3  0  -5;... 
    -10 -5 -2 -1 -2 -5 -10];... 
    
% % % -10 -5 -2 -1 -2 -5 -10 
% % % -5  0  3  4  3  0  -5 
% % % -2  3  6  7  6  3  -2 
% % % -1  4  7  8  7  4  -1 
% % % -2  3  6  7  6  3  -2 
% % % -5  0  3  4  3  0  -5
% % % -10 -5 -2 -1 -2 -5 -10
end
end



function  LoGstack = LaplacianOfGaussianStack(imgstack,dims,ksize)
        LoGstack = zeros(dims(1),dims(2),dims(3));
        for i = 1:size(imgstack,3)
        LoGstack(:,:,i) = logMasked(imgstack(:,:,i),ksize);
        end
        
end


function tsn = determineTimeFrame(foldername)
cd(foldername)
fflist = dir(strcat('*.tif'));
ffnames = {fflist.name};
[a,b,c,d] = regexp(ffnames,'_t[0-9]+');
[a,b,c,d] = regexp(ffnames,'_t[0-9]++');
tnames = cellfun(@(x) x{1},d,'UniformOutput',0);
[a,b,c,d] = regexp(tnames,'t[0-9]++');%added extra step because sometimes the time frame parsin was mistaken
tnames = cellfun(@(x) x{1},d,'UniformOutput',0);
tsn = sort(tnames);
cd ..
end



function PMthreshslider(source,callbackdata)
    global threshinput Imagez
    str = source.String;
     threshinput.(str) =source.Value;
     
    if strcmp(source.String,'channel1')
    elseif strcmp(source.String,'channel2')
        Imagez.NucSeg = segmentNucleus(Imagez.NucImage); %semgnet the image
    elseif strcmp(source.String,'channel3')
        Imagez.CellSeg = segmentCell(Imagez.CellImage);
    end
updatePMseg
updateImage
end



function keypress(fig_obj,~)
global  ImageDetails displaycomments
key = get(fig_obj,'CurrentKey');

switch key
    case '1'
        ImageDetails.Channel = 'EGFP';
        setSceneAndTime
    case '2'
        ImageDetails.Channel = '_Hoechst';
        setSceneAndTime    
    case '3'
        ImageDetails.Channel = 'mKate';
        setSceneAndTime
    case '4'
        ImageDetails.Channel = 'DIC';
        setSceneAndTime
    case '5'
        ImageDetails.Channel = 'BKGbinary';
        setSceneAndTime
    case '6'
        ImageDetails.Channel = 'overlay';
        setSceneAndTime
    case 'q'
        prevscenebutton_Callback([],[])
    case 'w'
        nextscenebutton_Callback([],[])
    case 'a'
        prevbutton_callback([],[])
    case 'f'
        nextbutton_callback([],[])
    case 'd'
        deletebutton_Callback([],[]);
    case 't'
        trackbutton_Callback([],[]);
    case 'e'
        eliminatebutton_Callback([],[]);
    case 'v'
        addareabutton_Callback([],[]);
    case 'r'
        linkCells_Callback([],[]);  
    case 'm'
        displayTrackingButton_Callback([],[])
    case 'g'
        finalbutton_callback([],[])
    case 'z'
        firstbutton_callback([],[])
    case 's'
        saveTrackingFileAs_callback([],[])
    case 'l'
        loadTrackingFile_callback([],[])
    case 'p'
        Plot_callback([],[])
    case 'o'
        labelCells;
    case 'u'
%         if displaycomments==1
%             displaycomments=0;
%         else
            displaycomments=1;
            xy = getxy([],[]);
%         end
    case 'c'
        contrast_Callback([],[])
    case 'k'
        comment_Callback([],[])
    case 'j'
        comment_CallbackJ([],[])
    case 'n'
        PlotCFPnorm_callback([],[])
    case 'b'
        PlotCFPnotnorm_callback([],[])
    case '0'
        displaycomments=1;
        xy = getxy([],[]);
        [~,comments,commentpos,cellidx]=updatecomments(xy);
        setcommentsTracking(comments,commentpos)
        dispxy(xy)
end

end

%choose frames
function nextbutton_callback(~,~)
global framesForDir ImageDetails 

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = framesForDir{1};
end

Idx = strcmp(ImageDetails.Frame,framesForDir);
idx = find(Idx == 1);
if idx == length(framesForDir)
else
idx = idx + 1;
end
ImageDetails.Frame = framesForDir{idx};


setSceneAndTime
end
function prevbutton_callback(~,~) 
global framesForDir ImageDetails 

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = framesForDir{1};
end

Idx = strcmp(ImageDetails.Frame,framesForDir);
idx = find(Idx == 1);
if idx == 1
else
idx = idx - 1;
end
ImageDetails.Frame = framesForDir{idx};


setSceneAndTime
end
function finalbutton_callback(~,~)
global framesForDir ImageDetails 

idx = length(framesForDir);
ImageDetails.Frame = framesForDir{idx};

setSceneAndTime
end
function firstbutton_callback(~,~)
global framesForDir ImageDetails 
idx = 1;
ImageDetails.Frame = framesForDir{idx};

setSceneAndTime
end
function gotobutton_callback(~,~)
global framesForDir ImageDetails 

if isempty(ImageDetails.Frame)
    ImageDetails.Frame = framesForDir{1};
end

prompt = {'Go to which frame'};
dlg_title = 'Go to frame...';
idx = str2num(cell2mat(inputdlg(prompt,dlg_title)));

ImageDetails.Frame = framesForDir{idx};


setSceneAndTime
end


%choose scenes
function nextscenebutton_Callback(~,~) 
global   ImageDetails SceneList A SceneDirectoryPath




if isempty(ImageDetails.Scene)
    ImageDetails.Scene = SceneList{1};
end

Idx = strcmp(ImageDetails.Scene,SceneList);
idx = find(Idx == 1);
if idx == length(SceneList)
else
idx = idx + 1;
end
ImageDetails.Scene = SceneList{idx};



cd(A)
cd('flatfield_corrected')
SceneDirectory = dir (strcat('*',ImageDetails.Scene,'*'));
cd(SceneDirectory.name)
SceneDirectoryPath = pwd;
setSceneAndTime
end
function prevscenebutton_Callback(~,~) 
global   ImageDetails SceneList A SceneDirectoryPath

if isempty(ImageDetails.Scene)
    ImageDetails.Scene = SceneList{1};
end

Idx = strcmp(ImageDetails.Scene,SceneList);
idx = find(Idx == 1);
if idx ==1
else
idx = idx - 1;
end
ImageDetails.Scene = SceneList{idx};

cd(A)
cd('flatfield_corrected')
SceneDirectory = dir (strcat('*',ImageDetails.Scene,'*'));
cd(SceneDirectory.name)
SceneDirectoryPath = pwd;
loadTrackingFile_callback([],[])
setSceneAndTime
end


function popup_menu_Callback_channels(source,~)
global ImageDetails Tracked 

Tracked=[];

% Determine the selected data set.
 str = source.String;
 val = source.Value;
 channel = char(str{val});

ImageDetails.Channel = channel;
setSceneAndTime
end

function popup_menu_Callback(source,~) 
global ImageDetails Tracked 

Tracked=[];

% Determine the selected data set.
 str = source.String;
 val = source.Value;
 pvalue = char(str{val});

ImageDetails.Scene = pvalue;
setSceneAndTime

end

%removeArea
function removeArea_Callback(~,~)
global  ImageDetails framesForDir Tracked imgsize
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add

       If = 1;
       button=1;
       while button==1
      [polyx,polyy,button] = ginput();
      button = round(mean(button));
      
      if button ==1
          M = zeros(1,length(polyx)*2);
          M(1:2:end) = polyx;
          M(2:2:end) = polyy;
          zeroImage = zeros(imgsize);
          zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
          zerogray = rgb2gray(zeroImage);

          if isempty(Tracked{1}.Cellz)


          else  %if there exists segmenttracking already...then load that. 

          imagio = zeros(imgsize);
          imagio(zerogray>0)=1;
          cc = bwconncomp(imagio);
          px = cc.PixelIdxList;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   determine the frame to load
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t = strcmp(framesForDir,ImageDetails.Frame);
            t = find(t==1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    CC = Tracked{t}.Cellz;
    PX = CC.PixelIdxList;        

    idxs = cellfun(@(x) sum(ismember(x,px{1})),PX,'UniformOutput',1);
    index = find(idxs>1);
        if ~isempty(index)
            for abc = index
            oldMass = PX{abc};
            overlap = ismember(oldMass,px{1});
            oldMass(overlap)=[];
            imagio = zeros(imgsize);
            imagio(oldMass)=1;
            cc = bwconncomp(imagio);
            numcells = cc.NumObjects;
                if numcells>1
                    splitcells = cc.PixelIdxList;
                    PX(abc) = {NaN};
                        for nums = 1:numcells
                            PX{end+1} = splitcells{nums};
                        end
                else
                    PX{abc} = oldMass;
                end
            end
%             PX{min(index)} = unique(vertcat(oldMass,px{1}));
%                 if length(index)>1
%                 index(find(index == min(index)))=[];
%                 PX(index) = {NaN};
%                 end
            CC.PixelIdxList = PX;
%             end
        end
    CC.NumObjects = length(CC.PixelIdxList);
        S = regionprops(CC,'Centroid');
        Smat = vertcat(S.Centroid);
        CC.Centroid = Smat;
    Tracked{t}.Cellz = CC;
          end



        nextbutton_callback([],[]);
      end 
      end

end

%delete cells
function deletebutton_Callback(~,~) 
    global imgsize ImageDetails framesForDir Tracked
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
     [cellxx,cellyy] = ginput();
       
% for j = 1:length(cellxx)
%     cellx = cellxx(j);
%     celly = cellyy(j);
%     
    cellx = round(cellxx);
    celly = round(cellyy);
    
    cellind = sub2ind(imgsize,celly,cellx);
    
    for j = 1:length(cellxx)
    if j==1
        idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    else
        idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    end
    end
    
PX(~idxs) = {NaN};
CC.PixelIdxList = PX;
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{t}.Cellz = CC;
    
setSceneAndTime;
   
end
function eliminatebutton_Callback(~,~)

global ImageDetails framesForDir Tracked imgsize

button=1;
while button == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
     [cellxx,cellyy,button] = ginput(1);

     

       
% for j = 1:length(cellxx)
%     cellx = cellxx(j);
%     celly = cellyy(j);
%     
    cellx = round(cellxx);
    celly = round(cellyy);
    
    cellind = sub2ind(imgsize,celly,cellx);
    
    for j = 1:length(cellxx)
    if j==1
        idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    else
        idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
    end
    end
    
PX(~idxs) = {NaN};
CC.PixelIdxList = PX;
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{t}.Cellz = CC;

if button==1
    nextbutton_callback([],[])
elseif button == 3
    prevbutton_callback([],[])
else 
   setSceneAndTime
end

end
end
function destroybutton_Callback(~,~)
%delete a cell from all frames
global ImageDetails framesForDir Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,idxs,[],[]);      
Tracked = Trackedz;


   setSceneAndTime


end
%choose the cells you want
function chosenOnes_Callback(~,~)
%choose the cells you want
global ImageDetails framesForDir Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
% Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 
Tracked = Trackedz;


   setSceneAndTime


end
function chosenOnesEnd_Callback(~,~)
%choose the cells you want
global ImageDetails framesForDir Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
% Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 
Tracked = Trackedz;


   setSceneAndTime

end

function erodeOnes_Callback(~,~)
%choose the cells you want
global ImageDetails framesForDir Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
% Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
% Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 


initialframe = 1;
finalframe = length(Tracked);
Stacked = Tracked;



if isempty(initialframe)
    initialframe=1;
end

if isempty(finalframe)
    finalframe=length(Stacked);
end

se = strel('disk',1);
idxf = find((~idxs) ==1);
for i=initialframe:finalframe 
CC = Stacked{i}.Cellz;
PX = CC.PixelIdxList;   
    for jim = idxf
    pixxies = PX{jim};
        if ~isnan(pixxies)
        imdub = zeros(imgsize);
        imdub(pixxies) = 1;
        imdub = imerode(imdub,se);
        pxtwo = find(imdub==1);
            if ~isempty(pxtwo)
            PX{jim} = pxtwo;
            else
            PX{jim} = NaN;
            end
        end
    end
CC.PixelIdxList = PX;
CC.NumObjects = length(PX);
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Stacked{i}.Cellz = CC;

end
Trackedz=Stacked;








Tracked = Trackedz;


   setSceneAndTime


end



function dilateOnes_Callback(~,~)
%choose the cells you want
global ImageDetails framesForDir Tracked imgsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    
    
% Display mesh plot of the currently selected data.
    [cellxx,cellyy,~] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    
      cellind = sub2ind(imgsize,celly,cellx);
      
      for j=1:length(cellxx)
      if j==1
      idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      else
      idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
      end
      end
      
% Trackedz = crushThem(Tracked,~idxs,length(Tracked),length(Tracked));
% Trackedz = crushThem(Tracked,~idxs,1,length(Tracked)); 


initialframe = 1;
finalframe = length(Tracked);
Stacked = Tracked;



if isempty(initialframe)
    initialframe=1;
end

if isempty(finalframe)
    finalframe=length(Stacked);
end

se = strel('disk',1);
idxf = find((~idxs) ==1);
for i=initialframe:finalframe 
CC = Stacked{i}.Cellz;
PX = CC.PixelIdxList;   
    for jim = idxf
    pixxies = PX{jim};
    disp(pixxies)
        if ~isnan(pixxies)
        imdub = zeros(imgsize);
        imdub(pixxies) = 1;
        imdub = imdilate(imdub,se);
        pxtwo = find(imdub==1);
            if ~isempty(pxtwo)
            PX{jim} = pxtwo;
            else
            PX{jim} = NaN;
            end
        end
    end
CC.PixelIdxList = PX;
CC.NumObjects = length(PX);
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Stacked{i}.Cellz = CC;

end
Trackedz=Stacked;








Tracked = Trackedz;


   setSceneAndTime


end

%internal function for chosenOnes and destroy
function Trackedz = crushThem(Stacked,idxs,initialframe,finalframe)  

if isempty(initialframe)
    initialframe=1;
end

if isempty(finalframe)
    finalframe=length(Stacked);
end

for i=initialframe:finalframe 
CC = Stacked{i}.Cellz;
PX = CC.PixelIdxList;   
PX(~idxs) = {NaN};
CC.PixelIdxList = PX;
CC.NumObjects = length(PX);
Stacked{i}.Cellz = CC;
end
Trackedz=Stacked;
end




%add cells and link cells
function addareabutton_Callback(~,~) 
 global  ImageDetails framesForDir Tracked imgsize
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add

       If = 1;
       button=1;
       while button==1
      [polyx,polyy,button] = ginput();
      button = round(mean(button));
      
      if button ==1
      M = zeros(1,length(polyx)*2);
      M(1:2:end) = polyx;
      M(2:2:end) = polyy;
      zeroImage = zeros(imgsize);
      zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
      zerogray = rgb2gray(zeroImage);

      if isempty(Tracked{1}.Cellz)

          
      else  %if there exists segmenttracking already...then load that. 
        
      imagio = zeros(imgsize);
      imagio(zerogray>0)=1;
      cc = bwconncomp(imagio);
      px = cc.PixelIdxList;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   determine the frame to load
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t = strcmp(framesForDir,ImageDetails.Frame);
        t = find(t==1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;        
        
idxs = cellfun(@(x) sum(ismember(x,px{1})),PX,'UniformOutput',1);
index = find(idxs>1);
    if ~isempty(index)
    newMass = vertcat(PX{index});
    PX{min(index)} = unique(vertcat(newMass,px{1}));
        if length(index)>1
        index(find(index == min(index)))=[];
        PX(index) = {NaN};
        end
    CC.PixelIdxList = PX;
    else
    CC.PixelIdxList = horzcat(PX,px);    
    end
CC.NumObjects = length(CC.PixelIdxList);
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{t}.Cellz = CC;
% 
% makeIMG = cellfun(@(x) length(x)==1,PX,'UniformOutput',1); %choose only the cells without NAN
% newpxlist = horzcat(PX(~makeIMG),px);
% CC.PixelIdxList = newpxlist;
% CC.NumObjects = length(newpxlist);
% 
% segmentimgL = labelmatrix(CC);
% axes(MainAxes)
% imagesc(segmentimgL);
% segmentimg = zeros(size(segmentimgL));
% segmentimg(segmentimgL>0)=1;
% If = segmentimg;
    
      end
      
      
      
    nextbutton_callback([],[]);
      end 
      end
end
function linkCells_Callback(~,~)
global ImageDetails framesForDir Tracked TC imgsize

%this is a quick workaround to get linking working with tracking
%trajectories on
% TC =1;
% setSceneAndTime

button=1;
i=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while button==1
TC=0;
    [xx,yy,button] = ginput(1); %record each click
    
    if button ==1
    t = strcmp(framesForDir,ImageDetails.Frame);
    t = find(t==1);
    cellxx(i)  =    xx; 
    cellyy(i)  =    yy;
    timingkeeper(i) = t;
        if t==length(framesForDir)
            i;
        else
            i=i+1;
        end
    nextbutton_callback([],[]);
    end
    
end
TC=1;

cellx = round(cellxx);
celly = round(cellyy);
cellind = sub2ind(imgsize,celly,cellx);

for i = 1:length(timingkeeper)
t = timingkeeper(i);
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;  
idxlog = cellfun(@(x) isempty(find(x==cellind(i),1)),PX,'UniformOutput',1);
idx(i) = find(idxlog==0);
end

cellindx = 1:length(framesForDir);
timeindx = 1:length(framesForDir);
cellindx(1:timingkeeper(1))=idx(1);
cellindx(timingkeeper(end):end) = idx(end);

for i=timingkeeper
cellindx(i) = idx(i-(min(timingkeeper)-1));  
end


for i=1:length(framesForDir)

CC = Tracked{i}.Cellz;
PX = CC.PixelIdxList;
PX{cellindx(1)} = PX{cellindx(i)};
if cellindx(1)~=cellindx(i)
PX{cellindx(i)} = NaN;
end
CC.PixelIdxList = PX;
    S = regionprops(CC,'Centroid');
    Smat = vertcat(S.Centroid);
    CC.Centroid = Smat;
Tracked{i}.Cellz = CC;
end

setSceneAndTime

end

%plot your cells!
function psettings = PlotSettings_callback(~,~)
prompt = {'tgfbeta frame','last frame'};
dlg_title = 'frames where cells must be tracked...';
framesThatMustBeTracked = str2num(cell2mat(inputdlg(prompt,dlg_title)));
psettings.framesThatMustBeTracked = framesThatMustBeTracked;
end


function PlotCFPnotnorm_callback(~,~)
global toggleCFPnorm
% tcontrast =99;
% lcontrast =1;
toggleCFPnorm = 0;
Plot_callback([],[]);
end
function PlotCFPnorm_callback(~,~)
global toggleCFPnorm
% tcontrast =99;
% lcontrast =1;
toggleCFPnorm = 1;
Plot_callback([],[]);
end

function Plot_callback(~,~)
global toggleCFPnorm SecondPlotAxes Tracked ImageDetails A SceneDirectoryPath timeFrames framesForDir PlotAxes imgsize plottingON psettings cmaplz displaytracking cmap

if plottingON == 0
psettings = PlotSettings_callback([],[]);
plottingON=1;
end
framesThatMustBeTracked = psettings.framesThatMustBeTracked;

for jy = 1:length(framesThatMustBeTracked)
PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
%     makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
    makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
end

[comments,commentpos,cellidx,plotidx]=commentsforplot(Tracked);
makeIMG = makeIMG(1,:)& makeIMG(2,:); 
makeIMGidx = find(makeIMG==1);
if ~isempty(plotidx)
iidd = find(~ismember(makeIMGidx,plotidx));
else
    iidd=[];
end
% makeIMG(makeIMGidx(iidd))=0;
% makeIMGidx = find(makeIMG==1);


smooththat=0;
[Smad,Cfp,mkate,CfpFC,SmadFC,mkateFC,Smadbkg,Cfpbkg,mkatebkg] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirectoryPath,timeFrames,framesForDir,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
smooththat=toggleCFPnorm;

% if smooththat==1
% toplot = SmadFC./CfpFC;
% else
% toplot = SmadFC;    
% end

if smooththat==1
toplot = SmadFC;
else
toplot = SmadFC;    
end

        xmin = 1;
        if xmin <1
            xmin = 1
        end
        
h = plot(SecondPlotAxes,toplot');
    if displaytracking ==1
        for i=1:length(h)
            h(i).Color = cmaplz(i,:);
%             h(i).Color = [0 1 1];
        end
    colormap(cmap);    
    end
        for po = iidd  %fade out the commented cells
        h(po).LineStyle = ':';
        end

SecondPlotAxes.XLim = ([xmin size(toplot,2)]);
SecondPlotAxes.YLim = ([0 6]);

% if smooththat==1
% toplot = Smad./CfpFC;
% else
% toplot = Smad;    
% end

if smooththat==1
toplot = Smad;
else
toplot = Smad;    
end

h = plot(PlotAxes,toplot');
    if displaytracking ==1
        for i=1:length(h)
            h(i).Color = cmaplz(i,:);
        end
    colormap(cmap);    
    end
        for po = iidd  %fade out the commented cells
        h(po).LineStyle = ':';
        end
        

PlotAxes.XLim = ([xmin size(toplot,2)]);
PlotAxes.YLim = ([0 max(max(toplot)).*1.2]);
% PlotAxes.YLim = ([0 1.4e06]);
% PlotAxes.YScale = 'log';


end
function Plot_SpecificCell_callback(~,~)
global plottingTotalOrMedian ThirdPlotAxes toggleCFPnorm Tracked ImageDetails A SceneDirectoryPath timeFrames framesForDir PlotAxes imgsize plottingON psettings

if plottingON == 0
psettings = PlotSettings_callback([],[]);
plottingON=1;
end




framesThatMustBeTracked = psettings.framesThatMustBeTracked;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
      
% Display mesh plot of the currently selected data.
    [cellxx,cellyy] = ginput();
    cellx = round(cellxx);
    celly = round(cellyy);
    cellind = sub2ind(imgsize,celly,cellx);
        for j = 1:length(cellxx)
            if j==1
                idxs = cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
            else
                idxs = idxs & cellfun(@(x) isempty(find(x==cellind(j),1)),PX,'UniformOutput',1);
            end
        end
makeIMG = ~idxs;        
makeIMGidx = find(makeIMG==1);


smooththat=0;
[Smad,Cfp,mkate,CfpFC,SmadFC,mkateFC,Smadbkg,Cfpbkg,mkatebkg] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirectoryPath,timeFrames,framesForDir,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
smooththat=toggleCFPnorm;
if smooththat==1
% toplot = SmadFC./CfpFC;
toplot = SmadFC;
else
toplot = SmadFC;    
end
h = plot(ThirdPlotAxes,toplot','LineWidth',3);
ThirdPlotAxes.XLim = ([1 size(toplot,2)]);
ThirdPlotAxes.YLim = ([0 6]);

end

function [Smad,Cfp,mkate,CfpFC,SmadFC,mkateFC,Smadbkg,Cfpbkg,mkatebkg] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirectoryPath,timeFrames,framesForDir,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat)
global plottingTotalOrMedian


plotbc=1;
plotTracesCell = cell(length(makeIMGidx),length(Tracked));
    for i = 1:length(Tracked)
        PXX = Tracked{i}.Cellz.PixelIdxList;
        plotTracesCell(:,i) = PXX(makeIMG);
    end
    
cd(SceneDirectoryPath)

    if plotbc == 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%   open the image files   %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %  open smad img  %
            smadimgstack = zeros(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','EGFP','*'));
            cd(ChannelDirectory.name)
            for k=1:length(timeFrames)
            imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
            smadimgstack(:,:,k) = double(imread(char(imgfile.name)));
            end
            cd .. 
            
                     %  open mKate img  %
            mkateimgstack = zeros(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','_mKate','*'));
            cd(ChannelDirectory.name)
            for k=1:length(timeFrames)
            imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
            mkateimgstack(:,:,k) = double(imread(char(imgfile.name)));
            end
            cd .. 


                   %    open cfp img  %
            cfpimgstack = zeros(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','CFP','*'));
            cd(ChannelDirectory.name)
            for k=1:length(timeFrames)
            imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
            cfpimgstack(:,:,k) = double(imread(char(imgfile.name)));
            end
            cd ..
                    % open background Logical img  %
            bkglogimgstack = false(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','BKGbinary','*'));
            cd(ChannelDirectory.name)
            for k=1:length(timeFrames)
            imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
            bkgimg = imread(char(imgfile.name));
            bkgimg(bkgimg>0) = 1;
            bkglogimgstack(:,:,k) = ~logical(bkgimg);
            end
            cd ..
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%   open the image files   %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        open smad img  %
            cd(SceneDirectoryPath)
            cd('tiffs')
%             imgfile = dir('EGFP_flat_bleach_corr*');
            imgfile = dir('*EGFP_flat_bleach_corr*');
%             imgfile = dir('poop*');
            if isempty(imgfile)
            smadimgstack = [];
            else
            smadimgstack = double(loadUpTiffStack({imgfile.name}));
            end
                
            
            %                        open mkate img  %
%             cd(SceneDirectoryPath)
%             cd('tiffs')
%             imgfile = dir('_mKate_flat_bleach_corr*');
%             mkateimgstack = double(loadUpTiffStack({imgfile.name}));
            mkateimgstack = zeros(size(smadimgstack));

%                       open cfp img  %
            cd(SceneDirectoryPath)
            cd('tiffs') 
            imgfile = dir('_Hoechst_flat_bleach_corr*');
            if isempty(imgfile)
            cfpimgstack = [];
            else
            cfpimgstack = double(loadUpTiffStack({imgfile.name}));
            end
            cd ..

            if isempty(smadimgstack)
                smadimgstack = ones(size(mkateimgstack));
            end
            if isempty(cfpimgstack)
                cfpimgstack = ones(size(mkateimgstack));
            end
            
%                    open background Logical img  %
            cd(SceneDirectoryPath)        
            bkglogimgstack = false(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','BKGbinary','*'));
            cd(ChannelDirectory.name)
                for k=1:length(timeFrames)
                imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
                bkgimg = imread(char(imgfile.name));
                bkgimg(bkgimg>0) = 1;
                bkglogimgstack(:,:,k) = ~logical(bkgimg);
                end
            cd ..
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end


%perform bkg subtraction
Smadbkg = zeros(1,length(timeFrames));
Cfpbkg = zeros(1,length(timeFrames));
mkatebkg = zeros(1,length(timeFrames));
bkgstd = zeros(1,length(timeFrames));
for k=1:length(timeFrames)
    bkglog = bkglogimgstack(:,:,k);
    smadimg = double(smadimgstack(:,:,k));
    cfpimg = double(cfpimgstack(:,:,k));
    mkateimg = double(mkateimgstack(:,:,k));
    %background subtraction is just subtraction with a value
    Smadbkg(k) = nanmedian(smadimg(bkglog));
    Cfpbkg(k) = nanmedian(cfpimg(bkglog));
    mkatebkg(k) = nanmedian(mkateimg(bkglog));
    bkgstd(k) = nanstd(mkateimg(bkglog));
    smadimgstack(:,:,k) = smadimg-Smadbkg(k);
    cfpimgstack(:,:,k) = cfpimg-Cfpbkg(k);
    mkateimgstack(:,:,k) = mkateimg - mkatebkg(k);
    %background subtraction is subtraction with an interpolated image
%     smadbkgimg = regionfill(smadimg,~bkglog);
%     cfpbkgimg = regionfill(cfpimg,~bkglog); %fill in the regions where bkglog is 0
%     smadimgstack(:,:,k) = smadimgstack(:,:,k)-smadbkgimg;
%     cfpimgstack(:,:,k) = cfpimgstack(:,:,k)-cfpbkgimg;
end


%extract pixel intensities
smadpxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));
cfppxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));

for i = 1:size(plotTracesCell,2)
    smadimg = double(squeeze(smadimgstack(:,:,i)));
    cfpimg = double(squeeze(cfpimgstack(:,:,i)));
    mkateimg = double(squeeze(mkateimgstack(:,:,i)));
    for j=1:size(plotTracesCell,1)
    pxidx = plotTracesCell{j,i};
        if ~isnan(pxidx)
        smadpxls(j,i) = {smadimg(pxidx)};
        cfppxls(j,i) = {cfpimg(pxidx)};
        mkatepxls(j,i) = {mkateimg(pxidx)};
        else
        smadpxls(j,i) = {123456789};
        cfppxls(j,i) = {123456789};
        mkatepxls(j,i) = {123456789};
        end
    end
end

if strcmpi(plottingTotalOrMedian,'total')
    Smad = cellfun(@nansum,smadpxls,'UniformOutput',1);
    % Smad = cellfun(@nanmean,smadpxls,'UniformOutput',1);
    Smad(Smad==123456789) = NaN;
    Cfp = cellfun(@nansum,cfppxls,'UniformOutput',1);
    Cfp(Cfp==123456789) = NaN;
    mkate = cellfun(@nansum,mkatepxls,'UniformOutput',1);
    mkate(mkate==123456789) = NaN;
elseif strcmpi(plottingTotalOrMedian,'median')
    Smad = cellfun(@nanmedian,smadpxls,'UniformOutput',1);
    % Smad = cellfun(@nanmean,smadpxls,'UniformOutput',1);
    Smad(Smad==123456789) = NaN;
    Cfp = cellfun(@nanmedian,cfppxls,'UniformOutput',1);
    Cfp(Cfp==123456789) = NaN;
    mkate = cellfun(@nanmedian,mkatepxls,'UniformOutput',1);
    mkate(mkate==123456789) = NaN;
end

basalSUB = framesThatMustBeTracked(1)-7;
if basalSUB<1
    basalSUB = framesThatMustBeTracked-1;
end
basalcfp = nanmean(Cfp(:,framesThatMustBeTracked(1)-basalSUB:framesThatMustBeTracked(1)),2);

CfpFC = zeros(size(Cfp));
    for i = 1:size(Cfp,2)
       CfpFC(:,i) = Cfp(:,i)./basalcfp; 
    end
    
    if smooththat==1;
    %%%%%%%%%%%%%%%%%%%
    Smad = Smad./CfpFC;
    %%%%%%%%%%%%%%%%%%%
    end
basalsmad = nanmean(Smad(:,framesThatMustBeTracked(1)-basalSUB:framesThatMustBeTracked(1)),2);
SmadFC = zeros(size(Smad));
    for i = 1:size(Smad,2)
       SmadFC(:,i) = Smad(:,i)./basalsmad; 
    end
    
    basalmkate = nanmean(mkate(:,framesThatMustBeTracked(1)-basalSUB:framesThatMustBeTracked(1)),2);
mkateFC = zeros(size(mkate));
    for i = 1:size(mkate,2)
       mkateFC(:,i) = mkate(:,i)./basalmkate; 
    end
    
    
end
function plotStruct = plotthemfunctionToStructure(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirectoryPath,timeFrames,framesForDir,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat)
plotStruct = struct();

plotTracesCell = cell(length(makeIMGidx),length(Tracked));
    for i = 1:length(Tracked)
        PXX = Tracked{i}.Cellz.PixelIdxList;
        plotTracesCell(:,i) = PXX(makeIMG);
    end
    
cd(SceneDirectoryPath)

plotbc=1; %0 or 1 to chose bleach corr or not bleach corr
    if plotbc == 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%   open the image files   %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %  open smad img  %
            smadimgstack = zeros(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','EGFP','*'));
            cd(ChannelDirectory.name)
            for k=1:length(timeFrames)
            imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
            smadimgstack(:,:,k) = double(imread(char(imgfile.name)));
            end
            cd .. 
            
                     %  open mKate img  %
            mkateimgstack = zeros(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','_mKate','*'));
            cd(ChannelDirectory.name)
            for k=1:length(timeFrames)
            imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
            mkateimgstack(:,:,k) = double(imread(char(imgfile.name)));
            end
            cd .. 

                   %    open cfp img  %
            cfpimgstack = zeros(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','CFP','*'));
            cd(ChannelDirectory.name)
            for k=1:length(timeFrames)
            imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
            cfpimgstack(:,:,k) = double(imread(char(imgfile.name)));
            end
            cd ..
                    % open background Logical img  %
            bkglogimgstack = false(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','BKGbinary','*'));
            cd(ChannelDirectory.name)
            for k=1:length(timeFrames)
            imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
            bkgimg = imread(char(imgfile.name));
            bkgimg(bkgimg>0) = 1;
            bkglogimgstack(:,:,k) = ~logical(bkgimg);
            end
            cd ..
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    else 
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%   open the image files   %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        open smad img  %
            cd(SceneDirectoryPath)
            cd('tiffs')
            imgfile = dir('*EGFP_flat_bleach_corr*');
            if isempty(imgfile)
                smadimgstack = [];
            else
               smadimgstack = double(loadUpTiffStack({imgfile.name}));
            end
                
  

%                       open cfp img  %
            cd(SceneDirectoryPath)
            cd('tiffs') 
            imgfile = dir('_mKate_flat_bleach_corr*');
            if isempty(imgfile)
                cfpimgstack = [];
            else
                cfpimgstack = double(loadUpTiffStack({imgfile.name}));
            end
            cd ..

%                       open mKate img  %            
            mkateimgstack = zeros(size(smadimgstack));
            
            
%                    open background Logical img  %
            cd(SceneDirectoryPath)        
            bkglogimgstack = false(imgsize(1),imgsize(2),length(Tracked));
            ChannelDirectory = dir(strcat('*','BKGbinary','*'));
            cd(ChannelDirectory.name)
                for k=1:length(timeFrames)
                imgfile = dir(strcat('*',framesForDir{k},'*.tif'));
                bkgimg = imread(char(imgfile.name));
                bkgimg(bkgimg>0) = 1;
                bkglogimgstack(:,:,k) = ~logical(bkgimg);
                end
            cd ..
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end


%perform bkg subtraction
Smadbkg = zeros(1,length(timeFrames));
Cfpbkg = zeros(1,length(timeFrames));
mkatebkg = zeros(1,length(timeFrames));
bkgstd = zeros(1,length(timeFrames));
for k=1:length(timeFrames)
    bkglog = bkglogimgstack(:,:,k);
    smadimg = double(smadimgstack(:,:,k));
    cfpimg = double(cfpimgstack(:,:,k));
    mkateimg = double(mkateimgstack(:,:,k));
    %background subtraction is just subtraction with a value
        Smadbkg(k) = nanmedian(smadimg(bkglog));
        Cfpbkg(k) = nanmedian(cfpimg(bkglog));
        mkatebkg(k) = nanmedian(mkateimg(bkglog));
        bkgstd(k) = nanstd(mkateimg(bkglog));
        smadimgstack(:,:,k) = smadimg-Smadbkg(k);
        cfpimgstack(:,:,k) = cfpimg-Cfpbkg(k);
        mkateimgstack(:,:,k) = mkateimg - mkatebkg(k);
    %background subtraction is subtraction with an interpolated image
%         smadbkgimg = regionfill(smadimg,~bkglog);
%         cfpbkgimg = regionfill(cfpimg,~bkglog); %fill in the regions where bkglog is 0
%         smadimgstack(:,:,k) = smadimgstack(:,:,k)-smadbkgimg;
%         cfpimgstack(:,:,k) = cfpimgstack(:,:,k)-cfpbkgimg;
end


%extract pixel intensities%%%%%%%%%%%%%%%%
%intialize pxl arrays
smadpxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));
cfppxls = cell(size(plotTracesCell,1),size(plotTracesCell,2));

for i = 1:size(plotTracesCell,2)
    smadimg = double(squeeze(smadimgstack(:,:,i)));
    cfpimg = double(squeeze(cfpimgstack(:,:,i)));
    mkateimg = double(squeeze(mkateimgstack(:,:,i)));
    for j=1:size(plotTracesCell,1)
    pxidx = plotTracesCell{j,i};
        if ~isnan(pxidx)
        smadpxls(j,i) = {smadimg(pxidx)}; %these are arrays of pixel intensity values
        cfppxls(j,i) = {cfpimg(pxidx)};
        mkatepxls(j,i) = {mkateimg(pxidx)};
        else
        smadpxls(j,i) = {123456789};
        cfppxls(j,i) = {123456789};
        mkatepxls(j,i) = {123456789};
        end
    end
end

%determine median pxl intensities
    Smad = cellfun(@nanmedian,smadpxls,'UniformOutput',1);
    Smad(Smad==123456789) = NaN;
    Cfp = cellfun(@nanmedian,cfppxls,'UniformOutput',1);
    Cfp(Cfp==123456789) = NaN;
    mkate = cellfun(@nanmedian,mkatepxls,'UniformOutput',1);
    mkate(mkate==123456789) = NaN;
for i = 1:size(smadpxls,1)
    plotStruct(i).medianNucEGFP = Smad(i,:);
    plotStruct(i).medianNucRFP = Cfp(i,:);
end

%determine total pxl intensities
    Smad = cellfun(@nansum,smadpxls,'UniformOutput',1);
    Smad(Smad==123456789) = NaN;
    Cfp = cellfun(@nansum,cfppxls,'UniformOutput',1);
    Cfp(Cfp==123456789) = NaN;
    mkate = cellfun(@nansum,mkatepxls,'UniformOutput',1);
    mkate(mkate==123456789) = NaN;
for i = 1:size(smadpxls,1)
    plotStruct(i).totalNucEGFP = Smad(i,:);
    plotStruct(i).totalNucRFP = Cfp(i,:);
end

%determine mean pxl intensities
    Smad = cellfun(@nanmean,smadpxls,'UniformOutput',1);
    Smad(Smad==123456789) = NaN;
    Cfp = cellfun(@nanmean,cfppxls,'UniformOutput',1);
    Cfp(Cfp==123456789) = NaN;
    mkate = cellfun(@nanmean,mkatepxls,'UniformOutput',1);
    mkate(mkate==123456789) = NaN;
for i = 1:size(smadpxls,1)
    plotStruct(i).meanNucEGFP = Smad(i,:);
    plotStruct(i).meanNucRFP = Cfp(i,:);
end



end



%%%%comments!!!
function xy = getxy(~,~)
global Tracked plottingON psettings imgsize 
    
    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
    t=length(Tracked);
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
        PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
        makeIMG = false(length(framesThatMustBeTracked),length(PX));
            for jy = 1:length(framesThatMustBeTracked)
                PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
%                 makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
                makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
            end
        makeIMG = makeIMG(1,:)&makeIMG(2,:);
        makeIMGidx = find(makeIMG==1);
        
         PXX = Tracked{t}.Cellz.PixelIdxList;
    PX = PXX(makeIMG);
    mx = nan(1,length(PX));
    my = nan(1,length(PX));
        for j = 1:length(PX)
        px = PX{j};
        y = rem(px-1,imgsize(1))+1; %these two lines replace ind2sub
        x = (px-y)/imgsize(2) + 1;  %these two lines replace ind2sub
        %%%% sort and select middle index is faster than median %%%
        sx = sort(x);   
        sy = sort(y);
        pseudomean = round(length(sx)./2);
            if pseudomean == 0
                mx(j) = NaN; 
                my(j) = NaN;
            else    
                mx(j) = sx(pseudomean);  %use the make Centroids index to keep the centroids the same color when plotting
                my(j) = sy(pseudomean);  
            end
        end
    xy = horzcat(mx',my');
     

end
function dispxy(xy)    
global framesForDir displaycomments Tracked
t = length(Tracked);
for i = 1:size(xy,1)
text(xy(i,1)+1,xy(i,2)+1,num2str(i),'Color',[0 0 0],'FontSize',12);hold on
text(xy(i,1),xy(i,2),num2str(i),'Color',[1 1 0],'FontSize',12);hold on
end

   xxyy=xy;
   xxyy(:,1) = xxyy(:,1)+10;
   xxyy(:,2) = xxyy(:,2)+10;
if t==length(framesForDir)
    if displaycomments==1
        comments = Tracked{t}.comments;
        commentpos = Tracked{t}.commentpos;
%         for i = 1:size(xxyy,1)
%         for i = 1:length(comments)           
            lll = length(comments);
            lllx = size(xxyy,1);
            if lll <lllx
                lcc = lll;
            else
                lcc =lllx;
            end
        for i = 1:lcc
        text(xxyy(i,1)+1,xxyy(i,2)+1,comments{i},'Color',[0 0 0],'FontSize',8);hold on
        text(xxyy(i,1),xxyy(i,2),comments{i},'Color',[1 1 0],'FontSize',8);hold on
        end
    end
end
end


function [comments,commentpos,cellidx,plotidx]=commentsforplot(Tracked)
%choose the cells you want to comment on
global ImageDetails framesForDir imgsize 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    

      
      ix = cellfun(@length,PX,'UniformOutput',1);
      idxs = ix>1;
       
idx = find(idxs==1);


fnames = fieldnames(Tracked{t});
if sum(strcmp(fnames,'comments'));
oldcomments = Tracked{t}.comments;
oldcommentpos = Tracked{t}.commentpos;
comments=[];
commentpos =[];
else
oldcomments = [];
oldcommentpos = [];
comments = [];
commentpos =[];
end

comments = cell(1,length(oldcommentpos));
commentpos = zeros(1,length(oldcommentpos));
cellidx =zeros(1,length(oldcommentpos));

indies = oldcommentpos;

if ~isempty(oldcommentpos)
    for jim = idx
        cycle=0;
        for i = 1:length(oldcommentpos)
            cycle=cycle+1;
       px = PX{jim};
       alreadycommented = ismember(oldcommentpos(i),px);
       indiidx = find(ismember(indies,px)==1);
           if alreadycommented ==1
            comments{indiidx} = oldcomments{i};
            commentpos(indiidx) = oldcommentpos(i);
            cellidx(indiidx) = jim;
           else
               
           end
        end
    end
else


for i=idx
    px = PX{i};
    indiidx = find(ismember(indies,px)==1);
    commentpos(indiidx) = indies(indiidx);  
end

end  


for i=1:length(comments)
    cellnums{i} = strcat('cell#',num2str(i));
    if isempty(comments{i})
        comments{i}='';
    end
end

stringsarray = {'overdriver','no cfp','nocfp','dimmer',...
            'overlap','remove','doublenuc','saturated','big'};
        for ubi = 1:length(stringsarray)
            if ubi==1
            stringer = strcat('(',stringsarray{ubi},'|');
            elseif ubi==length(stringsarray)
                stringer = horzcat(stringer,strcat(stringsarray{ubi},')'));
            else
                stringer = horzcat(stringer,strcat(stringsarray{ubi},'|'));
            end
        end
        
[o,oo,~,~] = regexp(comments,stringer);
overlay = cellfun(@length,o,'UniformOutput',1);
plotidx = cellidx(~overlay);
stophere=1;

end
function setcommentsTracking(comments,commentpos)
global Tracked
t=length(Tracked);
Tracked{t}.comments= comments;
Tracked{t}.commentpos = commentpos;
end
function [cellnums,comments,commentpos,cellidx]=updatecomments(xy)
%choose the cells you want to comment on
global ImageDetails framesForDir Tracked imgsize 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = Tracked{t}.Cellz;
PX = CC.PixelIdxList;   
    

      
      ix = cellfun(@length,PX,'UniformOutput',1);
      idxs = ix>1;
       
idx = find(idxs==1);


fnames = fieldnames(Tracked{t});
if sum(strcmp(fnames,'comments'));
oldcomments = Tracked{t}.comments;
oldcommentpos = Tracked{t}.commentpos;
comments=[];
commentpos =[];
else
oldcomments = [];
oldcommentpos = [];
comments = [];
commentpos =[];
end

comments = cell(1,size(xy,1));
commentpos = zeros(1,size(xy,1));
cellidx =zeros(1,size(xy,1));

indies = sub2ind([imgsize],xy(:,2),xy(:,1));

if ~isempty(oldcommentpos)
    for jim = idx
        cycle=0;
        for i = 1:length(oldcommentpos)
            if i==5
               stophere=1;
            end
            cycle=cycle+1;
       px = PX{jim};
       alreadycommented = ismember(oldcommentpos(i),px);
       indiidx = find(ismember(indies,px)==1);
           if alreadycommented ==1
            comments{indiidx} = oldcomments{i};
            commentpos(indiidx) = oldcommentpos(i);
            cellidx(indiidx) = jim;
           else
               
           end
        end
    end
    
    
    
else



for i=idx
    px = PX{i};
    indiidx = find(ismember(indies,px)==1);
    commentpos(indiidx) = indies(indiidx);  
end

end  

%try this out 2016_02_17
aa = find(commentpos==0);
commentpos(aa) = indies(aa);


for i=1:length(comments)
    cellnums{i} = strcat('cell#',num2str(i));
    if isempty(comments{i})
        comments{i}='';
    end
end

end

function comments = nowaddcomments(xy,cellnums,comments,commentpos)
global Tracked
t=length(Tracked);



prompt = cellnums;
dlg_title = 'Comment for selected cells:';
num_lines = 1;
defaultans = comments;
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ~isempty(answer)
comments = answer;

Tracked{t}.comments= comments;
Tracked{t}.commentpos = commentpos;
end
% displaycomments =1;

end

function comment_Callback(~,~)
global displaycomments Tracked

    fnames = fieldnames(Tracked{length(Tracked)});
    if sum(strcmp(fnames,'comments'))
        displaycomments =1;
    else
        displaycomments =0;
    end
   xy = getxy([],[]);
   dispxy(xy)
   [cellnums,comments,commentpos,cellidx]=updatecomments(xy);
   comments = nowaddcomments(xy,cellnums,comments,commentpos);
   displaycomments =1;
   dispxy(xy)   
%    displaycomments =0;


end

function comments = nowaddcommentsJ(xy,cellnums,comments,commentpos)
global Tracked
t=length(Tracked);



prompt = cellnums;
dlg_title = 'Comment for selected cells:';
num_lines = 1;
defaultans = comments;
% answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
answer = defaultans;

if ~isempty(answer)
comments = answer;

Tracked{t}.comments= comments;
Tracked{t}.commentpos = commentpos;
end
% displaycomments =1;

end

function comment_CallbackJ(~,~)
global ExportNameKey ExportName OGExpDate displaycomments SceneList Tracked ImageDetails A SceneDirectoryPath timeFrames framesForDir PlotAxes imgsize plottingON psettings

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(SceneDirectoryPath)
cd ..
    for scenenumber = 1:length(SceneList)
        cd(SceneDirectoryPath)
        cd ..
        sceneN = char(SceneList{scenenumber});
        disp(sceneN)
        scenedir = dir(strcat('*',sceneN,'*'));
        scenedirname = char({scenedir.name});
        cd(scenedirname)
        SceneDirPath = char({pwd});
        
        
        trackfile = dir(strcat(ExportNameKey,ExportName,'.mat'));
%         trackfile = dir('finalfricktrack.mat');
        trackfilename = char({trackfile.name});
        
            if ~isempty(trackfilename)
                load(trackfilename)
                fnames = fieldnames(Tracked{length(Tracked)});
                    if sum(strcmp(fnames,'comments'))
                        displaycomments =1;
                    else
                        displaycomments =0;
                    end
                xy = getxy([],[]);
                dispxy(xy)
                [cellnums,comments,commentpos,cellidx]=updatecomments(xy);
                comments = nowaddcommentsJ(xy,cellnums,comments,commentpos);
                displaycomments =1;
                dispxy(xy)   
                saveTrackingFileAs_callbackJ([],[],SceneDirPath)
            else
                stopmehere=1;
            end
    end
      
end

function exportCells(~,~)
global ExportNameKey ExportName exportdir OGExpDate SceneList Tracked ImageDetails A SceneDirectoryPath timeFrames framesForDir PlotAxes imgsize plottingON psettings
exportStruct = struct();

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(SceneDirectoryPath)
cd ..
    for scenenumber = 1:length(SceneList)
        cd(SceneDirectoryPath)
        cd ..
        sceneN = char(SceneList{scenenumber});
        disp(sceneN)
        scenedir = dir(strcat('*',sceneN,'*'));
        scenedirname = char({scenedir.name});
        cd(scenedirname)
        SceneDirPath = char({pwd});
        

        trackfile = dir(strcat(ExportNameKey,ExportName,'.mat'));
%         trackfile = dir('finalfricktrack.mat');
        trackfilename = char({trackfile.name});
        
            if ~isempty(trackfilename)
                load(trackfilename)
                PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
                makeIMG = false(length(framesThatMustBeTracked),length(PX));
                   
                    for jy = 1:length(framesThatMustBeTracked)
                    PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
%                     makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
                    makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
                    end
                    
                makeIMG = makeIMG(1,:)&makeIMG(2,:);
                makeIMGidx = find(makeIMG==1);
                smooththat=0;
                [Smad,Cfp,mkate,CfpFC,SmadFC,mkateFC,Smadbkg,Cfpbkg,mkatebkg] = plotthemfunction(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirPath,timeFrames,framesForDir,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
                plotStruct = plotthemfunctionToStructure(framesThatMustBeTracked,Tracked,A,ImageDetails,SceneDirPath,timeFrames,framesForDir,PlotAxes,imgsize,plottingON,psettings,makeIMG,makeIMGidx,smooththat);
                
                fnames = fieldnames(plotStruct);
                if isempty(fieldnames(exportStruct)) %if exportStruct is empty 
                    idx = 0;
                    for i = 1:length(plotStruct)
                        for j = 1:length(fnames)
                            exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});
                            exportStruct(idx+i).scene = sceneN;
                            exportStruct(idx+i).cellID = i;
                        end
                    end
                else    %if fields are defined, append the cell data to the next available index
                    idx = length(exportStruct);
                    for i = 1:length(plotStruct)
                        for j = 1:length(fnames)
                            exportStruct(idx+i).(fnames{j}) = plotStruct(i).(fnames{j});
                            exportStruct(idx+i).scene = sceneN;
                            exportStruct(idx+i).cellID = i;
                        end
                    end
                end
                
                
%                 exportStruct(:).scene = sceneN;
        
stophere=1;
%         xlswrite(filename,eggcell,sceneN);
                end
    end
        cd(exportdir)
        filename = strcat(OGExpDate,'_tracking_export.mat'); 
        save(filename,'exportStruct');
end


function xy = labelCells(~,~)

global ExportNameKey ExportName displaycomments OGExpDate SceneList Tracked ImageDetails A SceneDirectoryPath timeFrames framesForDir PlotAxes imgsize plottingON psettings

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(SceneDirectoryPath)
cd ..

sceneN = char(ImageDetails.Scene);
disp(sceneN)
scenedir = dir(strcat('*',sceneN,'*'));
scenedirname = char({scenedir.name});
cd(scenedirname)
SceneDirPath = char({pwd});

trackfile = dir(strcat(ExportNameKey,ExportName,'.mat'));
% trackfile = dir('finalfricktrack.mat');
trackfilename = char({trackfile.name});

    if ~isempty(trackfilename)
        load(trackfilename)
        
        PX = Tracked{framesThatMustBeTracked(1)}.Cellz.PixelIdxList;
        makeIMG = false(length(framesThatMustBeTracked),length(PX));
            for jy = 1:length(framesThatMustBeTracked)
                PX = Tracked{framesThatMustBeTracked(jy)}.Cellz.PixelIdxList;
%                 makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)==1,PX,'UniformOutput',1)); %choose only the cells without NAN
                makeIMG(jy,:) = ~logical(cellfun(@(x) length(x)<2,PX,'UniformOutput',1)); %choose only the cells without NAN
            end
        makeIMG = makeIMG(1,:)&makeIMG(2,:);
        makeIMGidx = find(makeIMG==1);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PXX = Tracked{t}.Cellz.PixelIdxList;
    PX = PXX(makeIMG);
    mx = nan(1,length(PX));
    my = nan(1,length(PX));
        for j = 1:length(PX)
        px = PX{j};
        y = rem(px-1,imgsize(1))+1; %these two lines replace ind2sub
        x = (px-y)/imgsize(2) + 1;  %these two lines replace ind2sub
        %%%% sort and select middle index is faster than median %%%
        sx = sort(x);   
        sy = sort(y);
        pseudomean = round(length(sx)./2);
            if pseudomean == 0
                mx(j) = NaN; 
                my(j) = NaN;
            else    
                mx(j) = sx(pseudomean);  %use the make Centroids index to keep the centroids the same color when plotting
                my(j) = sy(pseudomean);  
            end
        end
    xy = horzcat(mx',my');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(xy,1)
text(xy(i,1)+1,xy(i,2)+1,num2str(i),'Color',[0 0 0],'FontSize',12);hold on
text(xy(i,1),xy(i,2),num2str(i),'Color',[1 1 1],'FontSize',12);hold on
end

   xxyy=xy;
   xxyy(:,1) = xxyy(:,1)+10;
   xxyy(:,2) = xxyy(:,2)+10;
fnames = fieldnames(Tracked{t});
if sum(strcmp(fnames,'comments'))
    if t==length(framesForDir)
        if displaycomments==1
            comments = Tracked{t}.comments;
            commentpos = Tracked{t}.commentpos;
%             for i = 1:size(xxyy,1)
%             for i = 1:length(comments)
            lll = length(comments);
            lllx = size(xxyy,1);
            if lll <lllx
                lcc = lll;
            else
                lcc =lllx;
            end
            for i = 1:lcc
            text(xxyy(i,1)+1,xxyy(i,2)+1,comments{i},'Color',[0 0 0],'FontSize',8);hold on
            text(xxyy(i,1),xxyy(i,2),comments{i},'Color',[1 1 0],'FontSize',8);hold on
            end
        end
    end     
end

end

function exportFrames(~,~)
global ExportNameKey ExportName OGExpDate imgsize displaytracking SceneList Tracked ImageDetails A SceneDirectoryPath timeFrames framesForDir PlotAxes imgsize plottingON psettings

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
    
%   determine the frame to load
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);


framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(SceneDirectoryPath)
cd ..
CENTROID = struct();
    
        cd(SceneDirectoryPath)
        cd ..
%         sceneN = char(SceneList{scenenumber});
%         disp(sceneN)
%         ImageDetails.Scene = sceneN;
        sceneN = ImageDetails.Scene;
        scenedir = dir(strcat('*',sceneN,'*'));
        scenedirname = char({scenedir.name});
        cd(scenedirname)
        trackfile = dir(strcat(ExportNameKey,ExportName,'.mat'));
%         trackfile = dir('finalfricktrack.mat');
        trackfilename = char({trackfile.name});

        %%% save figure image and centroids %%%
        if ~isempty(trackfilename)
            load(trackfilename)
            setSceneAndTime;
            displaytracking =1;
%             finalbutton_callback([],[]);
            for i=1:length(framesForDir)
                if i==1
                    firstbutton_callback([],[]);
                else
                    nextbutton_callback([],[]);
                end
                
%        xy = labelCells([],[]);
            set(gcf,'Units','Pixels');
            set(gca,'Units','Pixels');
            P = get(gca,'pos');
            F = getframe(gcf,P);
            [X,Map] = frame2im(F);
            set(gcf,'Units','normalized');
            set(gca,'Units','normalized');
            cd(SceneDirectoryPath)
            cd ..
            cd ..
            cd('ANNOTATIONS')
            imwrite(X,strcat('s',sceneN,'-t',num2str(i),'.jpg'),'JPEG'); %save image
            cd .. 
            cd('CENTROIDS')
            cd ..
            
%             CENTROIDS.imgsize = imgsize;
%             CENTROIDS.centroids = xy;
%             CENTROIDS.scene = sceneN;
%             save(strcat('CENTROIDS-',sceneN,'.mat'),'CENTROIDS'); %save CENTROIDS
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        

        
        
        
 
end
function exportLabels(~,~)
global ExportNameKey ExportName OGExpDate imgsize displaytracking SceneList Tracked ImageDetails A SceneDirectoryPath timeFrames framesForDir PlotAxes imgsize plottingON psettings adjuster

    if plottingON == 0
    psettings = PlotSettings_callback([],[]);
    plottingON=1;
    end
    
framesThatMustBeTracked = psettings.framesThatMustBeTracked;
cd(SceneDirectoryPath)
cd ..
CENTROID = struct();
    for scenenumber = 1:length(SceneList)
        adjuster=1;
        cd(SceneDirectoryPath)
        cd ..
        sceneN = char(SceneList{scenenumber});
        disp(sceneN)
        ImageDetails.Scene = sceneN;
        scenedir = dir(strcat('*',sceneN,'*'));
        scenedirname = char({scenedir.name});
        cd(scenedirname)
        
        trackfile = dir(strcat(ExportNameKey,ExportName,'.mat'));
%         trackfile = dir('finalfricktrack.mat');
        trackfilename = char({trackfile.name});

        %%% save figure image and centroids %%%
        if ~isempty(trackfilename)
            load(trackfilename)
            setSceneAndTime;
            displaytracking =1;
            finalbutton_callback([],[]);
       xy = labelCells([],[]);
            set(gcf,'Units','Pixels');
            set(gca,'Units','Pixels');
            P = get(gca,'pos');
            F = getframe(gcf,P);
            [X,Map] = frame2im(F);
            set(gcf,'Units','normalized');
            set(gca,'Units','normalized');
            cd(SceneDirectoryPath)
            cd ..
            cd ..
            cd('ANNOTATIONS')
            imwrite(X,strcat(sceneN,'.jpg'),'JPEG'); %save image
            cd .. 
            cd('CENTROIDS')
            
            CENTROIDS.imgsize = imgsize;
            CENTROIDS.centroids = xy;
            CENTROIDS.scene = sceneN;
            save(strcat('CENTROIDS-',sceneN,'.mat'),'CENTROIDS'); %save CENTROIDS
        end
        %%%%%%%%%%%%%%%%%%%%%%%
        

        
        
        
    end
    adjuster=0;
end




%% Image Display functions
function setSceneAndTime
global TC A  framesForDir ImageDetails  Tracked SceneList  SceneDirectoryPath imgfile imgsize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the channel directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(A)
cd('flatfield_corrected')

    if isempty(ImageDetails.Scene)
        ImageDetails.Scene = SceneList{1};
    end
SceneDirectory = dir (strcat('*',ImageDetails.Scene,'*'));
cd(SceneDirectory.name)
SceneDirectoryPath = pwd;

    if isempty(ImageDetails.Channel)
        ImageDetails.Channel = 'EGFP';
    end
ChannelDirectory = dir(strcat('*',ImageDetails.Channel,'_*'));
    if isempty(ChannelDirectory) && ~strcmp(ImageDetails.Channel,'overlay')
        ImageDetails.Channel = 'mKate';
        ChannelDirectory = dir(strcat('*',ImageDetails.Channel,'_*'));
    elseif isempty(ChannelDirectory)
        ChannelDirectory = dir(strcat('*','mKate','_*'));
    end
cd(ChannelDirectory.name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ImageDetails.Frame)
       ImageDetails.Frame = framesForDir{1};
       t = strcmp(framesForDir,framesForDir{1});
    else
        t = strcmp(framesForDir,ImageDetails.Frame);
    end
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose the channel images
%options are overlay of background
%overlay of fluorescent channels
%normal image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgfile = dir(strcat('*',ImageDetails.Frame,'*.tif'));
    if strcmp(ImageDetails.Channel,'BKGbinary')                    %overlay background
        bkgimg = false(512,512,1);
        bkgimg = imread(char(imgfile.name));
        bkgimg(bkgimg>0) = 1;
        channelimg = ~logical(bkgimg);


%         ChannelDirectory = dir(strcat('*','EGFP_','*'));
%         cd(ChannelDirectory.name)
        cd ..        
        cd('tiffs')
%         ff = dir(strcat('EGFP','*'));
        ff = dir(strcat('*','EGFP','*'));
            if isempty(ff)
                ff = dir(strcat('*','mKate','*'));
            end
        channelimg = double(loadUpTiffStackFrame(char(ff.name),t));
% 
%         cimgfile = dir(strcat('*',ImageDetails.Frame,'*.tif'));
%         cimg = imread(char(cimgfile.name));
%         channelimg = double(cimg);
        prim = imdilate(bwperim(~logical(bkgimg)),strel('square',2));
        channelimg(prim) = max(max(channelimg));
    elseif strcmp(ImageDetails.Channel,'overlay')
           
        
         cd ..        
        cd('tiffs')
        ff = dir(strcat('EGFP','*'));
            if isempty(ff)
               channelimg = zeros(512,512);
            else
            channelimg = double(loadUpTiffStackFrame(char(ff.name),t));
            end
        imgone = channelimg;
        cd ..
        
        cd('tiffs')
        ff = dir(strcat('*','Hoechst','*'));
            if isempty(ff)
               channelimg = zeros(512,512);
            else
            channelimg = double(loadUpTiffStackFrame(char(ff.name),t));
            end
        imgtwo = channelimg;
        cd .. 
        
        
        cd('tiffs')
        ff = dir(strcat('*','mKate_','*'));
        channelimg = double(loadUpTiffStackFrame(char(ff.name),t));
        imgthree = channelimg;
    
        ff = dir(strcat('*','DIC','*'));
        channelimg = double(loadUpTiffStackFrame(char(ff.name),t));
        imgfour = channelimg;
        
        
        channelimg = zeros(size(imgone,1),size(imgone,2),3);
%         imgone(imgone<imgfour) = imgfour(imgone<imgfour);
%         imgtwo(imgtwo<imgfour) = imgfour(imgtwo<imgfour);
%         imgthree(imgthree<imgfour) = imgfour(imgthree<imgfour);
        channelimg(:,:,1) = imgone;
        channelimg(:,:,2) = imgtwo;
        channelimg(:,:,3) = imgthree;
        channelimg(:,:,4) = imgfour;
        
%         
%         channelimg = uint8(channelimg);
%         %make uint8?
        
        

    else
        cd .. 
        cd('tiffs')
        ff = dir(strcat('*',ImageDetails.Channel,'*'));
%         ff = dir(strcat(ImageDetails.Channel,'*'));
%         channelspacing = round(linspace(1,length(framesForDir),9));
        channelspacing = [1 2 3 4 40 41 42 44 45];

        channelimgstack = zeros([imgsize(1) imgsize(2) length(channelspacing)]);
        for i = 1:length(channelspacing)
            t = channelspacing(i);
            channelimgstack(:,:,i) = double(loadUpTiffStackFrame(char(ff.name),t));
        end
        
        channelimg = double(loadUpTiffStackFrame(char(ff.name),t));
    %     channelimg = double(imread(char(imgfile.name)));    %load normal image
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


If = zeros(size(channelimg));

%     subdirname = char(subdir);
%         [sceneinfo,b] = regexp(subdirname,'s[0-9]+');
%         scenename = subdirname(sceneinfo:b);
%         cd(subdirname)

FinalImage = channelimgstack;
subdirname = SceneDirectoryPath;
scenename = ImageDetails.Scene;
filename = imgfile;
filename = char(filename.name);
channel = ImageDetails.Channel;

if strcmp(ImageDetails.Channel,'EGFP')
IfStack = segmentationEGFP(FinalImage,subdirname,scenename,filename,channel);   
elseif strcmp(ImageDetails.Channel,'mKate')
IfStack = segmentationNucleus(FinalImage,subdirname,scenename,filename,channel);
end


displayImageFunct(IfStack,channelimgstack,channelspacing);
updateSliders
end

function contrast_Callback(~,~)
global tcontrast lcontrast adjuster
prompt = {'High Contrast Limit','Low Contrast Limit'};
dlg_title = 'Contrast limits from 0 to 100';
num_lines = 1;
defaultans = {num2str(tcontrast),num2str(lcontrast)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
tcontrast = str2double(answer{1});
lcontrast = str2double(answer{2});

adjuster =1;
setSceneAndTime
adjuster =0;
end

function displayImageFunct(IfStack,channelimgstack,channelspacing)
global subaxes displaycomments lprcntlt prcntlt tcontrast lcontrast MainAxes displaytracking ImageDetails framesForDir prcntlz lprcntlz prcntlk lprcntlk prcntl lprcntl D ExpDate cmap cmaplz adjuster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   determine the frame to load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = strcmp(framesForDir,ImageDetails.Frame);
t = find(t==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   display the images overlayed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axes(MainAxes);
children = findobj(MainAxes,'Type','image');
delete(children);



for i=1:size(channelimgstack,3)
    If = IfStack(:,:,i);
    channelimg = channelimgstack(:,:,i);
    ax = subaxes(i);
    axes(ax);
    
%scale the image brightness
if strcmp(ImageDetails.Channel,'overlay') %when overlay display is desired
    imgone = channelimg(:,:,1);
    imgtwo = channelimg(:,:,2);
    imgthree = channelimg(:,:,3);
    imgfour = channelimg(:,:,4);
    channelimg = zeros(size(channelimg,1),size(channelimg,2),3);
        if ifCHANGEofCHANNELorSCENE==1
        
        cimgline = reshape(imgone,[1 size(imgone,1).*size(imgone,2)]);
        lprcntlz = prctile(cimgline,lcontrast);
%         prcntlz = prctile(cimgline,tcontrast) - lprcntlz;
        prcntlz = prctile(cimgline-lprcntlz,tcontrast);
        
        cimgline = reshape(imgtwo,[1 size(imgtwo,1).*size(imgtwo,2)]);
        lprcntl = prctile(cimgline,lcontrast);
%         prcntl = prctile(cimgline,tcontrast)-lprcntl;
        prcntl = prctile(cimgline-lprcntl,tcontrast);
        
        cimgline = reshape(imgthree,[1 size(imgtwo,1).*size(imgtwo,2)]);
        lprcntlk = prctile(cimgline,lcontrast);
%         prcntlk = prctile(cimgline,tcontrast)-lprcntlk;
        prcntlk = prctile(cimgline-lprcntlk,tcontrast);
        
        cimgline = reshape(imgfour,[1 size(imgtwo,1).*size(imgtwo,2)]);
        lprcntlt = prctile(cimgline,lcontrast);
%         prcntlk = prctile(cimgline,tcontrast)-lprcntlk;
        prcntlt = prctile(cimgline-lprcntlt,tcontrast);
        ifCHANGEofCHANNELorSCENE=0;
        end
    imgone = uint8(((imgone-lprcntlz)./prcntlz).*255);
    imgtwo = uint8(((imgtwo-lprcntl)./prcntl).*255);
    imgthree = uint8(((imgthree-lprcntlk)./prcntlk).*255);
    imgfour = uint8(((imgfour-lprcntlt)./prcntlt).*255);
    
        imgone(imgone<imgfour) = imgfour(imgone<imgfour);
        imgtwo(imgtwo<imgfour) = imgfour(imgtwo<imgfour);
        imgthree(imgthree<imgfour) = imgfour(imgthree<imgfour);
    
    channelimg = uint8(channelimg);
    channelimg(:,:,2) = imgone;
    If = bwperim(If);
    imgtwo(If) = 255;
    channelimg(:,:,3) = imgtwo;
    imgthree(If) = 255;
    channelimg(:,:,1) = imgthree;
    
%     If = bwperim(If);
%     channelimg(:,:,3)=uint8(If.*255);

else  %under normal circumstances
        cimgline = reshape(channelimg,[1 size(channelimg,1).*size(channelimg,2)]);
        lprcntl = prctile(cimgline,lcontrast);
%         prcntl = prctile(cimgline,tcontrast)-lprcntl;
        prcntl = prctile(cimgline-lprcntl,tcontrast);


    channelimg = uint8(((channelimg-lprcntl)./prcntl).*255);
    channelimg(channelimg == 255) =254;
    colormap(cmap);
    If = bwperim(If);
    If = imdilate(If,strel('disk',1));
    channelimg(If>0)=255;
end




himg = imagesc(channelimg);
himgax = get(himg,'Parent');
himgax.CLim = [0 256];
ttl = get(himgax,'Title');
set(ttl,'String',strcat(num2str(i),'-',ExpDate,'...',ImageDetails.Scene,'...frame ',num2str(t),' out of', num2str(length(framesForDir))));
set(ttl,'FontSize',12);
himgax.YTick = [];
himgax.XTick = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% saveChannelFiveImages
end

%% saving functions
function saveTrackingFileAs_callback(~,~)
global  SceneDirectoryPath Tracked ExportName
cd(SceneDirectoryPath)
prompt = 'filename of tracking structure to be saved?';
dlg_title = 'save tracking structure as...specific filename';
filename = char(inputdlg(prompt,dlg_title));
save(strcat(filename,ExportName,'.mat'),'Tracked')
end
function saveTrackingFileAs_callbackJ(~,~,SceneDirPath)
global   Tracked
cd(SceneDirPath)
% prompt = 'filename of tracking structure to be saved?';
% dlg_title = 'save tracking structure as...specific filename';
% filename = char(inputdlg(prompt,dlg_title));
disp(iixixixi)
%why am i saving a specific file name without asking if I want to
%overwrite?
save(strcat('finalfricktrack.mat'),'Tracked')
end
function saveTrackingFile
global  SceneDirectoryPath
cd(SceneDirectoryPath)
save('generalTrackingSavedfricktrack.mat','Tracked')
end
function saveChannelFiveImages
global If imgfile SceneDirectoryPath channelinputs
cd(SceneDirectoryPath)
fname = imgfile.name;
[a,b] = regexp(fname,channelinputs);
fname(a:b) = 'c5';

cd('c5_flat')
imwrite(If,fname,'Tiff');
cd ..

end

%% functions for determining variables
function [timeFrames,framesForDir] = determineTimeFrames(spec_directory)
dirlist = dir('_Hoechst_flat');
if isempty(dirlist)
    foldername = '_mKate_flat';
else
    foldername = '_Hoechst_flat';
end
cd (foldername)

files = dir('*.tif');
[~,~,~,chlist] = regexp([files.name],'_t[0-9]+');
numbrsCell = cellfun(@(x) str2double(x(3:end)),chlist,'UniformOutput',0);
numbrsMat = cellfun(@(x) str2double(x(3:end)),chlist,'UniformOutput',1);
numJump = min(numbrsMat)-1;
maxnumJump = max(numbrsMat);
timeFrames = cellfun(@(x) doTime(x,numJump,maxnumJump,1,chlist),numbrsCell,'UniformOutput',0);
framesForDir = cellfun(@(x) doTime(x,numJump,maxnumJump,0,chlist),numbrsCell,'UniformOutput',0);
timeJump = numJump;
end
function time = doTime(x,numJump,maxnumJump,opt,chlist)
if opt ==1
digitz = length(chlist{1})-1;%determineNumber of digits in max time number;
time = num2str(zeros(digitz,1))';
time(1) = 't';
if length(num2str(x-numJump))>1
% time(2:end) = num2str(x-numJump);
time(length(time)-length(num2str(x-numJump))+1:end) = num2str(x-numJump);
else
   time(end) = num2str(x-numJump);
end

else
digitz = length(chlist{1})-1;%determineNumber of digits in max time number;
time = num2str(zeros(digitz,1))';
time(1) = 't';
if length(num2str(x))>1
% time(2:end) = num2str(x);
time(length(time)-length(num2str(x))+1:end) = num2str(x);
else
   time(end) = num2str(x);
end

end

end
function ImageDetails = InitializeImageDetails

ImageDetails.Scene=[];
ImageDetails.Channel=[];
ImageDetails.Frame=[];

end


%% load up images
function FinalImage=loadUpTiffStack(filenames)
cfile = filenames;
FileTif = char(cfile);
InfoImage=imfinfo(FileTif);

mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
    for i = 1:NumberImages;
       TifLink.setDirectory(i);
       FinalImage(:,:,i)=TifLink.read();
    end
end
function FinalImage=loadUpTiffStackFrame(filenames,frame)
cfile = filenames;
FileTif = char(cfile);
InfoImage=imfinfo(FileTif);

mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,1,'uint16');
 
TifLink = Tiff(FileTif, 'r');
    i = frame;
       TifLink.setDirectory(i);
       FinalImage=TifLink.read();
    
end
