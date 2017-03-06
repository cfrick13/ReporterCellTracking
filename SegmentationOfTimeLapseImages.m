function  SegmentationOfTimeLapseImages(~,datename,BACKGROUND)
global nucleus_seg foldernameglobal

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
exportdir = strcat(parentdir,'Tracking/Export/');




pStruct = struct();
parameterDefaults = [106 100 40;...  % nuc diameter [dic egfp mkate]
                    0.5 0.5 1;... % threshfactor [dic egfp mkate]
                    20 15 2;... % sigmaScaledToParticle [dic egfp mkate]
                    1 1 1]; % noparameter currently [dic egfp mkate]
parameterStrings = {'nucDiameter','threshFactor','sigmaScaledToParticle','noparametercurrently'};
channelList = {'DIC','EGFP','mKate'};
for p = 1:length(parameterStrings)
    pString = char(parameterStrings{p});
    for c = 1:length(channelList)
        cString = char(channelList{c});
        pStruct.(cString).(pString) = parameterDefaults(p,c); 
    end
end
pStruct = loadSegmentParameters(pStruct,datename,exportdir); %loads saved value of pStruct





%find directories with images
    for expdircell = explist
        expdirname = char(expdircell);


    primarydir = strcat(parentdir,datename,'/','flatfield_corrected');
    cd(primarydir)
    


    primarylist = dir('*_s*');
    % primarylist = dir('*_s02');
    subd = {primarylist.name};
        for i = 1:length(BACKGROUND)
            bkstr = num2str(BACKGROUND(i)); 
            if length(bkstr)>1
                bkarray{i} = strcat('s',bkstr);
            else
                bkarray{i} = strcat('s0',bkstr); 
            end
        end
    bkinputs =channelregexpmaker(bkarray);
    [~,~,~,d] = regexp(subd,bkinputs);
    subdirlog = cellfun(@isempty,d,'UniformOutput',1);
    subd = subd(subdirlog);

    % parfor i=1:length(subd)
        for i=1:length(subd)
            subdir = subd{i};
            subdirname = char(subdir);
                [sceneinfo,b] = regexp(subdirname,'s[0-9]+');
                scenename = subdirname(sceneinfo:b);
                cd(subdirname)
            finaldir = dir('*tiffs*');
                finaldirname = char({finaldir.name});
                cd(finaldirname)
            % file = dir('*EGFP_flat.tif*');

            cd .. 
            folderlist = dir('*mKate*');
            foldernameglobal=folderlist.name;
            nucleus_seg = foldernameglobal;
            
            
            
            dirlist = dir(nucleus_seg);
            if isempty(dirlist)
                dirlist = dir(foldernameglobal);
                foldername = foldernameglobal;
            else
                foldername = nucleus_seg;
            end
            
                cd(finaldirname)
                file = dir(strcat('*',foldername,'*tif'));
                filename = char({file.name});
                FinalImage = double(loadStack(filename));
                cd ..
                cd (foldername)
                filelist = dir('*.tif');
                fname = filelist(1).name;
                cd ..
                disp(scenename)
                segmentationNucleus(FinalImage,'mKate',scenename,fname,'NucleusBinary_flat',pStruct);
                

cd('NucleusBinary_flat')
filelist = dir('*.tif');
filenames = {filelist.name};
FinalImage = zeros(size(FinalImage,1),size(FinalImage,2),length(filenames));
for jkl = 1:length(filenames)
    fname = char(filenames{jkl});
    FinalImage(:,:,jkl) = double(imread(fname));
end
cd ..
% segmentationREPORTERBKG(FinalImage,'mKate',scenename,fname,'BKGbinary_flat',pStruct);
                
            dirlist = dir(foldernameglobal);
            if isempty(dirlist)
                dirlist = dir(nucleus_seg);
                foldername = nucleus_seg;
            else
                foldername = foldernameglobal;
            end   
            cd(finaldirname)
                file = dir(strcat('*',foldername,'*tif'));
                filename = char({file.name});
                FinalImage = double(loadStack(filename));
                    cd ..
                    cd (foldername)
                    filelist = dir('*.tif');
                    fname = filelist(1).name;
                    cd ..
            %     segmentationRFP(FinalImage,subdirname,scenename,fname,'mKatebinary_flat');


            dirlist = dir('_EGFP_flat');
            if isempty(dirlist)
                dirlist = dir(foldernameglobal);
                foldername = foldernameglobal;
            else
                foldername = '_EGFP_flat';
            end
            cd(finaldirname)

            file = dir(strcat('*',foldername,'*tif'));
            filename = char({file.name});
            FinalImage = double(loadStack(filename));
            cd ..
            cd (foldername)
            filelist = dir('*.tif');
            fname = filelist(1).name;
            cd ..
%             segmentationMNG(FinalImage,subdirname,scenename,fname,'EGFPbinary_flat',pStruct);
            stophere=1;
%             segmentationBKGsecond(FinalImage,subdirname,scenename,fname,'BKGbinary_flat');
            segmentationEGFP(FinalImage,'EGFP',scenename,fname,'BKGbinary_flat',pStruct);
            
            cd ..
        end

    end
end


function pStruct = loadSegmentParameters(pStruct,datename,exportdir)


cd(exportdir)
filename = strcat('*',datename,'*segmentParameters*');

filelist = dir(filename);
if ~isempty(filelist)
loadname = char((filelist.name));
A = load(loadname); %load pstruct values
pStruct = A.pStruct;
end

    

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






function If = segmentationNucleus(FinalImage,channel,scenename,filename,segchannel,pStruct)
global foldernameglobal
% cd(subdirname)
dirlist = dir(strcat(segchannel));
if isempty(dirlist)
mkdir(strcat(segchannel));
end



foldername = foldernameglobal;
% foldername = foldernameglobal;
tsn = determineTimeFrame(foldername);

% channel ='mKate';
% parameters
nucDiameter = pStruct.(channel).nucDiameter;
threshFactor = pStruct.(channel).threshFactor;
sigmaScaledToParticle = pStruct.(channel).sigmaScaledToParticle;
finalerode=2;



% prepareCcodeForAnisotropicDiffusionDenoising(denoisepath)

%start
IfFinal = zeros(size(FinalImage));
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
    weinerP=5;
    imgRawDenoised = wiener2(imgRaw,[weinerP weinerP]);


    
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
    lowperc = prctile(vec,1);
    highperc = prctile(vec,100);
    [numbers,bincenters] = hist(vec,lowperc:(highperc-lowperc)/1000:highperc);
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
    Ihe = imerode(Ih,strel('disk',2));
    Ihed = imdilate(Ihe,strel('disk',2));
    Ihc = imclose(Ihed,strel('disk',2));
    Ihcf = imfill(Ihc,'holes');
    Im=Ihcf;

    

    
    
%%%%% this is the ultimate addition for watershed segmentation!!!
    see = strel('disk',1);
    Isum = Im;
    Ier = Isum;
%     figure(2)
    for i=1:round((nucDiameter/2))
        Ier = imerode(Ier,see);
        Isum = Isum+Ier;
    end
    Isum(Isum>nucDiameter) = nucDiameter;
    waterBoundary = imerode(Im,strel('disk',1));
    
    
    
    
% I = imgRawDenoised;
% I = gaussianBlurz(rawMinusLPScaled,sigma./4,kernelgsize);
I = rawMinusLPScaledContrasted;

%gradmag
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

%Smoothing
I = Isum;
width = round(nucDiameter./10);
se = strel('disk', width);
Io = imopen(I, se);
Ie = imerode(Io, se);
Ieg = gaussianBlurz(Ie,round(sigma./2),round(kernelgsize./2));
%     width = round(nucDiameter./10);
%     Ime = imerode(Ihcf,strel('disk',width));
%     Imeo = imopen(Ime,strel('disk',width));
%     Ieg(~Imeo)=0;
fgm = imregionalmax(Ieg);
width = round(nucDiameter./10);
fgm4 = imdilate(fgm,strel('disk',width));
% fgm4 =fgm;
bw = Im;
D = bwdist(bw);
DL = watershed(D,4);
bgm = DL == 0;
gradmag2 = imimposemin(gradmag, bgm | fgm4);

L = watershed(gradmag2,8);
L(waterBoundary<1) = 0;
If = L>1;


CellObjects = bwconncomp(If,8);
PX = CellObjects.PixelIdxList;
pxl = cellfun(@length,PX,'UniformOutput',1);
pxlog = pxl>((pi.*((nucDiameter).*2)));
PXX = PX(~pxlog);
% CellObjects.PixelIdxList = PXX;
% CellObjects.NumObjects = length(PXX);
If(vertcat(PXX{:})) = 0;




%%%%%%%%%remove segmentation if it overlaps with the border of theimage%%%%%%%%%%%%%%%%%%%
    Cellstruct = bwconncomp(If);
    PX = Cellstruct.PixelIdxList;
    testimage = zeros(size(If));
    edgesize=0;
    dim = size(FinalImage);
    testimage(1:(1+edgesize),1:dim(2)) =1; testimage((dim(1)-edgesize):dim(1),1:dim(2)) =1; testimage(1:dim(1),(1:1+edgesize)) =1; testimage(1:dim(1),(dim(2)-edgesize):dim(2)) =1;
    borderstruct = bwconncomp(testimage);
    pxborder = borderstruct.PixelIdxList;
    borderpixels  = pxborder{1};
%     stophere=1;
    for pidx = 1:length(PX)
        px = PX{pidx};
        testlog = ismember(px,borderpixels);
%         sumtest(pidx) = sum(testlog);
        if sum(testlog)>20
            If(px) = 0;
        end
        
    end


    % figure(3)
    % 
    % subplot(1,7,1);imagesc(I)
    % subplot(1,7,2);imagesc(Io)
    % subplot(1,7,3);imagesc(Ie)
    % subplot(1,7,4);imagesc(Ieg);
    % subplot(1,7,5);imagesc(bgm | fgm4);
    % subplot(1,7,6);imagesc(gradmag2);
    % subplot(1,7,7);imagesc(If);




        time = tsn{frames};
        savethatimage(scenename,time,If.*255,frames,filename,segchannel)
    end

stophere=1;
end

function If = segmentationMNG(FinalImage,subdirname,scenename,filename,channel,pStruct)
global nucleus_seg foldernameglobal
fig=1;

dirlist = dir(strcat(channel));
if isempty(dirlist)
mkdir(strcat(channel));
end

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
    dirlist = dir(foldernameglobal);
    foldername = foldernameglobal;
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

function If = segmentationEGFP(FinalImage,channel,scenename,filename,segchannel,pStruct)
global   nucleus_seg foldernameglobal

dirlist = dir(strcat(segchannel));
if isempty(dirlist)
mkdir(strcat(segchannel));
end

foldername = foldernameglobal;
% foldername = foldernameglobal;
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

%initial segmentation to determine how much of image is covered by cells
    img = FinalImage(:,:,1); 
    imgRaw = gaussianBlurz(double(img),ceil(sigma./10),ceil(kernelgsize./10));

    imgW = wiener2(img,[1 20]);
    imgWW = wiener2(imgW,[20 1]);
    imgWWW = wiener2(imgWW,[5 5]);
    imgRawDenoised = imgWWW;
    denoiseVec = double(reshape(imgRawDenoised,size(imgRawDenoised,1)^2,1));
    highpoints = prctile(denoiseVec,95);
    imgRawDenoised(imgRawDenoised>highpoints) = highpoints;
    %
    imgLowPass = gaussianBlurz(double(imgRawDenoised),sigma,kernelgsize);
    rawMinusLP = double(imgRawDenoised) -double(imgLowPass);%%%%%%% key step!
    rawMinusLPvec = reshape(rawMinusLP,size(rawMinusLP,1)^2,1);
    globalMinimaValues = prctile(rawMinusLPvec,0.01);
    globalMinimaIndices = find(rawMinusLP < globalMinimaValues);
    LPscalingFactor = imgRawDenoised(globalMinimaIndices)./imgLowPass(globalMinimaIndices);
    imgLPScaled = imgLowPass.*nanmedian(LPscalingFactor);
    rawMinusLPScaled = double(imgRawDenoised) - double(imgLPScaled);
    %
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
    Ih = imclose(Ih,strel('disk',20));
    areaOfSegmentation = sum(sum(Ih));
    %
    percentageOfImageSegmented = round(100*(areaOfSegmentation./(size(img,1)*size(img,2))));
    disp(percentageOfImageSegmented);
    if percentageOfImageSegmented==100
        percentageOfImageSegmented=99;
    end
%     percentageOfImageSegmented=10;
%     nucDiameter = nucDiameter.*(percentageOfImageSegmented/50);
% percentageOfImageSegmented=90;


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
    highpoints = prctile(denoiseVec,percentageOfImageSegmented);
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
    width = 10;
    Ihc = imclose(Ih,strel('disk',width));
%     Ihcf = imfill(Ihc,'holes');
%     Ihcf = Ihc;
%     Ihcd = imdilate(Ihc,strel('disk',width));
%     Ihcfd = Ihcf;
    Im=Ihc;
    If =Im;


If = imgRawDenoised;
mmIf = max(max(If)) ;
If(If<mmIf)=0;
If(If == mmIf)=1;
If = logical(If);
stophere=1;
arealimit = (100-percentageOfImageSegmented)./8;
imgarea = (size(If,1).*size(If,2));

   a = length(If==0);
   width = 10;
   Ig= If;
   while  1
       se = strel('disk',width);
       a = ((imgarea-sum(sum(Ig)))./imgarea).*100;
       if a<arealimit
           break
       else
            Ig = imdilate(Ig,se);
       end
%        disp(a)
   end
%     imagesc(Ig)

    If=Ig;





    
    time = tsn{frames};
    tim = time(2:end);
%     IfFinal(:,:,frames)=If;
    savethatimage(scenename,time,If.*255,frames,filename,segchannel)
end


stophere=1;
end


function If = segmentationREPORTERBKG(FinalImage,channel,scenename,filename,segchannel,pStruct)
global   nucleus_seg foldernameglobal
mkdir(strcat(segchannel));

foldername = foldernameglobal;
% foldername = foldernameglobal;
tsn = determineTimeFrame(foldername);

% parameters
nucDiameter = pStruct.(channel).nucDiameter;
threshFactor = pStruct.(channel).threshFactor;
sigmaScaledToParticle = pStruct.(channel).sigmaScaledToParticle;
kernelgsize = nucDiameter; %set kernelgsize to diameter of nuclei at least
sigma = nucDiameter./sigmaScaledToParticle; %make the sigma about 1/5th of kernelgsize


for frames=1:size(FinalImage,3)
img = FinalImage(:,:,frames);
If = logical(img);
se = strel('disk',20);
Ifd = imdilate(If,se);
If = ~Ifd;


    
    time = tsn{frames};
    tim = time(2:end);
%     IfFinal(:,:,frames)=If;
    savethatimage(scenename,time,If.*255,frames,filename,segchannel)
end


stophere=1;
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
