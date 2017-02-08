function  SegmentationOfTimeLapseImages(~,datename,BACKGROUND)
global pStruct

pStruct = struct();
parameterDefaults = [106 106 40;...
                    0.5 0.5 1.3;...
                    20 20 5;...
                    1 1 1];
                    
parameterStrings = {'nucDiameter','threshFactor','sigmaScaledToParticle','noparametercurrently'};
channelList = {'DIC','EGFP','mKate'};

for p = 1:length(parameterStrings)
    pString = char(parameterStrings{p});
    for c = 1:length(channelList)
        cString = char(channelList{c});
        pStruct.(cString).(pString) = parameterDefaults(p,c); 
    end
end




choosefiles(datename,BACKGROUND)
end

function choosefiles(datename,BACKGROUND)
global nucleus_seg 
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


    primarydir = strcat(parentdir,datename,'/','flatfield_corrected');
    cd(primarydir)

    nucleus_seg = '_mKate_flat';
    % nucleus_seg = 'CFP';
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
            dirlist = dir(nucleus_seg);
            if isempty(dirlist)
                dirlist = dir('mKate_flat');
                foldername = 'mKate_flat';
            else
                foldername = nucleus_seg;
            end
            cd(finaldirname)

                file = dir(strcat('*',foldername,'*tif'));
                filename = char({file.name});
                FinalImage = loadStack(filename);
                cd ..
                cd (foldername)
                filelist = dir('*.tif');
                fname = filelist(1).name;
                cd ..
                disp(scenename)
                segmentationNucleus(FinalImage,'mKate',scenename,fname,'NucleusBinary_flat');

            dirlist = dir('mKate_flat');
            if isempty(dirlist)
                dirlist = dir(nucleus_seg);
                foldername = nucleus_seg;
            else
                foldername = 'mKate_flat';
            end   
            cd(finaldirname)
                file = dir(strcat('*',foldername,'*tif'));
                filename = char({file.name});
                FinalImage = loadStack(filename);
                    cd ..
                    cd (foldername)
                    filelist = dir('*.tif');
                    fname = filelist(1).name;
                    cd ..
            %     segmentationRFP(FinalImage,subdirname,scenename,fname,'mKatebinary_flat');


            dirlist = dir('_EGFP_flat');
            if isempty(dirlist)
                dirlist = dir('mKate_flat');
                foldername = 'mKate_flat';
            else
                foldername = '_EGFP_flat';
            end
            cd(finaldirname)

            file = dir(strcat('*',foldername,'*tif'));
            filename = char({file.name});
            FinalImage = loadStack(filename);
            cd ..
            cd (foldername)
            filelist = dir('*.tif');
            fname = filelist(1).name;
            cd ..
            segmentationMNG(FinalImage,subdirname,scenename,fname,'EGFPbinary_flat');
            stophere=1;
        %     segmentationBKGsecond(FinalImage,subdirname,scenename,fname,'BKGbinary_flat');
            segmentationEGFP(FinalImage,'EGFP',scenename,fname,'BKGbinary_flat');
            cd ..
        end

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






function If = segmentationNucleus(FinalImage,channel,scenename,filename,segchannel)
global  pStruct 
% cd(subdirname)
mkdir(strcat(segchannel));

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
    imgRawDenoised = double(imgRaw);


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
            Isum = Isum+Ier;
        end
        Isum(Isum>nucDiameter) = nucDiameter;
         waterBoundary = imerode(Im,strel('disk',1));




    % I = imgRawDenoised;
    % I = gaussianBlurz(rawMinusLPScaled,sigma./4,kernelgsize);
    I = rawMinusLPScaled;

    %gradmag
    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);

    %Smoothing
    I = Isum;
    se = strel('disk', 10);
    Io = imopen(I, se);
    % se = strel('disk', ceil(nucDiameter./6));
    Ie = imerode(I, se);
    Ieg = gaussianBlurz(Ie,sigma./2,kernelgsize);
    Iobr = imreconstruct(Ie, I);
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    fgm = imregionalmax(Ieg);
    % fgm = imregionalmax(Iobrcbr);
    se2 = strel(ones(2,2));
    fgm2 = imclose(fgm, se2);
    fgm3 = imerode(fgm2, se2);
    % fgm4 = bwareaopen(fgm3, 2);
    fgm4 = imdilate(fgm,strel('disk',5));
    % bw = imbinarize(Iobrcbr);
    bw = Im;
    D = bwdist(bw);
    DL = watershed(D,4);
    bgm = DL == 0;
    gradmag2 = imimposemin(gradmag, bgm | fgm4);

    L = watershed(gradmag2,8);
    L(waterBoundary<1) = 0;
    If = L>1;

    
    %remove segmentation if it overlaps with the border of the image
    Cellstruct = bwconncomp(If);
    PX = Cellstruct.PixelIdxList;
    testimage = zeros(size(If));
    edgesize=0;
    testimage(1:(1+edgesize),1:512) =1; testimage((512-edgesize):512,1:512) =1; testimage(1:512,(1:1+edgesize)) =1; testimage(1:512,(512-edgesize):512) =1;
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


        time = tsn{frames};
%         tim = time(2:end);
%         IfFinal(:,:,frames)=If;
        savethatimage(scenename,time,If.*255,frames,filename,segchannel)
    end

stophere=1;
end

function If = segmentationNucleusOld(FinalImage,subdirname,scenename,filename,channel)
global nucleus_seg
fig=1;
mkdir(strcat('NucleusBinary_flat'));
mkdir(strcat('c5_flat'));

% parameters
dimdiff = 2048./size(FinalImage(:,:,1),1);
zerostrel = 5;
firststrel = round(50./(dimdiff.^2));
sigmafirst = firststrel.*5;
kernelgsizefirst = firststrel.*10;
% fracsmoothing = 0.5.*dimdiff;
fracsmoothing = 0.5;
weiner2p = 20;

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


img = wiener2(img,[weiner2p weiner2p]);
se =strel('disk',zerostrel);
Ie = imerode(imgorig,se);
Iobr = imreconstruct(Ie,img);
Iobrone = Iobr;
Iobrd = imdilate(Iobr,se);
Iobrdone = Iobrd;
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
IobrcbrF = imcomplement(Iobrcbr);
gaus = double(IobrcbrF);

% figure(99)
% subplot(3,3,1);
% imagesc(Iobr);
% subplot(3,3,2);
% imagesc(Iobrd);
% subplot(3,3,3);
% imagesc(Iobrcbr);
% subplot(3,3,4);
% imagesc(Iobrcbr);
% subplot(3,3,5);
% imagesc(IobrcbrF);
% subplot(3,3,6);
% imagesc(gaus);

se =strel('disk',firststrel);
Ie = imerode(gaus,se);
Iobr = imreconstruct(Ie,gaus);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
gaus = double(Iobrcbr);
% 
% figure(999)
% subplot(3,3,1);
% imagesc(Iobr);
% subplot(3,3,2);
% imagesc(Iobrd);
% subplot(3,3,3);
% imagesc(Iobrcbr);
% subplot(3,3,4);
% imagesc(Iobrcbr);
% subplot(3,3,5);
% imagesc(IobrcbrF);
% subplot(3,3,6);
% imagesc(gaus);
% subplot(3,3,7)
% imagesc(imgorig)


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
left=0.5*mf;
slopedown=0.4*mf;
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

subtracted = sub_scale_corr-subtractionthreshold.*1.1;
% subtracted = sub_scale_corr-subtractionthreshold.*1.1;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;





Ih = Ie>0;
Ihd = imdilate(Ih,strel('disk',1));
Ihdc = imclose(Ihd,strel('disk',4));
Ihdcf = imfill(Ihdc,'holes');
Im = Ihdcf;



%%%%% this is the ultimate addition for watershed segmentation!!!
see = strel('disk',1);
seo = strel('disk',8);
Ier = Im;
Isum = Ier;
for i=1:30
    Ier = imerode(Ier,see);
    Iero = imopen(Ier,seo);
%     Isum = Isum+(Iero.*i);
    Isum = Isum+(Iero);
    Ier=Iero;
end
Isum(Isum>130) = 130;
figure(66)
imagesc(Isum);
gausshed = gaussianBlurz(Isum,round(sigmafirst./dimdiff),round(kernelgsizefirst./dimdiff));
imgt = -double(gausshed);
waterBoundary = imerode(Im,strel('disk',1));
imgt(~(waterBoundary>0)) = -Inf;
L=watershed(imgt);

L(waterBoundary<1) = 0;
If = L>1;
If = imerode(If,strel('disk',2));




time = tsn{frames};
tim = time(2:end);
if frames==10
figure(1)
imagesc(If)
stophere=1;
end
savethatimage(scenename,time,If.*255,frames,filename,channel)
end


end

function If = segmentationRFP(FinalImage,subdirname,scenename,filename,channel)
fig=1;
mkdir(strcat(channel));

% parameters
left = 0.004;
slopedown = 0.003;


dimdiff = 2048./size(FinalImage(:,:,1),1);

zerostrel = 2;
firststrel = round(30./(dimdiff.^2));
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;
% fracsmoothing = 0.5.*dimdiff;
fracsmoothing = 0.5;

foldername = 'mKate_flat';
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
left=0.5*mf;
slopedown=0.4*mf;
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

subtracted = sub_scale_corr-subtractionthreshold.*1.04;
% subtracted = sub_scale_corr-subtractionthreshold.*1.1;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;


if frames==10
figure(2)
imagesc(Ie)
stophere=1;
end

%%%%%%%%%%%%%%%
Ih = Ie>0;
se = strel('disk',5);
Ihe = imerode(Ih,se);
Ihed = imdilate(Ihe,se);
If = Ihed;



stophere=1;
% time = settimecharacter(frames);
time = tsn{frames};
tim = time(2:end);
savethatimage(scenename,time,If.*255,frames,filename,channel)
end


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

function If = segmentationEGFP(FinalImage,channel,scenename,filename,segchannel)
global  pStruct nucleus_seg
mkdir(strcat(segchannel));

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
    imgRawDenoised = double(imgWWW);
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
%     Ihcf = imfill(Ihc,'holes');
%     Ihcf = Ihc;
    Ihcd = imdilate(Ihc,strel('disk',20));
%     Ihcfd = Ihcf;
    Im=Ihcd;
    If =Im;





    
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
function segmentation(FinalImage,subdirname,scenename,filename,channel)
fig=1;
mkdir(strcat('NucleusBinary_flat'));
mkdir(strcat('c5_flat'));
% mkdir(strcat('NucleusBinary_flat'));
for frames = 1:size(FinalImage,3)
% for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 

%% parameters
strelsize           =    3;      %3
peakthresh          =    200;    %set low...will be corrected by the segthresh.
sigma               =    20;          %80
kernelgsize         =    15;   %15 for 512x512
meanregion          =    2;
segthresh           =    2200;
strelsizegaus       =    3;  %10
strelsizesub        =    2;
subtractionthreshold =   1400;
Upz = 0.07;
Downz = 0.06;

mask_em         = zeros(size(img));

%% start

%smooth initial image
se = strel('disk',strelsize);
Ie = imerode(img,se);
Iobr = imreconstruct(Ie,img);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
img = Iobrcbr;


%     figure(fig)
%     fig = fig+1;
%     imagesc(img);
%     title('smoothed input image');


            %ONLY NECESSARY IF TRYING TO FIND THE PEAKS IN THE SMOOTHED GAUS 
            
            gaus = gaussianblur(img);
           
    % figure(fig)
    % fig = fig+1;
    % hold off
    % bar(bincenters,fraction)
    % ylim([0 0.01])
    % ylim([0 0.003])
    % xlim([0 6000])
    % hold on
    % stem(segthresh,1,'g');
    % hold off


% findpeaksgaus(gaus) % homemade function to find peaks after applying gaussian blur 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub = double(img) -double(gaus);%%%%%%% key step!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div = double(gaus)./double(img);
maxlocale = max(max(div(10:500,10:500)));
d = find(div == maxlocale);
div(d);
positiond = d(1);
[e,f] = ind2sub(size(img),positiond);


 minlocale = min(min(sub(10:500,10:500)));
a = find(sub == minlocale);
sub(a);

position = a(1);
[b,c] = ind2sub(size(img),position);

%     figure(fig)
%     fig = fig+1;
%     imagesc(sub);
%     title('subtracted');
%     hold on
%     plot(c,b,'y+');


displaynormalizationfactor(e,f,frames,scenename)
%normalizationfactors = double(gaus(b,c))./double(img(b,c));
normalizationfactors = double(gaus(e,f))./double(img(e,f));
normalizationfactor = mean(mean(normalizationfactors));

% figure(fig)
% fig = fig+1;
% imagesc(subgauss);hold on
% plot(c,b,'y+')

scaledgaus = double(gaus)./normalizationfactor;
sub_scale_corr = double(img)-double(scaledgaus);

 
%     figure(fig)
%     fig = fig+1;
%     imagesc(sub_scale_corr);
%     title('subtracted_scaled');


subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:2:10000);
numbers = medfilt1(numbers, 15); %smooths curve
fraction = numbers./sum(numbers);
% rightedge = find(fraction > 0.003,1,'last');
% slopedown = find(fraction(rightedge:end) <0.002,1,'first');
mf = max(fraction);
Upz = mf.*0.8;
Downz = mf.*0.7;
rightedge = find(fraction > Upz,1,'last');
slopedown = find(fraction(rightedge:end) <Downz,1,'first');
threshlocation=[];
if ~isempty(slopedown)
thresh = bincenters(rightedge+slopedown);
threshlocation = thresh - 0.05*thresh;
else
    if ~isempty(threshlocation)
    threshlocation = threshlocation;
    else
        threshlocation = 300;
    end
end

%  figure(fig)
% fig = fig+1;
% bar(bincenters,fraction);hold on
% ylim([0 0.01])
% xlim([0 2000]);hold on
% stem(threshlocation,1,'r');
%     hold off
%     
    subtractionthreshold = threshlocation;


subtracted = sub_scale_corr-ones(size(sub_scale_corr)).*subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);

%     figure(fig)
%     fig = fig+1;
%     imagesc(subtractedzero,[0 20]);
%     title('subtracted_scaled_zeroed');


se = strel('disk',strelsizesub);
Ie = imerode(subtractedzero,se);
a = find(Ie>1);
submax = zeros(size(Ie));
Ie(a)=50;
% 
%     figure(fig)
%     fig=fig+1;
%     imagesc(Ie);
%     title('subtractedzeromax');
%     fig =1;


time = settimecharacter(frames);
        
stophere=1;
savethatimage(scenename,time,Ie,frames,filename,channel)
% imwrite(uint8(Ie),strcat(scenename,'_','t',time,'_NucleusBinary_flat.tif'));
end
stophere=1;
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
