function makeNiceMovie(background)
dirpath = cell(1,2);
dirpath{1} = 'H:\2016_09_28 screen';
folder = 'flatfield_corrected\2016_09_28 screen_scene_s01';
dirpath{2} = strcat(dirpath{1},'\',folder);
cd(dirpath{2});

filelist = dir;
fnames = {filelist.name};
[a,b,c,d] = regexpi(fnames,'egfpstack');
filenames = fnames{~cellfun(@isempty,d)};
egfpstack = loadUpTiffStack(filenames);

egfpstacksub = zeros(size(egfpstack));
for i=1:size(egfpstacksub,3)
egfpstacksub(:,:,i) = egfpstack(:,:,i)-background(i);
end

mkdir('background subtracted')
cd('background subtracted')
for i=1:size(egfpstacksub,3)
    imwrite(uint16(egfpstacksub(:,:,i)),strcat('egfp_bsub_',num2str(i),'.tif'),'Tiff');
end

end

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
