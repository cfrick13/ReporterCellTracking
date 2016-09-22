function renamemWRONGtoRIGHT(A,B)
C = B;

wrongname = {'HEX','Dic'};
rightname = {'EGFP','DIC'};

for indi =1:length(wrongname)
wrongnamedir{indi} = strcat('*',wrongname{indi},'*.tif');
rightnamedir{indi} = strcat('*',rightname{indi},'*.tif');
end


E={B}; %beginning of file name for saving
channelstoinput = {'_mKate','_EGFP','_CFP','DIC'};
% channelstoinput = {'mKate','_EGFP','_CFP','_DIC'};
channelinputs = '(';
for i=1:length(channelstoinput)
    if i ==1
    channelinputs = strcat(channelinputs,channelstoinput{i});
    elseif i < length(channelstoinput)
        channelinputs = strcat(channelinputs,'|',channelstoinput{i});
    else
        channelinputs = strcat(channelinputs,'|',channelstoinput{i},')');
    end
end



BACKdir = strcat(A,B,'\back');mkdir (BACKdir);
DIRone = strcat(A,B,'\',C);

cd(DIRone);
for indi = 1:length(wrongnamedir)
    wnamedir = char(wrongnamedir{indi});
dirlist = dir(wnamedir);
cfile = {dirlist.name};
    wname = char(wrongname{indi});
    rname = char(rightname{indi});
% parfor i = 1:length(cfile)
for i = 1:length(cfile)
    oldfilename = cfile(i);
    oldfile = oldfilename{1};
    [a,b] = regexp(oldfile,wname);
    newfile = strcat(oldfile(1:a-1),rname,(oldfile(b+1:end)));
    movefile(oldfile,newfile);
    disp(i)
end
end


end