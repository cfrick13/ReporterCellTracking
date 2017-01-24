function SortBKGandFLATcorrectedImagesIntoFolders(~,datename,channelstoinput)


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
    cd ('flatfield_corrected')

    %%
    folderlist = dir('*_s*');
        for folder = {folderlist.name} %set the scene directory for all
        cd(folder{1})
        mkdir('tiffs');
        %check if directories have been made previously
        %make individual directories
            for i = 1:length(channelstoinput)
                dd = dir(strcat(channelstoinput{i},'_flat'));
                    if isempty(dd) 
                        nodir = 1;
                        filelist = dir(strcat('*_',channelstoinput{i},'*.tif')); 
                            if isempty(filelist) %if no files are present do not make the directory
                            else
                            mkdir(strcat(channelstoinput{i},'_flat'))
                            end
                    else
                        nodir = 0;
                    end

                foldername = strcat(channelstoinput{i},'_flat');
                filelist = dir(strcat('*_',channelstoinput{i},'*.tif'));  
                    if isempty(filelist) && nodir==0
                        cd(foldername)
                        filelist = dir(strcat('*_',channelstoinput{i},'*.tif'));  

                        if isempty(filelist) 
                            cd .. 
    %                         disp(foldername)
                            ddd = dir(foldername);
                            if isempty(ddd)
                            else
                                disp(foldername)
                                disp(ddd)
                                pause(0.5)
                            rmdir(foldername);
                            end
                        else
                            cd ..
                        end

                    else

                        for cfile={filelist.name}
                           filepath = which(cfile{1});
                           disp(filepath);
                           movefile(filepath,foldername)
    %                        cd(foldername)
    %                        movefile(filepath)
    %                        cd .. 
                        end
                    end
            end
            cd ..
        end
end
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


