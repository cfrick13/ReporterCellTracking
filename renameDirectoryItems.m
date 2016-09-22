% currentdir = pwd;
cd (currentdir)
clear
currentdir = pwd;
folderlist = dir('*scene*');
foldernames = {folderlist.name};

for folder = foldernames
    folder =char(folder);
    cd (folder)
    fdir = pwd;
    cd('c1_flat')
    filelist = dir('*tif');
    filenames = {filelist.name};
    
    for fname = filenames
        fname = char(fname);
        [a,b] = regexp(fname,'c2');
        c = fname(a:b);
        
        if ~isempty(c)
            c(2) = '1';
            fnamenew = fname;
            fnamenew(a:b) = c;
            disp(fnamenew)
            movefile(fname,fnamenew);
        end
        
    end
    
    clear c
    
    cd (fdir)
    cd ('c2_flat')
    
    filelist = dir('*tif');
    filenames = {filelist.name};
    
    for fname = filenames
        fname = char(fname);
        [a,b] = regexp(fname,'c1');
        c = fname(a:b);
        
        if ~isempty(c)
            c(2) = '2';
            fnamenew = fname;
            fnamenew(a:b) = c;
            disp(fnamenew)
            movefile(fname,fnamenew);
        end
        
    end
    clear c
    cd(currentdir)
    
end