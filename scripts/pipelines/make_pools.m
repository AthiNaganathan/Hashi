function make_pools(pools, no_of_microstates)
    % Pool output structures from RANCH into user-defined pools
    foldername = fullfile(pwd, "..", "results", "microstates");
    folders = dir(foldername);
    folders = folders(3:end);
    n = length(folders);

    pool_dir = fullfile(pwd, "..", "results", "pools");
    mkdir(pool_dir);
    rmdir(pool_dir, 's');
    mkdir(pool_dir);

    % Move ranch files into pools
    folderindex = 0;
    for i = 1:size(pools, 1)
        curpoolpath = fullfile(pool_dir, "pool_"+string(i));
        mkdir(curpoolpath);
        rmdir(curpoolpath, 's');
        mkdir(curpoolpath);

        for j = 1:no_of_microstates
            folderindex = folderindex + 1;
            curfolder = folders(folderindex);

            ranchPath = fullfile(curfolder.folder, curfolder.name, 'ranch');
            cifFiles = dir(fullfile(ranchPath, '*.cif'));
            for k = 1:numel(cifFiles)
                srcFile = fullfile(cifFiles(k).folder, cifFiles(k).name);
                copyfile(srcFile, curpoolpath);
            end
        end
    end
    
    % Move microstate info files for future reference
    for i = 1:n
        folders_list = dir(fullfile(foldername, sprintf("%05d", i)+"_*"));
        req_path = fullfile(folders_list(1).folder, folders_list(1).name);
        pool_path = fullfile(pool_dir, "pool_" + string(floor((i-1)/no_of_microstates) + 1));
        copyfile(fullfile(req_path, "microstate_info.mat"), fullfile(pool_path, string(i)+".mat"))
    
        if mod(i, no_of_microstates) == 0
            fileList = dir(fullfile(pool_path, "*.mat"));
            numFiles = numel(fileList);
            collectedData = cell(numFiles, 1);
            
            for k = 1:numFiles
                temp = load(fullfile(pool_path, fileList(k).name));
                collectedData{k} = temp.mis; 
            end
            
            % 5. Save the combined variable(s) to a new .mat file
            % 'finalMatrix' will be the variable name in the new file
            save(fullfile(pool_dir, "pool_"+string(floor((i-1)/no_of_microstates) + 1)+".mat"), 'collectedData');
            delete(fullfile(pool_path, '*.mat'));
        end
    end
    disp("Done")
end