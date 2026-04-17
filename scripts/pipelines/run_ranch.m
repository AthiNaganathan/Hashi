function run_ranch(pname, reps)     
    % Calls and runs RANCH from terminal
    foldername = fullfile("..", "results", "microstates");
    folders = dir(foldername);
    n = length(folders);
    reps = string(reps);
    
    for i = 3:n
        folder_path = fullfile(foldername, folders(i).name);
        ranch_path = fullfile(folder_path, "ranch");
        mkdir(ranch_path);
        rmdir(ranch_path, 's');
        mkdir(ranch_path);
        
        mis = load(fullfile(folder_path, "microstate_info.mat"));
        mis = mis.mis;
        
        % Check if SSA, DSA, DSAw/L corresponding to mis(7) = 1, 2, 3
        if mis(7) == 1  % SSA
            seq = readlines(fullfile(folder_path, "sequence.seq"));
            seq = strlength(seq(2));
            ast = readlines(fullfile(folder_path, "assignment.txt"));
            ast = split(ast(2));
        
            % Check that microstate is not fully folded
            if ~(str2double(ast(2)) == 1 && str2double(ast(3)) == seq)
                % Build command
                command = 'ranch --repetitions '+reps+' --model-format=cif --prefix '+fullfile(ranch_path, sprintf('%d_%d_rank%d_', mis(1), mis(7), mis(8)))+' '+fullfile(folder_path, 'assignment.txt')+' '+fullfile(folder_path, 'sequence.seq')+' '+fullfile(folder_path, 'dom1.pdb');
                % Execute RANCH from terminal
                [status, cmdout] = system(command);
            else
                % If the microstate is fully folded, RANCH cannot process
                % it (throws an error). Therefore, just copy the input .cif
                % structure to mimic a fully folded output.
                copyfile(fullfile("..", "data", pname+".cif"), ranch_path) 
                [status, cmdout] = movefile(fullfile(ranch_path, pname+".cif"), fullfile(ranch_path, sprintf('%d_%d_rank1.cif', mis(1) ,mis(7))));
                if status == 1
                    status = 0;
                end
            end
        else  % DSA or DSAw/L; discrimination between these two is built into assignment.txt  
            command = 'ranch --repetitions '+reps+' --model-format=cif --prefix '+fullfile(ranch_path, sprintf('%d_%d_rank%d_', mis(1), mis(7), mis(8)))+' '+fullfile(folder_path, 'assignment.txt')+' '+fullfile(folder_path, 'sequence.seq')+' '+fullfile(folder_path, 'dom1.pdb')+' '+fullfile(folder_path, 'dom2.pdb');
            [status, cmdout] = system(command);
        end
    
        % Check for successful execution of RANCH
        if status == 0
            disp(sprintf('%d of %d', i-2, n-2) + " done with status " + sprintf('%d', status))
        else
            error("EnsembleWSME - Ranch exited with non-zero status %d \n %s", status, cmdout)
        end
    end
    
    disp("run_ranch done")
end