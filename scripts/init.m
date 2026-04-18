%%% Initializiation
%   Creates missing directories and/or checks for input files

disp("Beginning initialization")

thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);

% Due to difference in launching of terminal, it is required to manually
% set the ATSAS environment variable on Mac.
if ismac
    error("EnsembleWSME : Running on MacOS. Delete line 12 and replace line 13 in init.m with setenv('ATSAS' , '/path/to/atsas')")
    % setenv('ATSAS' , '/path/to/atsas')
end

if ~isfile(fullfile(pwd, "main_RANCH.m")) || ~isfile(fullfile(pwd, "main_WSME.m"))
    error("Main scripts missing. Check if working directory is set to \scripts\")
end

% Check and create /results/ directory
if ~isfolder(fullfile(pwd, "..", "results"))
    mkdir(fullfile(pwd, "..", "results"))
    mkdir(fullfile(pwd, "..", "results", "microstates"))
    mkdir(fullfile(pwd, "..", "results", "pools"))
    disp("Created 'results' directory.")
end

% If the /data/ directory is not present (first initialization), creates it
if ~isfolder(fullfile(pwd, "..", "data"))
    mkdir(fullfile(pwd, "..", "data"))
    mkdir(fullfile(pwd, "..", "data", "WSME_outputs"))
    disp("Created 'data' directory.")
% If /data/ is present (subsequent runs) verifies input files
else
    folderpath = fullfile(pwd, "..", "data");
    if isempty(dir(fullfile(folderpath, "*.pdb")))
        error("Error : No input .pdb structure found!")
    elseif isempty(dir(fullfile(folderpath, "*.cif"))) 
        error("Error : No input .cif structure found!")
    else
        disp("All input files present")
    end
end

% Modifies MATLAB path variable
addpath(genpath(fullfile(thisDir, 'pipelines')))
addpath(genpath(fullfile(thisDir, '..', 'src')))
addpath(genpath(fullfile(thisDir, '..', 'data')))
addpath(genpath(fullfile(thisDir, '..', 'results')))
disp("Set up path")

fclose('all');
disp("Initialization done")