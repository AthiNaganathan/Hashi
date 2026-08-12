%%% Run the bWSME model
%   Run the bWSME model to obtain 1-D microstate probabilities

if ~exist("init.m", "file")
    error("init function not found. Set pwd to /scripts/")
end

init

% Set protein name

% Villin
pname  = 'Villin';
BlockSize = 1;
is_long_protein = 0;
struct_file = 'struct.txt';
pH = 7;
srcutoff = 5.0;
ene = -87.5/1000;
DS = -14.5/1000;
DCp = -0.36/1000;
T = 310;
IS = 0.15;

% % Pertactin
% pname = 'Pertactin';
% BlockSize = 5;
% is_long_protein = 1;
% struct_file = 'struct.txt';
% pH = 7;
% srcutoff = 5.0;
% ene = -79/1000;
% DS = -14.5/1000;
% DCp = -0.36/1000;
% T = 298;
% IS = 0.05;

%%% Some verifications

% Check if struct.txt is present
folderpath = fullfile(pwd, "..", "data");
if isempty(dir(fullfile(folderpath, "*.txt")))
    error("Error : No stride output file (.txt) found!")
end
% Empty WSME_outputs from previous run
confirm = input('WARNING : About to reset contents of \\data\\WSME_outputs\\. Do you want to proceed? (y/[n]): ', 's');
if ~strcmpi(confirm, 'Y')
    error('Operation cancelled by user.');
end
rmdir(fullfile(pwd, "..", "data", "WSME_outputs"),'s');
mkdir(fullfile(pwd, "..", "data", "WSME_outputs"))

%%% Refer to individual function scripts within /scripts/pipelines/ for
%%% complete documentation of input parameters

cmapCalcElecBlock(pname, struct_file, pH, srcutoff, BlockSize)
tic
if ~is_long_protein
    disp("Not a long protein! Running FesCalc_Block...")
    FesCalc_Block(pname, struct_file, ene, DS, DCp, T, IS)
else
    disp("Long protein! Running FesCalc_Block_gen...")
    FesCalc_Block_gen(pname, struct_file, ene, DS, DCp, T, IS)
end
toc


% Move files used specifically for RANCH
copyfile(fullfile(pwd, '..', 'data', 'WSME_outputs', [pname '_pepval.mat']), fullfile(pwd, '..', 'data', [pname '_pepval.mat']))
copyfile(fullfile(pwd, '..', 'data', 'WSME_outputs', [pname '_BlockDet.dat']), fullfile(pwd, '..', 'data', [pname '_BlockDet.mat']))

disp("WSME Done")