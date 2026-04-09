%%% Run the bWSME model
%   Run the bWSME model to obtain 1-D microstate probabilities

init

% Set protein name

% Villin
pname  = 'Villin';

% % Pertactin
% pname = 'Pertactin';

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

% Inputs = protein name, stride output file, pH, srcutoff, BlockSize
cmapCalcElecBlock(pname, 'struct.txt', 7, 5.0, 1)
% Inputs = protein name, stride output file, ene, DS, DCp, T, IS
FesCalc_Block(pname, 'struct.txt', -98/1000, -14.5/1000, -0.3579/1000, 310, 0.1)

% Move files used specifically for RANCH
copyfile(fullfile(pwd, '..', 'data', 'WSME_outputs', [pname '_pepval.mat']), fullfile(pwd, '..', 'data', [pname '_pepval.mat']))
copyfile(fullfile(pwd, '..', 'data', 'WSME_outputs', [pname '_BlockDet.dat']), fullfile(pwd, '..', 'data', [pname '_BlockDet.mat']))

disp("WSME Done")