%%% Run RANCH to generate 3D structural ensembles

init

% Run disp_blocks(pname) in console to visualize graph to assign splits

%%% Setup Input Parameters

% % Example 1 : Villin
% pname = 'Villin'; % Protein name
% macrostate_pools = [25 34]; % Split 1D landscape into pools wherein to visualize ensemble
% no_of_microstates = 20; % No. of microstates considered per pool of macrostates
% no_of_conformations = 5; % No. of conformations to be generated per microstate

% Example 2 : Pertactin
pname = 'Pertactin'; % Protein name
macrostate_pools = [68 76; 107 111; 128 137]; % Split 1D landscape into pools wherein to visualize ensemble
no_of_microstates = 5; % No. of microstates considered per pool of macrostates
no_of_conformations = 5; % No. of conformations to be generated per microstate


% Run disp_blocks(pname, macrostate_pools) in console to visualize assigned splits

% Confirm intent to run
confirm = input('WARNING : About to reset contents of \\results\\. Do you want to proceed? (y/[n]): ', 's');
if ~strcmpi(confirm, 'Y')
    error('Operation cancelled by user.');
end

%%% Run RANCH

% Generate the required input files to call RANCH
% Refer RANCH documentation for full details
generate_ranch_files(pname, no_of_microstates, macrostate_pools);

tic
% Executes RANCH via terminal
run_ranch(pname, no_of_conformations)
toc

% Group RANCH outputs into user-defined pools
make_pools(macrostate_pools, no_of_microstates);

disp("Done.")

