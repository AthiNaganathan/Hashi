%%% Run RANCH to generate 3D structural ensembles

if ~exist("init.m", "file")
    error("init function not found. Set pwd to /scripts/")
end

init

% Run disp_blocks(pname) in console to visualize graph to assign splits

%%% Setup Input Parameters

timeoutSec = 300; % How long to run RANCH for (longer explanation towards end of code, and in documentation)

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
run_ranch(pname, no_of_conformations, no_of_microstates, timeoutSec)
toc

% Group RANCH outputs into user-defined pools
make_pools(macrostate_pools, no_of_microstates);

% Display skipped microstates
load(fullfile("..", "results", "skipped_mis.mat"))
if size(skip_mis, 1)
    fprintf('\n\n--- WARNING ---\n\n%d microstates were timed out.\n', size(skip_mis, 1))
    for i = unique(skip_mis(:, end)).'
        fprintf('%d microstate(s) from pool_%d were skipped\n', length(find(skip_mis(:, end) == i)), i);
    end

    [seq_len, ~, ~] = gen_seq(pname + ".pdb");
    
    % Plot the microstates that fail. Refer documentation for complete
    % details; there appears to be an edge case where if RANCH is asked to
    % model a very small disordered region with a particular sequence
    % between two structured regions in the DSAw/L approximation, the
    % algorithm gets stuck in a loop and does not quit out even after
    % arbitrarily long times. 
    % 
    % To avoid this, the program automatically kills the RANCH process 
    % after timeoutSec number of seconds, skips that microstate, stores it,
    % and displays it after the run is complete. Resolving this behaviour 
    % completely would require deeper understanding of the implementation 
    % of the RANCH algorithms.
    %
    % Thankfully, this behaviour is extremely anamalous; in our testing,
    % only one protein out of >30 tested showed the algorithm getting stuck
    % in such loops, and that too for very specific microstates.

    plot_fails(skip_mis, seq_len)
else
    disp("All microstates were successful.")
end

disp("Done.")

