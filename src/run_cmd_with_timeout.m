function [status, cmdout] = run_cmd_with_timeout(args, timeoutSec, outFile)
    import java.util.concurrent.TimeUnit

    if ~iscell(args)
        error('run_cmd_with_timeout:args', ...
            ['ARGS must be a cell array of tokens, e.g. ' ...
             '{''ranch'', ''--repetitions'', ''5'', path1, path2, ...}. ' ...
             'A single command string cannot be used any more because it ' ...
             'is no longer parsed by a shell (that shell was the source ' ...
             'of the orphaned-process bug).']);
    end

    argsChar = cellfun(@char, args, 'UniformOutput', false);
    pb = java.lang.ProcessBuilder(argsChar);

    if nargin < 3 || isempty(outFile)
        outFile = [tempname, '.txt'];
    end
    pb.redirectErrorStream(true);
    pb.redirectOutput(java.io.File(char(outFile)));

    proc = pb.start();
    finished = proc.waitFor(timeoutSec, TimeUnit.SECONDS);

    if ~finished
        proc.destroyForcibly();
        proc.waitFor(5, TimeUnit.SECONDS);   % give the OS a moment to reap it
        status = 124;                        % convention for timeout errors
    else
        status = proc.exitValue();
    end

    if isfile(outFile)
        cmdout = fileread(outFile);
        delete(outFile);
    else
        cmdout = '';
    end
end