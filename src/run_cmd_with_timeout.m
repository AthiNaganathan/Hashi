function [status, cmdout] = run_cmd_with_timeout(command, timeoutSec)
    import java.lang.ProcessBuilder
    import java.util.concurrent.TimeUnit

    if ispc
        pb = ProcessBuilder({'cmd.exe', '/c', char(command)});
    else
        pb = ProcessBuilder({'/bin/sh', '-c', char(command)});
    end

    % Merge stderr into stdout and send both to a temp file, so we never
    % block on a full pipe buffer while waiting.
    outFile = fullfile("..", "results", "log1.txt");
    pb.redirectErrorStream(true);
    pb.redirectOutput(java.io.File(outFile));

    proc = pb.start();
    finished = proc.waitFor(timeoutSec, TimeUnit.SECONDS);

    if ~finished
        proc.destroyForcibly();
        proc.waitFor(5, TimeUnit.SECONDS);
        status = 124; % convention for timeout errors
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