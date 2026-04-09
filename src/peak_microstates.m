function fin_mis = peak_microstates(filename, no_of_microstates, pools)
    %%% Return the top microstates from each pool of macrostates
    raw_data = load(filename);
    raw_data = raw_data.pepval;
    
    fin_mis = [];

    if ~isscalar(pools)
        for i = 1:size(pools, 1)
            reps_ens = raw_data(find(raw_data(:, 1) >= pools(i, 1) & raw_data(:, 1) <= pools(i, 2)), :);
            reps_ens = sortrows(reps_ens, 2, 'descend');
            
            % Select the top microstates
            fin_mis = [fin_mis; reps_ens(1:no_of_microstates, :)];
        end
        [~, index] = sortrows(fin_mis, 2, 'descend');
        fin_mis = [fin_mis index];
    else
        error("Pools not defined properly")
    end
end