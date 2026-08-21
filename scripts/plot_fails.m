function plot_fails(skip_mis, seq_len)
    clf("reset")

    vis = zeros(size(skip_mis, 1), seq_len);
    for i = 1:size(skip_mis, 1)
        row = zeros(1, seq_len);
        cur = skip_mis(i, 3:6);
    
        row(cur(1):cur(1)+cur(2)-1) = repelem(1, cur(2));
        row(cur(3):cur(3)+cur(4)-1) = repelem(1, cur(4));
        vis(i, :) = row;
    end

    n = ceil(sqrt(size(skip_mis, 1)));
    tiledlayout(n, n)

    for i = 1:size(skip_mis, 1)
        nexttile
        plot(vis(i, :))
        xlabel("Residue index")
        xlim([0 seq_len])
        ylabel("Structure state")
        ylim([0 1])
        title(sprintf('Pool %d | A : %d', skip_mis(i, 8), skip_mis(i, 7)))
    end

    disp("Failures plotted")
end