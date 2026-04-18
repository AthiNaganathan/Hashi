function disp_blocks(pname, pools)
    % Display 1D macroscape landscape of the protein
    
    arguments
        pname char          % Protein name is required
        pools double = 0    % User-defined pools (optional)
    end

    T = input("Enter Temperature (in K) : ");

    filename = pname + "_pepval.mat";
    raw_data = load(filename);
    raw_data = raw_data.pepval;

    x = min(raw_data(:, 1)) : max(raw_data(:, 1));
    y = zeros(1, length(x));
    for i = x
        temp = find(raw_data(:, 1) == x(i));
        y(i) = sum(raw_data(temp, 2));
    end
    y = y ./ sum(y);
    clf('reset')
    hold on

    % Right axis (area)
    yyaxis right
    area(x, y,'LineStyle',':', 'EdgeColor','b', 'FaceColor', [0.9 0.9 0.9], 'FaceAlpha', 0.5);
    ylabel("Probability")
    set(gca, 'YColor', 'b');
    
    % Left axis (line)
    yyaxis left
    y_data = -8.314.*T.*log(y)/1000;
    plot(x, y_data,'k');
    ylabel("FE (kJ mol^{-1})")
    yticks(0:10:max(y_data))
    ylim([min(y_data) max(y_data)+0.1*range(y_data)])
    set(gca, 'YColor', 'k');

    xlabel("No. of Structured Blocks")
    xticks(0:10:max(x))

    % If user-defined pools are provided, highlight the same
    if ~isscalar(pools)
        for i = 1:size(pools, 1)
            y_lim = ylim;
            rectangle('Position', [pools(i, 1) 0 pools(i, 2)-pools(i, 1) y_lim(2)], 'FaceColor', [1 0.5 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        end
    end
end