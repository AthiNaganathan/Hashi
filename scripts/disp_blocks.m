function disp_blocks(pname, pools)
    % Display 1D macroscape landscape of the protein
    
    arguments
        pname char          % Protein name is required
        pools double = 0    % User-defined pools (optional)
    end

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

    % Plot free energy as a function of number of structured blocks
    yyaxis left
    plot(x, -8.314.*298.*log(y),'LineWidth',3)
    ylabel("Delta G (kJ/mol)")
    yticklabels({})

    % If user-defined pools are provided, highlight the same
    if ~isscalar(pools)
        for i = 1:size(pools, 1)
            y_lim = ylim;
            rectangle('Position', [pools(i, 1) 0 pools(i, 2)-pools(i, 1) y_lim(2)], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end

    % Plot probability of different macrostates
    yyaxis right
    plot(x, y, '--','LineWidth',3)
    ylabel("Probability")
    yticklabels({})

    xlabel("No. of structured blocks")
end