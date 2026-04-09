function res_mis = strb2res_mis(mis, sb2r_ref)
    %%% Convert structured blocks indices to residue indices using
    %%% reference (extracted from BlockDet)
    
    res_mis = zeros(size(mis));
    res_mis(:, [1 2]) = mis(: , [1 2]);
    res_mis(:, 7:end) = mis(: , 7:end);

    for i = 1:size(mis, 1)
        mi = mis(i, :);
        if mi(3) ~= 1
            res_mis(i, 3) = sum(arrayfun(@(x) length(cell2mat(x)), sb2r_ref(1:mi(3)-1))) + 1;
        else
            res_mis(i, 3) = 1;
        end
        res_mis(i, 4) = sum(arrayfun(@(x) length(cell2mat(x)), sb2r_ref(mi(3):mi(3)+mi(4)-1)));
        if mi(5) ~= 0
            res_mis(i, 5) = sum(arrayfun(@(x) length(cell2mat(x)), sb2r_ref(1:mi(5)-1))) + 1;
            res_mis(i, 6) = sum(arrayfun(@(x) length(cell2mat(x)), sb2r_ref(mi(5):mi(5)+mi(6)-1)));
        end
    end
end