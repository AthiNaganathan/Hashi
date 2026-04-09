function strb_res_ref = strb2res_ref(filename)
    %%% Extract the structured blocks to residue reference from BlockDet
    
    data = load(filename);
    strb_res_ref = cell(1, max(data(:, 2)));
    for i = 1:size(data, 1)
        entry = data(i, :);
        strb_res_ref(1, entry(2)) = {[cell2mat(strb_res_ref(1, entry(2))) entry(1)]};
    end
end
