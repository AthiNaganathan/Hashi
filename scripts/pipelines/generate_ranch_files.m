function generate_ranch_files(pname, no_of_microstates, pools) 
    % Generate all input files required to run RANCH

    % Read input .pdb structure
    pdb = pdbread(pname + ".pdb");
    atoms = pdb.Model.Atom;
    chainIDs = {atoms.chainID};    
    chainID = unique(chainIDs);

    % Check that .pdb has only one chain
    if length(chainID) > 1
        error("ERROR : RANCH cannot work with multiple-chain .pdb files")
    end

    % Extract and set chain name from the structure
    sub_name = chainID{1};
    fprintf("Assigning Chain ID %s\n", string(sub_name))

    % Protein sequence | A2AA_ref = atoms to amino acids reference
    [len, seq, A2AA_ref] = gen_seq(pname + ".pdb");

    % MIS_strb = microstates to structural blocks;
    MIS_strb = peak_microstates(pname + "_pepval.mat", no_of_microstates, pools);

    % structured blocks to residues (sb2r) reference
    SB2R_REF = strb2res_ref(pname + "_BlockDet.dat");

    % convert strb mis to residue (res) mis
    MIS_res = strb2res_mis(MIS_strb, SB2R_REF);

    % generate assignment texts
    ASSMTS = res2assmts(MIS_res, len, sub_name);

    % native pdb structure
    pdb_data = pdbread(pname + ".pdb");

    % Generate all ranch-ready files
    foldername = fullfile(pwd, "..", "results", "microstates");
    mkdir(foldername);
    rmdir(foldername, 's');
    mkdir(foldername);

    for i = 1:length(ASSMTS)
        % Read microstate information and create folder for RANCH files
        mis = MIS_res(i, :);
        path = fullfile(foldername, sprintf('%05d_%d_%d_%0.4f',i, mis(1), mis(7), mis(8)));
        mkdir(path);
        mkdir(fullfile(path, "ranch"));

        % Write the protein sequence
        fid = fopen(fullfile(path, 'sequence.seq'), 'w');
        fprintf(fid, '> %s\n%s\n', sub_name, seq);
        fclose(fid);

        % Write the assignment of ordered/disordered regions
        fid = fopen(fullfile(path, 'assignment.txt'), 'w');
        fprintf(fid, '%s\n', ASSMTS(i));
        fclose(fid);

        % Write the domain structural information
        generate_dom_pdbs(mis, pdb_data, A2AA_ref, path);

        % Save the microstate information for future reference
        save(fullfile(path, 'microstate_info.mat'), "mis")
    end
    
    disp("generate_ranch_files done")
end