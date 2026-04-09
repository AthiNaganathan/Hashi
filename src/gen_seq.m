function [len, seq, ref] = gen_seq(filename)
    %%% Extract the protein sequence from the pdb structure
    %   Outputs the length of sequence, the sequence itself, and a mapping
    %   from atoms to amino acids.

    pdb_data = pdbread(filename);
    resName = {pdb_data.Model.Atom.resName};
    resSeq = [pdb_data.Model.Atom.resSeq];
    resSeq = resSeq - min(resSeq) + 1;
    len = max(resSeq);
    ref = zeros(1, len);

    res = cell(1, len);
    cur = 0;
    for i = 1:length(resSeq)
        if resSeq(i) ~= cur
            res(resSeq(i)) = resName(i);
            ref(resSeq(i)) = i;
            cur = cur + 1;
        end
    end

    threeToOneMap = containers.Map({'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}, ...
        {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'});
    seq = '';
    
    for i = 1:length(res)
        aa = res{i};
        seq = [seq threeToOneMap(aa)];
    end
end