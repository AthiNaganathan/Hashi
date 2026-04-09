function generate_dom_pdbs(mi_dom, pdb_data, ref, path)
   %%% Generate the structural information for individual domains by
   %%% splicing it out of the input pdb structure.

   %   Does not return any values, directly saves the domains in the
   %   required folder.

   atoms = pdb_data.Model.Atom;
   mi = mi_dom(1, :);

   ref = [ref size(atoms, 2)+1];
   
   mi(3) = ref(mi_dom(3));
   mi(4) = ref(mi_dom(3)+mi_dom(4))-ref(mi_dom(3));
   if mi_dom(5) ~= 0
       mi(5) = ref(mi_dom(5));
       mi(6) = ref(mi_dom(5)+mi_dom(6))-ref(mi_dom(5));
   end

   dom1 = atoms(mi(3):mi(3)+mi(4)-1);
   dom1_ser = [dom1.AtomSerNo];
   dom1_ser = dom1_ser - min(dom1_ser) + 1;

   for i = 1:size(dom1, 2)
       dom1(1, i).AtomSerNo = dom1_ser(1, i);
   end

   new_dom1 = pdb_data;
   new_dom1.Model.Atom = dom1;
   new_dom1.Model.Terminal = {};
   pdbwrite(fullfile(path, "dom1.pdb"), new_dom1);

   if mi(5) ~= 0 % Indirect check for if DSA or DSAw/L
       dom2 = atoms(mi(5):mi(6)+mi(5)-1);
       dom2_ser = [dom2.AtomSerNo];
       dom2_ser = dom2_ser - min(dom2_ser) + 1;
    
       for i = 1:size(dom2, 2)
           dom2(1, i).AtomSerNo = dom2_ser(1, i);
       end
    
       new_dom2 = pdb_data;
       new_dom2.Model.Atom = dom2;
       new_dom2.Model.Terminal = {};
       pdbwrite(fullfile(path, "dom2.pdb"), new_dom2);
   end
end