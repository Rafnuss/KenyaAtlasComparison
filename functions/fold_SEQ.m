function seq = fold_SEQ(seq)
%FOLD_SEQ Join legacy external SEQs through the curated historical mapping.
opts_mapping = detectImportOptions('data/taxonomy/atlas_mapping.csv', TextType='string');
opts_mapping = setvartype(opts_mapping, 'atlas_seq', 'string');
mapping = readtable('data/taxonomy/atlas_mapping.csv', opts_mapping);
[found, idx] = ismember(string(seq), mapping.atlas_seq);
keep = false(size(seq));
keep(found) = mapping.external_use(idx(found)) == "true";
seq(~keep) = 0;
seq(keep) = mapping.comparison_seq(idx(keep));
end
