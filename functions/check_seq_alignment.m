function check_seq_alignment(sp_base, varargin)
%CHECK_SEQ_ALIGNMENT  Fail loudly if the map layers disagree on their species key.
%
%   check_seq_alignment(sp_base, "data/ebirdatlas.mat", "data/kbmatlas.mat")
%
% The three species x grid layers (map_old, map_ebird, map_kbm) are combined
% by array position - map_kbm | map_ebird, and sp_base(i) naming slice i - so
% all three must be indexed by the same SEQ vector. Re-running
% A_import_old_atlas without also re-running B and C leaves them silently
% misaligned, and every species then carries some other species' data.
%
% That is not hypothetical: the 2026-09 taxonomy update rebuilt oldatlas.mat
% and ebirdatlas.mat (SEQ-sorted) but left kbmatlas.mat from 2024 (unsorted),
% so every species on the exported website carried an unrelated species' KBM
% records. Nothing detected it. Hence this check.

for i = 1:numel(varargin)
    f = varargin{i};
    d = load(f, "seq_map");
    assert(isfield(d, "seq_map"), ...
        "%s has no seq_map: it predates the alignment check. Re-run A, B and C together.", f);
    assert(isequal(d.seq_map, sp_base.SEQ), ...
        "%s was built from a different species list than data/oldatlas.mat. Re-run A, B and C together.", f);
end
end
