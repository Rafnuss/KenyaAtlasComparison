function check_seq_alignment(sp_base, varargin)
%CHECK_SEQ_ALIGNMENT Verify the taxon IDs indexing each map layer.
for i = 1:numel(varargin)
    d = load(varargin{i}, "avibase_map");
    assert(isequal(d.avibase_map, sp_base.avibase_id), ...
        "%s uses a different taxon list. Re-run A, B and C together.", varargin{i});
end
end
