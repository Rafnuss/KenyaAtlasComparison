[status, message] = system('python3 data/taxonomy/update_taxonomy.py');
assert(status == 0, message);
A_import_old_atlas
B_import_KBM
C_import_ebird
D_correction
E_export_website
close all
