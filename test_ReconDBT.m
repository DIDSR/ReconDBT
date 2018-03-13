%test ReconDBT

%Create DBT projection views
disp 'Create DBT projection views'
run gen_dbtproj_example;

%Test circular arc DBT reconstruction
disp 'Test circular arc DBT reconstruction'
run dbt_recon_circular_example;

%Test non-arc DBT reconstruction
disp 'Test non-circular DBT reconstruction'
run dbt_recon_noncircular_example;
