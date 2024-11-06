module load conda
conda activate lefse

lefse_format_input.py micro_all_signature_group_lefse_count.txt micro_all_signature_group_lefse_count.in -c 2 -s -1 -u 1
lefse_run.py micro_all_signature_group_lefse_count.in micro_all_signature_group_lefse_count.res 

lefse_plot_res.py micro_all_signature_group_lefse_count.res  micro_all_signature_group_lefse_count.pdf --format pdf
lefse_plot_cladogram.py data/micro_all_signature_group_lefse_count.res fig/micro_all_signature_group_lefse_count.pdf --format pdf --title "" 