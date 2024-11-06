picrust2_pipeline.py -s temp/pircust2/{1}/{1}.final.fa -i temp/pircust2/{1}/final.tsv -o temp/pircust2/{1}/final -p 20 --stratified
pathway="/share/data0/UserData/beile/miniconda3/envs/picrust2/lib/python3.8/site-packages/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv"
pathway_pipeline.py -i temp/pircust2/$sample/results/KO_metagenome_out/pred_metagenome_unstrat.tsv -o temp/pircust2/$sample/results/KEGG_pathways_out --no_regroup --map $pathway

