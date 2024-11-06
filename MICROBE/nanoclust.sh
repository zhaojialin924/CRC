mkdir db db/taxdb
wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz && tar -xzvf 16S_ribosomal_RNA.tar.gz -C db
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -xzvf taxdb.tar.gz -C db/taxdb

db="/share/data0/UserData/beile/db/nanoclust/speciedb/16S_ribosomal_RNA"
tax="/share/data0/UserData/beile/db/nanoclust/taxdb/"
fastq="/share/data0/UserData/beile/project/fl16s/qc"
path="/share/data0/UserData/beile/software/NanoCLUST"

nextflow run $path/main.nf -profile conda --reads $fastq/"$sample".final.fastq.gz --db $db --tax $tax --min_read_length 1000 -resume

