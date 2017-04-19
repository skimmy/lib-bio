# Reads statistic
bio-tk readstats -r reads.fastq -d out_dir/ -X file-prefix
# Generate data accoding to configuration file gen.btconf
bio-tk generate -C gen.btconf
