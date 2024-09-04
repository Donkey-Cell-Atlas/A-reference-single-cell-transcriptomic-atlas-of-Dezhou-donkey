## DNBC4tools 2.1.0 (https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software.git)
### Build index for reference genome
DNBC4tools mkref --action mkgtf --ingtf Equus_asinus.ASM1607732v2.109.gtf --outgtf gene.filter.gtf \
----attribute gene_biotype:protein_coding \
                        gene_biotype:lncRNA \
                        gene_biotype:IG_C_gene \
                        gene_biotype:IG_D_gene \
                        gene_biotype:IG_J_gene \
                        gene_biotype:IG_LV_gene \
                        gene_biotype:IG_V_gene \
                        gene_biotype:IG_V_pseudogene \
                        gene_biotype:IG_J_pseudogene \
                        gene_biotype:IG_C_pseudogene \
                        gene_biotype:TR_C_gene \
                        gene_biotype:TR_D_gene \
                        gene_biotype:TR_J_gene \
                        gene_biotype:TR_V_gene \
                        gene_biotype:TR_V_pseudogene \
                        gene_biotype:TR_J_pseudogene
						
DNBC4tools mkref --action mkref --ingtf gene.filter.gtf \
            --fasta Equus_asinus.ASM1607732v2.dna.toplevel.fa \
            --genomeDir . \
            --thread 40
			
### Running the main workflow
DNBC4tools run \
	--cDNAfastq1 /test/data/test_R1.fastq.gz \
	--cDNAfastq2 /test/data/test_R2.fastq.gz \
	--oligofastq1 /test/data/test_oligo_R1.fq.gz \
	--oligofastq2 /test/data/test_oligo_R2.fq.gz \
	--genomeDir /genome/genomeDir \
	--gtf /genome/gene.filter.gtf \
	--name sntestis \
	--species Equus_asinus \
	--thread 40