grep "gi" ant_warner.map | sed 's/.*ref\|//' | sed 's/\|.*//' | sort | uniq > placed_scaff
grep ">" GCF_000980195.1_M.pharaonis_V2.0_genomic.fna | sed 's/>//' | sed 's/ .*//' | sort > all_scaff
sed 's/ .*//' GCF_000980195.1_M.pharaonis_V2.0_genomic.fna > Mphar.fna
comm -23 all_scaff placed_scaff > unplaced_scaff
xargs samtools faidx Mphar.fna < unplaced_scaff > unplaced.fa
blastn -db ant_warner.fasta -query unplaced.fa -out results.out -outfmt 6 -max_target_seqs 1
cat results.out | awk '{{print$1}}' | uniq > suc_blast
comm -23 unplaced_scaff suc_blast > still_missing
grep "exon" GCF_000980195.1_M.pharaonis_V2.0_genomic.gff | gff2bed > exon.bed
grep "gene=" exon.bed | sed 's/\([A-Z_a-z\.0-9]*\).*gene=\([^;]*\).*/\1 \2/' | sort | uniq > all_genes
