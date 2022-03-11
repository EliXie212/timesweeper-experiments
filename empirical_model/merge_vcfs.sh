proc_vcf () {
    bgzip -c $1 > $1.gz
    bcftools index $1.gz
}

for i in vcfs/*.vcf 
do
    proc_vcf $i &
done
wait

bcftools merge -Oz --threads True --force-samples -0 vcfs/S11696.Y1.E2.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I11698_green.hg19.vcf.gz vcfs/S13179.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I13698_v39.hg19.vcf.gz vcfs/S11697.Y1.E2.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S14000.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6358_aln.sort.mapped.rmdupse_adna_v2.md.vcf.gz vcfs/I7021_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6221_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S13180.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S12977.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S13957.Y1.E2.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6361_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S12978.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S12957.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S12955.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S13173.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6348_green.hg19.vcf.gz vcfs/S13964.Y1.E2.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S12960.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S12976.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S14037.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6347_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6264_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6262_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S12973.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S13766.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S13767.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6367_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S12975.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6363_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I13768_green.hg19.vcf.gz vcfs/I7039_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I7033_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S12969.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S13505.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6352_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6364_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6362_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6353_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6351_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S12971.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I7032_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S14194.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S13963.Y1.E2.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6349_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6365_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6359_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6369_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S13504.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S12970.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I7027_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S13965.Y1.E2.L1.1240k_plus.hg19.v1.vcf.gz vcfs/I6356_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I7022_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I7030_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I7024_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6232_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I7023_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6233_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6224_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I7029_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6263_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6226_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6357_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6230_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/I6228_aln.sort.mapped.rmdupse_adna_v2.md.1.vcf.gz vcfs/S13175.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz vcfs/S12974.Y1.E1.L1.1240k_plus.hg19.v1.vcf.gz > ts_merged.vcf.gz