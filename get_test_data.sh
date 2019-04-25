#!/bin/bash

set -eou pipefail

# CRAM FILES ---------------------------------------------------------------------------
# wget https://storage.googleapis.com/gatk-test-data/wgs_cram/NA12878_20k_hg38/NA12878.cram
# wget https://storage.googleapis.com/gatk-test-data/wgs_cram/NA12878_20k_hg38/NA12878.crai

# BAM FILES ---------------------------------------------------------------------------
wget https://storage.googleapis.com/gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_med.hg38.bam
wget https://storage.googleapis.com/gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_med.hg38.bai

# VCF FILES ---------------------------------------------------------------------------
wget https://storage.googleapis.com/gatk-test-data/wgs_vcf/PlatinumTrio_b37/PlatinumTrio_b37.vcf.gz
wget https://storage.googleapis.com/gatk-test-data/wgs_vcf/PlatinumTrio_b37/PlatinumTrio_b37.vcf.gz.tbi

# bcftools view -Ob PlatinumTrio_b37.vcf.gz > PlatinumTrio_b37.bcf
# bcftools index PlatinumTrio_b37.bcf
