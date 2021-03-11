#!/usr/bin/env bash
set -euxo pipefail
/data/user/yangzz/worktools/Plink/plink --bcf ~/mapping/08.mergeGVCF/field_cultivar/chr1B.trimed.snp.bcf.gz --make-bed --allow-extra-chr --out chr1B_setname --set-missing-var-ids @:#
