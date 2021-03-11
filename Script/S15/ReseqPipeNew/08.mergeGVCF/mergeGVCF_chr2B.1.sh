#!/usr/bin/env sh
java  \
  -Djava.io.tmpdir=/home/wangzh/tmp/ \
  -jar /home/wangzh/bin/GenomeAnalysisTK.jar \
  -R /home/wangzh/Project/reseq/indexSplit/IWGSC_MP/WholeGenomeSplit_MP.fa \
  -T GenotypeGVCFs \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164642/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164657/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164652/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/LM5-8/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164686/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164681/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164661/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164639/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164680/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164682/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164643/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164683/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR7164689/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S10/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S64/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S67/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S135/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S136/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S137/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S142/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S185/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S186/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S188/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S211/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S212/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S213/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S214/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S215/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S216/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S217/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S218/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S223/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S225/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S227/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/S228/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/C36/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/C38/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/C43/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/YM02/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/YM03/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/YM05/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/YM07/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/YM08/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/NZ02/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/NZ05/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/NZ10/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/NZ11/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/NZ38/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH07/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH08/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH10/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/HZ01/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/HZ02/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/HZ03/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/HZ04/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/HZ05/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/HZ06/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH42/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH43/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH46/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH47/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH48/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH49/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH50/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH51/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH52/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH53/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH54/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH55/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH56/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH57/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH58/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH59/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH60/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH61/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH62/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH63/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH64/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH65/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH66/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH67/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH68/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH69/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH70/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH71/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH72/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH73/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH74/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH75/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH76/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH77/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH78/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH79/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH80/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH81/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH82/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH83/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH84/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH85/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH86/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH87/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH88/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH89/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH90/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH91/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH92/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH93/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH94/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH95/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH96/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH97/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH98/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH99/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH100/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH101/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH102/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH103/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH104/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH105/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH106/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH107/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH108/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH109/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH110/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH111/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH112/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH113/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH114/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH115/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH116/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH117/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH118/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH119/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH120/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH121/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH122/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH123/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH124/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH125/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH126/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH127/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH128/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH129/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH130/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH131/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH132/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH133/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH134/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH135/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH136/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH137/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH138/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH139/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH140/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH141/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH142/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH143/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH144/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH145/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH146/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH147/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH148/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH149/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH150/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH151/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH152/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/PH153/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/Ak58Simu/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766513/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766527/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766532/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766565/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766609/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766620/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766623/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766631/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR512/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766546/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766626/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766574/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766601/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766605/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766610/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766611/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766493/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766496/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766525/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766562/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766570/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766503/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766530/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766612/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766613/06.callsnp/chr2B.1.g.vcf.gz \
  --variant /data/user/yangzz/mapping/Reseq_data/SRR10766614/06.callsnp/chr2B.1.g.vcf.gz \
  -o chr2B.1.raw.vcf.gz
