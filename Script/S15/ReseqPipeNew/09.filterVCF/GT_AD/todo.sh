#11:954072 06:jimai44 08:jinan17
#11:954072 07:ji200040919 09:jimai229
#./bisample_compare_pipline.sh -c CP07 -s CP09 -d CP07_CP09_CP11 -f 8 -n 11 &
#./bisample_compare_pipline.sh -c CP07 -s CP11 -d CP07_CP09_CP11 -f 8 -n 5 &
#./bisample_compare_pipline.sh -c CP09 -s CP11 -d CP07_CP09_CP11 -f 11 -n 5 &
#lis=(a{_b,_c})
#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${1}_${2}" -s ".${1}_${2}unmatch_homo_snp_level"  > plotfile_${1}_${2}.pt
#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${1}_${3}" -s ".${1}_${3}unmatch_homo_snp_level"  > plotfile_${1}_${3}.pt
#./ptf_maker.py -p "/data/user/yangzz/mapping/09.filterVCF/GT_AD/tmp/${2}_${3}" -s ".${2}_${3}unmatch_homo_snp_level"  > plotfile_${2}_${3}.pt

#Shi4185-1 : S64 : ZY-64 | Fu4185-1 : NZ18 : ZY-65
#Shi4185-2 : NZ33 : Shi4185 | Fu4185-2 : NZ34 : Fu4185
#./bisample_compare_pipline.sh -c S64 -v Shi4185-1 -s NZ18 -m Fu4185-1 -d fu4185_shi4185 -f 11 -n 14  #Shi4185-1 and Fu4185-1
#./bisample_compare_pipline.sh -c S64 -v Shi4185-1 -s NZ33 -m Shi4185-2 -d fu4185_shi4185 -f 11 -n 8
#./bisample_compare_pipline.sh -c S64 -v Shi4185-1 -s NZ34 -m Fu4185-2 -d fu4185_shi4185 -f 11 -n 5
#./bisample_compare_pipline.sh -c NZ18 -v Fu4185-1 -s NZ33 -m Shi4185-2 -d fu4185_shi4185 -f 14 -n 8
#./bisample_compare_pipline.sh -c NZ18 -v Fu4185-1 -s NZ34 -m Fu4185-2 -d fu4185_shi4185 -f 14 -n 5
#./bisample_compare_pipline.sh -c NZ33 -v Shi4185-2 -s NZ34 -m Fu4185-2 -d fu4185_shi4185 -f 8 -n 5

#C4-2 NZ10 5 |  LX99-1-1 YM05 8 | S64 Shi1 11 | NZ18 Fu1 14 | jimai229 CP09 17 | jimai44 CP06 20 | jinan17 CP08 23
#./bisample_compare_pipline.sh -c S64 -v Shi4185-1 -s CP08 -m JiNan17 -d S64_NZ18_CP06_CP08_CP09_YM05_NZ10 -f 11 -n 23
#./bisample_compare_pipline.sh -c NZ18 -v Fu4185-1 -s NZ10 -m JiMai22 -d S64_NZ18_CP06_CP08_CP09_YM05_NZ10 -f 14 -n 5
#./bisample_compare_pipline.sh -c NZ18 -v Fu4185-1 -s YM05 -m LiangXing99 -d S64_NZ18_CP06_CP08_CP09_YM05_NZ10 -f 14 -n 8
#./bisample_compare_pipline.sh -c NZ18 -v Fu4185-1 -s CP09 -m JiMai229 -d S64_NZ18_CP06_CP08_CP09_YM05_NZ10 -f 14 -n 17
#./bisample_compare_pipline.sh -c NZ18 -v Fu4185-1 -s CP06 -m JiMai44 -d S64_NZ18_CP06_CP08_CP09_YM05_NZ10 -f 14 -n 20
#./bisample_compare_pipline.sh -c NZ18 -v Fu4185-1 -s CP08 -m JiNan17 -d S64_NZ18_CP06_CP08_CP09_YM05_NZ10 -f 14 -n 23

#while read line;do
#echo ${line} | awk '{print $1}'
#./bisample_compare_pipline.sh -c PH09 -v Fielder -f 5 -a Fielder -s `echo ${line} | awk '{print $1}'` -m `echo ${line} | awk '{print $3}'` -b `echo ${line} | awk '{print $2}'` -d field_cultivar -n 8 &
#name1=`echo ${line} | awk '{print $1}'`
#name2=`echo ${line} | awk '{print $3}'`

#Rscript ~/R/graph_CNV.R --sample1 PH09 --sample1_name Fielder --sample2 $name1 --sample2_name $name2  -s "`ls tmp | grep $name2 | cut -d _ -f3,4`" &
#wait_all
#done < /data/user/yangzz/mapping/09.filterVCF/GT_AD/field_cultivar/cultivar_list

#./bisample_compare_pipline.sh -c PH09 -v Fielder -f 5 -a Fielder -s C13 -m VA_America3 -b 20-40982 -d field _cultivar -n 8 &

#CP01:S742 | CP03:FuGou583 | C46:AK58
#./bisample_compare_pipline.sh -c CP01 -v S742 -s CP03 -m FuGou583 -d CP01_CP03_C46 -f 11 -n 8 &
#./bisample_compare_pipline.sh -c CP01 -v S742 -s C46 -m AK58 -d CP01_CP03_C46 -f 11 -n 5 &
#./bisample_compare_pipline.sh -c CP03 -v FuGou583 -s C46 -m AK58 -d CP01_CP03_C46 -f 8 -n 5 &

#S67:nongda3338:8 | S6554:23 |  S3331:nongda3331:20 | Zang1817:z1817:14 | YM02:nongda981:5 | YM03:nongda3097:17 | S7:alegold:11
#./bisample_compare_pipline.sh -c S67 -v ND3338 -s S6554 -m S6554 -d S67_S6554_S3331_Zang1817_YM02_YM03_S7 -f 8 -n 23 &
#./bisample_compare_pipline.sh -c S67 -v ND3338 -s S7 -m Alegold -d S67_S6554_S3331_Zang1817_YM02_YM03_S7 -f 8 -n 11 &
#./bisample_compare_pipline.sh -c S3331 -v ND3331 -s Zang1817 -m z1817 -d S67_S6554_S3331_Zang1817_YM02_YM03_S7 -f 20 -n 14 &
#./bisample_compare_pipline.sh -c YM02 -v ND981 -s YM03 -m ND3097 -d S67_S6554_S3331_Zang1817_YM02_YM03_S7 -f 5 -n 17 &

#Jing411:NZ11:8 | YuMai8679:NZ12:11 | JiMai20:PH10:14 | JiMai22:NZ10:5
#./bisample_compare_pipline.sh -c NZ11 -v Jing411 -s NZ12 -m YuMai8679 -d NZ11_NZ12_PH10_NZ10 -f 8 -n 11 &
#./bisample_compare_pipline.sh -c PH10 -v JiMai20 -s NZ10 -m JiMai22 -d NZ11_NZ12_PH10_NZ10 -f 14 -n 5 &
#sample_list=S64.11.Shi4185-1,NZ18.14.Fu4185-1,NZ10.5.JiMai22,YM05.8.LiangXing99,CP09.17.JiMai229,CP06.20.JiMai44,CP08.23.JiNan17
#for i in `echo ${sample_list} | tr "," " "`; do
#    spec_sample=`echo ${i} | tr "." " "`
#    my_array=(${spec_sample})
    #for num in {1..7}; do
    #    for a in {A,B,D}; do
    #        ./muti_func_snp_compare.py --density_graph_raw_vcf_one on -i S64_NZ18_CP06_CP08_CP09_YM05_NZ10/chr${num}${a}_snp.gt -1 ${my_array[1]} --chromosome tmp/${my_array[0]}/chr${num}${a}.mask_CNV --DP_low 3 --DP_high 20 --GQ_sample 8 -b 1000000 --LEVEL "1.5" -o tmp/${my_array[0]}/chr${num}${a}.alt_homosnp_density &
    #    done
    #done
    #wait
#    cat tmp/${my_array[0]}/chr*.alt_homosnp_density | grep -v 'CHR'|grep -v 'CNV'> tmp/${my_array[0]}/combine.${my_array[0]}_alt_homosnp_density
#    Rscript ~/R/snp_density.R --sample ${my_array[0]} --sample_name ${my_array[2]} --sample_number 1
#done

# yannong15 S136 jimai22 C4-2 xiaoyanmai Elytrigia_repens
#./bisample_compare_pipline.sh -c NZ10 -v JiMai22  -a C4-2 -f 5 -s S136 -m YanNong15 -b S136 -d yn15_jm22_xym -n 8 &
#./bisample_compare_pipline.sh -c NZ10 -v JiMai22  -a C4-2 -f 5 -s XiaoYanMai -m XiaoYanMai -b Elytrigia_repens -d yn15_jm22_xym -n 8 &
#./bisample_compare_pipline.sh -c S136 -v YanNong15  -a S136 -f 5 -s XiaoYanMai -m XiaoYanMai -b Elytrigia_repens -d yn15_jm22_xym -n 8 &

#jimai44 jimai44 CP06 | jimai22 C4-2 NZ10 | jimai20 JM20 PH10
#./bisample_compare_pipline.sh -c NZ10 -v JiMai22  -a C4-2 -f 5 -s CP06 -m JiMai44 -b jimai44 -d nz10_cp06_ph10 -n 8 &
#./bisample_compare_pipline.sh -c NZ10 -v JiMai22  -a C4-2 -f 5 -s PH10 -m JiMai20 -b JM20 -d nz10_cp06_ph10 -n 8 &
#nongda5181 S140 S140 | nongda3097 s3097-1 YM03 | lunxuan987 LX987 NZ05
#./bisample_compare_pipline.sh -c S140 -v NongDa5181  -a S140 -f 5 -s YM03 -m NongDa3097 -b s3097-1 -d field_cultivar -n 8 &
#./bisample_compare_pipline.sh -c S140 -v NongDa5181  -a S140 -f 5 -s NZ05 -m LX987 -b LX987 -d field_cultivar -n 8 &

#BS_SDNingYang5 ZY-15 S15 | BS_HNPuYang21 ZY-60 S60 | BS_HNPuYang23 ZY-61 S61 | BS_SDChangQing49 ZY-62 S62
#LiangXing99  LX99-1-1 YM05 | JiMai22 C4-2 NZ10
#./bisample_compare_pipline_old.sh -c S15 -v BS_SDNingYang5   -a ZY-15 -f 5 -s NZ10 -m JiMai22 -b C4-2 -d BS_wheat -n 8 &
#wait_all
#./bisample_compare_pipline_old.sh -c S60 -v BS_HNPuYang21    -a ZY-60 -f 5 -s NZ10 -m JiMai22 -b C4-2 -d BS_wheat -n 8 &
#wait_all
#./bisample_compare_pipline_old.sh -c S61 -v BS_HNPuYang23    -a ZY-61 -f 5 -s NZ10 -m JiMai22 -b C4-2 -d BS_wheat -n 8 &
#wait_all
#./bisample_compare_pipline_old.sh -c S62 -v BS_SDChangQing49 -a ZY-62 -f 5 -s NZ10 -m JiMai22 -b C4-2 -d BS_wheat -n 8 &
#wait_all
#./bisample_compare_pipline_old.sh -c S15 -v BS_SDNingYang5   -a ZY-15 -f 5 -s YM05 -m LiangXing99 -b LX99-1-1 -d BS_wheat -n 8 &
#wait_all
#./bisample_compare_pipline_old.sh -c S60 -v BS_HNPuYang21    -a ZY-60 -f 5 -s YM05 -m LiangXing99 -b LX99-1-1 -d BS_wheat -n 8 &
#wait_all
#./bisample_compare_pipline_old.sh -c S61 -v BS_HNPuYang23    -a ZY-61 -f 5 -s YM05 -m LiangXing99 -b LX99-1-1 -d BS_wheat -n 8 &
#wait_all
#./bisample_compare_pipline_old.sh -c S62 -v BS_SDChangQing49 -a ZY-62 -f 5 -s YM05 -m LiangXing99 -b LX99-1-1 -d BS_wheat -n 8 &
#PH07 TAM107 | NZ11 C5-1
#./bisample_compare_pipline_old.sh -c PH07 -v TAM107 -a TAM107 -f 5 -s NZ11 -m Jing411 -b C5-1 -d GuanPF -n 8 &
#PH07 TAM107 |PH08 Luzi238
#./bisample_compare_pipline_old.sh -c PH07 -v TAM107 -a TAM107 -f 5 -s PH08 -m Luzi238 -b LZ238 -d GuanPF -n 8 &
#XC01 291 | XC02 291S
#./bisample_compare_pipline_new.sh -c XC01 -v 291 -a 291 -f 5 -s XC02 -m 291S -b 291S -d ZhangMH -n 8 &

#PH16 A21 PH17 A168
#./bisample_compare_pipline_new.sh -c PH16 -v A21 -a A21 -f 8 -s PH17 -m A168 -b A168 -d wangyf_PH -n 5 &
#S67 S7
#./bisample_compare_pipline_new.sh -c S67 -v S67 -a ZY-67 -f 5 -s S7 -m S7 -b ZY-7 -d wangyf_S -n 8 &

#while read line;do
#a=`echo $line | awk '{print $1}'`
#b=`echo $line | awk '{print $2}'`
#c=`echo $line | awk '{print $4}'`
#./bisample_compare_pipline.sh -c CS -v CS -a CS -f 5 -s ${a} -m ${c} -b ${b} -d field_cultivar -n 8 &
#wait_all
#done < field_cultivar/metadata_cultivar_CS.txt
#
#while read line;do
#a=`echo $line | awk '{print $1}'`
#b=`echo $line | awk '{print $2}'`
#c=`echo $line | awk '{print $4}'`
#./bisample_compare_pipline.sh -c C46 -v AK58 -a S11 -f 5 -s ${a} -m ${c} -b ${b} -d field_cultivar -n 8 &
#wait_all
#done < field_cultivar/metadata_cultivar_AK.txt
#./bisample_compare_pipline.sh -c CS -v CS -a CS -f 5 -s C46 -m AK58 -b S11 -d cs_ak_fi_yu -n 8 &
#./bisample_compare_pipline.sh -c CS -v CS -a CS -f 5 -s PH09 -m Fielder -b Fielder -d cs_ak_fi_yu -n 8 &
#./bisample_compare_pipline.sh -c C46 -v AK58 -a S11 -f 5 -s PH09 -m Fielder -b Fielder -d cs_ak_fi_yu -n 8 &

./bisample_compare_pipline.sh -c YM01  -v YC_Gao5 -a G5-1  -f 5 -s YM04 -m YC_5224 -b s5224-1 -d lijh -n 8 &
./bisample_compare_pipline.sh -c YM02  -v 981 -a N981-1  -f 5 -s YM03 -m 3097 -b s3097-1 -d lijh -n 8 &


