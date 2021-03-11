touch jm22_result
echo 'TGGAAAGACTATGGACCTCG_1' >> jm22_result
zcat C4-2_HHNG7ALXX_L5_1.fq.gz | grep 'TGGAAAGACTATGGACCTCG' >> jm22_result
echo 'TGGAAAGACTATGGACCTCG_2' >> jm22_result
zcat C4-2_HHNG7ALXX_L5_2.fq.gz | grep 'TGGAAAGACTATGGACCTCG' >> jm22_result
echo 'AGCAGGTAATCCACACCAAA_1' >> jm22_result
zcat C4-2_HHNG7ALXX_L5_1.fq.gz | grep 'AGCAGGTAATCCACACCAAA' >> jm22_result
echo 'AGCAGGTAATCCACACCAAA_2' >> jm22_result
zcat C4-2_HHNG7ALXX_L5_2.fq.gz | grep 'AGCAGGTAATCCACACCAAA' >> jm22_result
echo 'AATCCACACCAAATCCCCAT_1' >> jm22_result
zcat C4-2_HHNG7ALXX_L5_1.fq.gz | grep 'AATCCACACCAAATCCCCAT' >> jm22_result
echo 'AATCCACACCAAATCCCCAT_2' >> jm22_result
zcat C4-2_HHNG7ALXX_L5_2.fq.gz | grep 'AATCCACACCAAATCCCCAT' >> jm22_result
