for i in {1..281};do
for i in {1..281};do
echo 0
done |tr '\n' '\t'|sed s/'\t'$/'\n'/g
done > tmp.dist
if true;then
for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B chr1D chr2D chr3D chr4D chr5D chr6D  chr7D; do
num=$(wc -l < ~/mapping/fieldergenomecompare/1.diff_dev/${CHR}.1M.bed)
for n in $(eval echo {1..$num});do
awk 'NR==FNR{for(i=1; i<=NF;i++){a[FNR,i] = $i; } } NR!=FNR{for (i=1; i<=NF; i++){ printf "%d ",$i+a[FNR,i];} print "" } ' /data2/rawdata2/variant_density/5.cross_sample_200923/raw/${CHR}.${n}.rawdist.dist tmp.dist |sponge tmp.dist
done
done
fi
