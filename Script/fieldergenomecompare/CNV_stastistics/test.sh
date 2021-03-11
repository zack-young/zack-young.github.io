for CHR in  chr1A chr2A chr3A chr4A chr5A chr6A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
num=1
(for CHR1 in  chr1A chr2A chr3A chr4A chr5A chr6A chr6A chr7A chr1B chr2B chr3B chr4B chr5B chr6B chr7B chr1D chr2D chr3D chr4D chr5D chr6D chr7D; do
num=$(( $num + 1 ))
echo ${CHR} ${CHR1} ${num}
done) &
done
