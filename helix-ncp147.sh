#!/bin/bash

#SBATCH -J helix
#SBATCH -n 1
#SBATCH -N 1

. /etc/profile.d/modules.sh

cat > nucleic.in <<EOF
parm ncp147.prmtop
trajin md.ncdf 1 last 2
center :1-998
image center familiar
rms first :1-998
strip !(:1-294)
trajout nucleic.pdb pdb
go
EOF
mpirun -np 1 cpptraj -i nucleic.in
rm nucleic.in
awk '/ATOM/' < nucleic.pdb > nucleic1.pdb
rm nucleic.pdb
mkdir split_pdb
cd split_pdb
cp ../nucleic1.pdb .
rm ../nucleic1.pdb
awk -v count=9346 'BEGIN{i=1} {print $0 > sprintf("%d.pdb",i);if (NR>=i*count) {close(sprintf("%d.pdb",i));i++;}}' nucleic1.pdb
declare -i N=`cat nucleic1.pdb | wc -l`/9346
rm nucleic1.pdb
cd ..
mkdir helix 
cd helix
for ((i = 1; i <= $N; i = i + 1))
do
echo "Processing configuration $i... for step 1"
/home/blei/cur+/Cur+<<!
&inp file=../split_pdb/${i}.pdb, lis=${i}, lib=/home/blei/cur+/standard, isym=2, &end
2 1 -1 0 0
1:147
294:148
!
done
rm -f *cda*
rm -r split_pdb
for ((i = 1; i <= $N; i = i + 1))
do
echo "Processing configuration $i... for step 2"
awk 'NR >= 340 && NR <= 485' ${i}.lis | awk '{a+=$6}{b+=$7}{c+=$9}{d+=$10}{e+=$11}{f+=$12}END{printf "%f %f %f %f %f %f\n",a/NR,b/NR,c/NR,d/NR,e/NR,f/NR}' >> helix_data.dat
done
echo "Done!"
