
maindir=`pwd`
for j in $( find . -wholename "*R1*.fastq.gz" -type f); do
 cd `dirname $j`
 mv `basename $j` `dirname $j`.R1.fq.gz
 cd $maindir
 done



maindir=`pwd`
for j in $( find . -wholename "*R2*.fastq.gz" -type f ); do
cd `dirname $j`
mv `basename $j` `dirname $j`.R2.fq.gz
cd $maindir
done
