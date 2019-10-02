maindir=`pwd`;
for i in $( find . -maxdepth 2 -wholename "*_htseq.cnt" -type f ); do
cp $maindir/`dirname $i`/`basename $i` $maindir/HtSeqCounts
done;