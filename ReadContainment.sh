readsFn="/home/mkirsche/ecoli/oxford.fasta"
pafOutputFn="/home/mkirsche/ecoli/oxford.paf"
filteredReadsOutputFn="/home/mkirsche/ecoli/oxford_megamap.fasta"
groundTruthFn="/home/mkirsche/ecoli/oxford_against_ref.paf.filtered"

javac *.java
java -Xmx6G alignment -in $readsFn -out $pafOutputFn
java containment $pafOutputFn $readsFn $filteredReadsOutputFn
java FileIntersect $groundTruthFn $filteredReadsOutputFn
