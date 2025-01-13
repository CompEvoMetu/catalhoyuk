#!/bin/bash -l

######## Count kmers length of 15 for each fastq file
for i in $(cat fastqlist)
do

~/bin/jellyfish count -s 1000000000 -C -m 15 -t 30 <(zcat ${i}.fastq.gz)
~/bin/jellyfish dump mer_counts.jf | perl -pe 's/>(.*)/>\1\t/g; s/\n//g; s/>/\n/g' | grep -v '^$' | awk '{print $2, $1}' > kmer_counts/${i}_15mer
echo $i

done


######## Kmer count files into one
hdr() { awk 'FNR==1{ print "kmer", FILENAME }1' "$1"; }

join -a1 -a2 -e 0 -o auto <(hdr fastq1 ) <(hdr fastq2) > join.tmp

while read -r line; do
   join -a1 -a2 -e 0 -o auto join.tmp <(hdr "$line") > join.tmp.1
   mv join.tmp.1 kmer_15
   echo $line
done < fastq_list


######## Randomly sample kmer counts
python3 ./shuf-master/powershuf.py -n 1000 --file ../kmer_counts/join.tmp > kmer_1K
python3 ./shuf-master/powershuf.py -n 10000 --file ../kmer_counts/join.tmp > kmer_10K
python3 ./shuf-master/powershuf.py -n 50000 --file ../kmer_counts/join.tmp > kmer_50K
python3 ./shuf-master/powershuf.py -n 100000 --file ../kmer_counts/join.tmp > kmer_100K
python3 ./shuf-master/powershuf.py -n 500000 --file ../kmer_counts/join.tmp > kmer_500K
python3 ./shuf-master/powershuf.py -n 1000000 --file ../kmer_counts/join.tmp > kmer_1M
python3 ./shuf-master/powershuf.py -n 5000000 --file ../kmer_counts/join.tmp > kmer_5M
python3 ./shuf-master/powershuf.py -n 10000000 --file ../kmer_counts/join.tmp > kmer_10M
python3 ./shuf-master/powershuf.py -n 50000000 --file ../kmer_counts/join.tmp > kmer_50M