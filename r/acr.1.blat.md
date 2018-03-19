## Maize ACR analysis

1. Create a working directory and switch to it
  ```bash
  mkdir -p /home/springer/zhoux379/data/misc1/maize.acr/data
  cd /home/springer/zhoux379/data/misc1/maize.acr/data
  ```

2. Create an alias for the BED file
  ```bash
  ln -sf /home/springer/nosha003/wgbs_schmitz/ACR/B73_peaks/B73L_final_ACR.bed 01.bed
  ```

3. Extract sequences using the coordinates in the BED file from reference genome
  ```bash
  seqret.pl -d /home/springer/zhoux379/data/genome/Zmays_v4/11_genome.fas -b 01.bed -o 02.fas
  ```
  * [link to seqret.pl](https://github.com/orionzhou/luffy/blob/master/perl/seqret.pl)

4. Load the BLAT toolkit (in MSI)
  ```bash
  module load blat/34
  ```

5. Build genome index for BLAT, this step is optional since BLAT can also take a fasta file dire
  ```bash
  faToTwoBit Zmays_v4.fasta Zmays_v4.2bit
  ```
  * [link to faToTwoBit](https://genome.ucsc.edu/goldenpath/help/blatSpec.html#faToTwoBitUsage)

6. Blat against the (already indexed) PH207 genome
  ```bash
  blat /home/springer/zhoux379/data/genome/PH207/21.blat/db.2bit 02.fas 04.psl
  ```
  * [link to blat usage](https://genome.ucsc.edu/goldenpath/help/blatSpec.html#blatUsage)
  * blat output is a 21-column text file in [PSL format](https://useast.ensembl.org/info/website/upload/psl.html)

7. Convert Blat PSL output to more human-readable tab-separated file
  ```bash
  psl2tsv.pl -i 04.psl -o 05.tsv
  ```
  * [link to psl2tsv.pl](https://github.com/orionzhou/luffy/blob/master/perl/psl2tsv.pl)
