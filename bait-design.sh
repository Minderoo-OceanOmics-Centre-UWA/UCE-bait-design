#Download genomes

#The following genomes were downloaded with the ncbi-datasets-cli tool
#(ncbi-datasets-cli-15.11)
conda activate ncbi_datasets
datasets download genome accession 
GCF_010645085.2
GCA_029692085.1
GCF_021018895.1
GCF_017976325.1
GCF_013368585.1
GCA_020698575.1
GCF_016859285.1
GCA_020466265.1
GCF_018320785.1
GCF_029633865.1
GCA_013358685.1
GCF_022458985.1
GCA_951216825.1
GCF_008729295.1
GCA_029721515.1
GCF_027409825.1
GCA_027579695.1
GCA_949987555.1
GCA_024500275.1
GCA_009823395.1

#clean up genomes 
#Remove everything but the asseccion number information for each contig and scaffold. 
#Rename the genomes to contain only the first 3 letters of the Genus and the first 4 letters of the species. 

conda activate biopython 
#(biopython-1.81)

from Bio import SeqIO
with open("Name_of_Genome_File.fna", "rU") as infile:
with open("outfileName.fasta", "w") as outf:
    for seq in SeqIO.parse(infile, 'fasta'):
        seq.name = ""
        seq.description = ""
        outf.write(seq.format('fasta'))


#Convert genomes to 2bit format
conda fatotwobit
(ucsc-fatotwobit-377) 
for critter in *; do faToTwoBit $critter/$critter.fasta $critter/${critter%.*}.2bit; done


# simulate 100bp reads at 2X, insert size 200bp, SD 150. 
conda activate art 
#(art-2016.06.05)

art_illumina \
    --paired \
    --in ~/path/to/input/genome.fasta \
    --out prefix-of-output-file \
    --len 100 --fcov 2 --mflen 200 --sdev 150 -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na

# merge paired end reads

for critter in AciRuth BorAnta CheRost EleAnta LampInc MugCeph PhyEque ScoJapo SipTubi ThaSept AntMacu BosSine CheUndu GymJava MegCypr NeoGrae PleLeop SerAure SynVerr XipGlad;
    do
        echo "working on $critter";
        touch $critter-pe100-reads.fq;
        cat $critter-pe100-reads1.fq > $critter-pe100-reads.fq;
        cat $critter-pe100-reads2.fq >> $critter-pe100-reads.fq;
        rm $critter-pe100-reads1.fq;
        rm $critter-pe100-reads2.fq;
        gzip $critter-pe100-reads.fq;
    done;

#Prepare the Base Genome
#ScoJapo (Scomber japonicus)

## to get the the docker container Sorking and writing to the correct place (docker run -v $(pwd):$(pwd) -u $(id -u $USER):$(id -g $USER) -w $PWD)

docker run -v $(pwd):$(pwd) -u $(id -u $USER):$(id -g $USER) -w $PWD markderbyshire/stampy:latest stampy.py --species="Scomberjaponicus" --assembly="ScoJapo" -G ScoJapo /home/lhuet/bait_project/base/ScoJapo.fasta
docker run -v $(pwd):$(pwd) -u $(id -u $USER):$(id -g $USER) -w $PWD markderbyshire/stampy:latest stampy.py -g ScoJapo -H ScoJapo

#Align reads to base genome using stampy

export cores=32
export base=ScoJapo
export base_dir=$HOME/bait_project/alignments
for critter in AciRuth BorAnta CheRost EleAnta LampInc MugCeph PhyEque SipTubi ThaSept AntMacu BosSine CheUndu GymJava MegCypr NeoGrae PleLeop SerAure SynVerr XipGlad;
    do
        export reads=$critter-pe100-reads.fq.gz;
        mkdir -p $base_dir/$critter;
        cd $base_dir/$critter;
        docker run -v /home:/home -v $(pwd):$(pwd) -u $(id -u $USER):$(id -g $USER) -w $PWD markderbyshire/stampy:latest stampy.py --maxbasequal 93 -g /home/lhuet/bait_project/base/$base -h /home/lhuet/bait_project/base/$base \
        --substitutionrate=0.05 -t$cores --insertsize=400 -M \
        /home/lhuet/bait_project/reads/$reads | samtools sort -o $critter-to-$base.bam 
        cd $base_dir
    done;

# Remove unmapped reads from bam files 
(samtools 1.17)

for critter in AciRuth BorAnta CheRost EleAnta LampInc MugCeph PhyEque SipTubi ThaSept AntMacu BosSine CheUndu GymJava MegCypr NeoGrae PleLeop SerAure SynVerr XipGlad;
    do
        samtools view -h -F 4 -b $critter/$critter-to-ScoJapo.bam > $critter/$critter-to-ScoJapo-MAPPING.bam;
        rm $critter/$critter-to-ScoJapo.bam;
        ln -s ../$critter/$critter-to-ScoJapo-MAPPING.bam all/$critter-to-ScoJapo-MAPPING.bam;
    done;

## Conserved locus identificaion 

# BAM to bed 
# make a directory to hold the BED data
> mkdir bed
> cd bed

conda activate bedtools 
(bedtools-2.31)

# now, convert our *-MAPPING.bam files to BED format
for i in ../alignments/all/*.bam; do echo $i; bedtools bamtobed -i $i -bed12 > `basename $i`.bed; done

#Sort the beds by chromosome/scaffold/contig and position along that chromosome/scaffold/contig
for i in *.bed; do echo $i; bedtools sort -i $i > ${i%.*}.sort.bed; done

# Merge overlapping intervals incase of small gaps between proximate regions of conservation
for i in *.bam.sort.bed; do echo $i; bedtools merge -i $i > ${i%.*}.merge.bed; done


# How many merged putatively conserved regions in each exemplar taxon are shared with the base genome?

for i in *.bam.sort.merge.bed; do wc -l $i; done


## Remove repetitive intervals

conda activate phyluce-1.7.3

for i in *.sort.merge.bed;
    do
        phyluce_probe_strip_masked_loci_from_set \
            --bed $i \
            --twobit ../genomes/ScoJapo/ScoJapo.2bit \
            --output ${i%.*}.strip.bed \
            --filter-mask 0.25 \
            --min-length 80
    done;

phyluce_probe_get_multi_merge_table \
    --conf bed-files.conf \
    --base-taxon ScoJapo \
    --output Actinopteri-to-ScoJapo.sqlite

sqlite> select * from ScoJapo limit 10


## determine shared conserved loci 

phyluce_probe_query_multi_merge_table \
        --db Actinopteri-to-ScoJapo.sqlite  \
        --base-taxon ScoJapo

phyluce_probe_query_multi_merge_table \
        --db Actinopteri-to-ScoJapo.sqlite \
        --base-taxon ScoJapo \
        --output ScoJapo+19.bed \
        --specific-counts 19

Counter({'aciruth': 10341, 'boranta': 10341, 'cherost': 10341, 'eleanta': 10341, 'lampinc': 10341, 'mugceph': 10341, 'phyeque': 10341, 'siptubi': 10341, 'thasept': 10341, 'antmacu': 10341, 'bossine': 10341, 'cheundu': 10341, 'gymjava': 10341, 'megcypr': 10341, 'neograe': 10341, 'pleleop': 10341, 'seraure': 10341, 'synverr': 10341, 'xipglad': 10341})

## Design baits to target them

phyluce_probe_get_genome_sequences_from_bed \
    --bed ScoJapo+19.bed  \
    --twobit ../genomes/ScoJapo/ScoJapo.2bit \
    --buffer-to 160 \
    --output ScoJapo+19.fasta

Screened 10341 sequences.  Filtered 10 < 160 bp or with > 25.0% masked bases or > 0 N-bases. Kept 10331.

# Design temp bait sets 

phyluce_probe_get_tiled_probes \
    --input ScoJapo+19.fasta \
    --probe-prefix "uce-" \
    --design Actinopteri-v1 \
    --designer huet \
    --tiling-density 3 \
    --two-probes \
    --overlap middle \
    --masking 0.25 \
    --remove-gc \
    --output ScoJapo+19.temp.probes

## Remove duplicates

phyluce_probe_easy_lastz \
    --target ScoJapo+19.temp.probes \
    --query ScoJapo+19.temp.probes \
    --identity 50 --coverage 50 \
    --output ScoJapo+19.temp.probes-TO-SELF-PROBES.lastz

phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
    --fasta ScoJapo+19.temp.probes  \
    --lastz ScoJapo+19.temp.probes-TO-SELF-PROBES.lastz \
    --probe-prefix=uce-

Screened 20534 fasta sequences.  Filtered 2348 duplicates. Kept 15848.

# Run some insilico tests 

> mkdir Actinopteri-genome-lastz

phyluce_probe_run_multiple_lastzs_sqlite \
    --probefile ../bed/ScoJapo+19.temp-DUPE-SCREENED.probes \
    --scaffoldlist AciRuth BorAnta CheRost EleAnta LampInc MugCeph PhyEque SipTubi ThaSept AntMacu BosSine CheUndu GymJava MegCypr NeoGrae PleLeop SerAure SynVerr XipGlad ScoJapo \
    --genome-base-path ../genomes \
    --identity 50 \
    --cores 32 \
    --db ScoJapo+19.sqlite \
    --output Actinopteri-genome-lastz


# extract sequence around conserved loci from exemplar genomes

phyluce_probe_slice_sequence_from_genomes \
    --conf Actinopteri-genome.conf \
    --lastz Actinopteri-genome-lastz \
    --probes 180 \
    --name-pattern "ScoJapo+19.temp-DUPE-SCREENED.probes_v_{}.lastz.clean" \
    --output Actinopteri-genome-fasta

phyluce_probe_get_multi_fasta_table \
    --fastas ../Actinopteri-genome-fasta\
    --output multifastas.sqlite \
    --base-taxon ScoJapo

phyluce_probe_query_multi_fasta_table \
    --db multifastas.sqlite \
    --base-taxon ScoJapo

phyluce_probe_query_multi_fasta_table \
    --db multifastas.sqlite \
    --base-taxon ScoJapo \
    --output ScoJapo+19-back-to-15.conf \
    --specific-counts 15

#Final bait design

phyluce_probe_get_tiled_probe_from_multiple_inputs \
    --fastas ../Actinopteri-genome-fasta \
    --multi-fasta-output ScoJapo+19-back-to-15.conf \
    --probe-prefix "uce-" \
    --designer huet \
    --design actinopteri-v1 \
    --tiling-density 3 \
    --overlap middle \
    --masking 0.25 \
    --remove-gc \
    --two-probes \
    --output Actinopteri-v1-master-probe-list.fasta

2023-08-11 04:32:39,425 - phyluce_probe_get_tiled_probe_from_multiple_inputs - INFO - Conserved locus count = 3302
2023-08-11 04:32:39,425 - phyluce_probe_get_tiled_probe_from_multiple_inputs - INFO - Probe Count = 114972
2023-08-11 04:32:39,566 - phyluce_probe_get_tiled_probe_from_multiple_inputs - INFO - == Completed phyluce_probe_get_tiled_probe_from_multiple_inputs =

# remove dups 

phyluce_probe_easy_lastz \
    --target Actinopteri-v1-master-probe-list.fasta \
    --query Actinopteri-v1-master-probe-list.fasta \
    --identity 50 \
    --coverage 50 \
    --output Actinopteri-v1-master-probe-list-TO-SELF-PROBES.lastz



phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
    --fasta Actinopteri-v1-master-probe-list.fasta \
    --lastz Actinopteri-v1-master-probe-list-TO-SELF-PROBES.lastz \
    --probe-prefix=uce-

In silico tests 

phyluce_probe_run_multiple_lastzs_sqlite \
    --db ScoJapo-test.sqlite \
    --output Actinopteri-genome-lastz \
    --probefile ../probe-design/Actinopteri-v1-master-probe-list-DUPE-SCREENED.fasta \
    --scaffoldlist AciRuth BorAnta CheRost EleAnta LampInc MugCeph PhyEque SipTubi ThaSept AntMacu BosSine CheUndu GymJava MegCypr NeoGrae PleLeop SerAure SynVerr XipGlad ScoJapo \
    --genome-base-path ../genomes \
    --identity 50 \
    --cores 32
