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

    phyluce_probe_slice_sequence_from_genomes \
    --conf ~/bait_project/probe-design/Actinopteri-genome.conf \
    --lastz Actinopteri-genome-lastz\
    --output Actinopteri-genome-fasta \
    --flank 500 \
    --name-pattern "Actinopteri-v1-master-probe-list-DUPE-SCREENED.fasta_v_{}.lastz.clean"

# match contigs to baits
phyluce_assembly_match_contigs_to_probes \
    --contigs Actinopteri-genome-fasta \
    --probes ~/bait_project/probe-design/Actinopteri-v1-master-probe-list-DUPE-SCREENED.fasta \
    --output in-silico-genomes-lastz \
    --min-coverage 67 \
    --log-path log

# get match counts and extract fasta informaiton 

phyluce_assembly_get_match_counts \
    --locus-db in-silico-genomes-lastz/probe.matches.sqlite \
    --taxon-list-config taxon-set.conf \
    --taxon-group 'all' \
    --output taxon-sets/all/all-taxa-incomplete.conf \
    --log-path log \
    --incomplete-matrix

#extract fasta info for each  locus into a monolithic FASTA file

phyluce_assembly_get_fastas_from_match_counts \
    --contigs ~/bait_project/probe-design-test/Actinopteri-genome-fasta \
    --locus-db ../../in-silico-genomes-lastz/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log

#explode monolithic fasta file 

phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas-uce \
    --by-taxon

# get summary stats on the FASTAS
for i in exploded-fastas-uce/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done
aciruth.unaligned.fasta,188,212427,1129.9308510638298,4.981391304812315,759,1184,1157.0,175
antmacu.unaligned.fasta,2551,2932895,1149.7040376323012,1.097558139155687,660,1473,1161.0,2481
boranta.unaligned.fasta,2328,2619029,1125.0124570446735,2.092735197897589,137,1431,1160.0,2125
bossine.unaligned.fasta,2557,2922455,1142.9233476730544,1.2912783258060432,391,1442,1160.0,2471
cherost.unaligned.fasta,2547,2957275,1161.0816647035729,0.5140541884709552,760,1472,1162.0,2537
cheundu.unaligned.fasta,2577,2978979,1155.987194412107,0.6787568177166526,707,1473,1161.0,2556
eleanta.unaligned.fasta,2216,2462703,1111.3280685920577,2.4639269219274738,190,1472,1159.0,1959
gymjava.unaligned.fasta,1924,2194093,1140.380977130977,1.519437078350043,539,1255,1160.0,1843
lampinc.unaligned.fasta,2540,2854281,1123.7326771653543,1.8798676910468195,460,1473,1160.0,2301
megcypr.unaligned.fasta,1373,1578871,1149.9424617625637,1.5000766984325935,633,1631,1160.0,1347
mugceph.unaligned.fasta,2562,2965597,1157.532006245121,0.7067515776982285,460,1555,1161.0,2544
neograe.unaligned.fasta,2340,2510476,1072.8529914529915,2.9285891930762484,273,1402,1145.0,1828
phyeque.unaligned.fasta,2542,2896137,1139.3143194335169,1.2409871992284993,656,1551,1160.0,2445
pleleop.unaligned.fasta,2485,2856448,1149.4760563380282,1.1372456780089713,182,1950,1161.0,2439
scojapo.unaligned.fasta,2577,2969717,1152.3930927435001,1.2644450693981824,306,1500,1162.0,2518
seraure.unaligned.fasta,2552,2959482,1159.6716300940438,0.5639492743532605,659,1540,1161.0,2543
siptubi.unaligned.fasta,2338,2651470,1134.076133447391,1.60329798744617,337,1402,1160.0,2204
synverr.unaligned.fasta,2556,2946877,1152.9252738654147,0.9682932778755361,737,1472,1161.0,2489
thasept.unaligned.fasta,2568,2953024,1149.9314641744547,0.8355905669327105,682,1363,1160.0,2530
xipglad.unaligned.fasta,2563,2970101,1158.837690206789,0.60118901114323,746,1473,1161.0,2552

#concatenate the files
cat exploded-fastas-uce/* > insilico-incomplete.fasta

#explode the fastas 

phyluce_assembly_explode_get_fastas_file \
    --input insilico-incomplete.fasta \
    --output exploded-uce-loci

#make mafft directory 
mkdir mafft_alns log

# remove single seq from exploded_uce_loci 
grep -c ">" *.fasta | awk -F ':' '$2==1{print $1}' | xargs -I {} mv {} ../exploded_uce_loci_single_seq


for i in *.unaligned.fasta; do
mafft --auto --adjustdirection --thread 20 ${i} > ../mafft_alns/${i%.*}.aligned.fasta;
done

#  rename the alignments before cleaning and subsettings - remove "unaligned" from the file name 
## CONVERT FROM FASTA TO NEXUS

cd mafft_alns
for file in *unaligned.aligned.fasta ; do mv $file ${file//unaligned.aligned.fasta/fasta} ; done

# convert from fasta to nexus
phyluce_align_convert_one_align_to_another \
                --alignments mafft_alns \
                --output mafft_alns-nexus \
                --input fasta \
                --output-format nexus \
                --cores 20 \
                --log-path log 

### GET INFORATIVE SITES 
# merge exon + uce and do combo
phyluce_align_get_informative_sites \
                --alignments mafft_alns-nexus \
                --cores 20 \
                --input-format nexus \
                --log-path log \
                --output mafft_alns-nexus_total_informative_sites


phyluce_align_get_smilogram_from_alignments --alignments mafft_alns-nexus  --output-file smilogram_sebastes  --output-database data --log log --output-missing data




### testing removing long unformaitve baits 

awk '($3 - $2) <= 200' ScoJapo+19.bed > ScoJapo+19+200.bed


## Design baits to target them

phyluce_probe_get_genome_sequences_from_bed \
    --bed ScoJapo+19+200.bed  \
    --twobit ../genomes/ScoJapo/ScoJapo.2bit \
    --buffer-to 160 \
    --output ScoJapo+19+200.fasta

 Screened 6026 sequences.  Filtered 10 < 160 bp or with > 25.0% masked bases or > 0 N-bases. Kept 6016.

# Design temp bait sets 

phyluce_probe_get_tiled_probes \
    --input ScoJapo+19+200.fasta \
    --probe-prefix "uce-" \
    --design Actinopteri-v1 \
    --designer huet \
    --tiling-density 3 \
    --two-probes \
    --overlap middle \
    --masking 0.25 \
    --remove-gc \
    --output ScoJapo+19+200.temp.probes

2023-10-06 01:42:33,112 - phyluce_probe_get_tiled_probes - INFO - Conserved locus count = 6007
2023-10-06 01:42:33,112 - phyluce_probe_get_tiled_probes - INFO - Probe Count = 11993


## Remove duplicates

phyluce_probe_easy_lastz \
    --target ScoJapo+19+200.temp.probes \
    --query ScoJapo+19+200.temp.probes \
    --identity 50 --coverage 50 \
    --output ScoJapo+19+200.temp.probes-TO-SELF-PROBES.lastz

phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
    --fasta ScoJapo+19+200.temp.probes  \
    --lastz ScoJapo+19+200.temp.probes-TO-SELF-PROBES.lastz \
    --probe-prefix=uce-

Screened 11992 fasta sequences.  Filtered 777 duplicates. Kept 10440.

# Run some insilico tests 

> mkdir Actinopteri-genome-lastz

phyluce_probe_run_multiple_lastzs_sqlite \
    --probefile ../bed/ScoJapo+19+200.temp-DUPE-SCREENED.probes \
    --scaffoldlist AciRuth BorAnta CheRost EleAnta LampInc MugCeph PhyEque SipTubi ThaSept AntMacu BosSine CheUndu GymJava MegCypr NeoGrae PleLeop SerAure SynVerr XipGlad ScoJapo \
    --genome-base-path ../genomes \
    --identity 50 \
    --cores 32 \
    --db ScoJapo+19+200.sqlite \
    --output Actinopteri-genome-lastz


# extract sequence around conserved loci from exemplar genomes

phyluce_probe_slice_sequence_from_genomes \
    --conf ~/bait_project/probe-design/Actinopteri-genome.conf \
    --lastz Actinopteri-genome-lastz \
    --probes 180 \
    --name-pattern "ScoJapo+19+200.temp-DUPE-SCREENED.probes_v_{}.lastz.clean" \
    --output Actinopteri-genome-fasta

cd Actinopteri-genome-lastz

phyluce_probe_get_multi_fasta_table \
    --fastas ../Actinopteri-genome-fasta\
    --output multifastas.sqlite \
    --base-taxon ScoJapo

phyluce_probe_query_multi_fasta_table \
    --db multifastas.sqlite \
    --base-taxon ScoJapo

2023-10-06 03:19:58,414 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 0 taxa: 3,321.0
2023-10-06 03:19:58,414 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 1 taxa: 3,321.0
2023-10-06 03:19:58,414 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 2 taxa: 2,906.0
2023-10-06 03:19:58,414 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 3 taxa: 2,685.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 4 taxa: 2,544.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 5 taxa: 2,438.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 6 taxa: 2,344.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 7 taxa: 2,276.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 8 taxa: 2,207.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 9 taxa: 2,156.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 10 taxa:        2,110.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 11 taxa:        2,066.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 12 taxa:        2,021.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 13 taxa:        1,974.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 14 taxa:        1,916.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 15 taxa:        1,826.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 16 taxa:        1,688.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 17 taxa:        1,451.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 18 taxa:        1,094.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 19 taxa:        576.0
2023-10-06 03:19:58,415 - phyluce_probe_query_multi_fasta_table - INFO - Loci shared by 20 taxa:        78.0


12

phyluce_probe_query_multi_fasta_table \
    --db multifastas.sqlite \
    --base-taxon ScoJapo \
    --output ScoJapo+19-back-to-12.conf \
    --specific-counts 12

Counter({'seraure': 1980, 'cherost': 1975, 'cheundu': 1971, 'synverr': 1970, 'scojapo': 1968, 'mugceph': 1967, 'xipglad': 1961, 'antmacu': 1955, 'bossine': 1944, 'phyeque': 1938, 'thasept': 1930, 'pleleop': 1898, 'lampinc': 1870, 'siptubi': 1797, 'boranta': 1734, 'neograe': 1636, 'eleanta': 1612, 'gymjava': 1406, 'megcypr': 1102, 'aciruth': 241})
2023-10-06 03:21:07,976 - phyluce_probe_query_multi_fasta_table - INFO - Total loci = 2021



# 
#Final bait design

phyluce_probe_get_tiled_probe_from_multiple_inputs \
    --fastas ../Actinopteri-genome-fasta \
    --multi-fasta-output ScoJapo+19-back-to-12.conf \
    --probe-prefix "uce-" \
    --designer huet \
    --design actinopteri-v1 \
    --tiling-density 3 \
    --overlap middle \
    --masking 0.25 \
    --remove-gc \
    --two-probes \
    --output Actinopteri-v1-master-probe-list.fasta

  INFO - Conserved locus count = 2021
  Probe Count = 69085

# remove dups 

phyluce_probe_easy_lastz \
    --target Actinopteri-v1-master-probe-list.fasta \
    --query Actinopteri-v1-master-probe-list.fasta \
    --identity 50 \
    --coverage 50 \
    --output Actinopteri-v2-master-probe-list-TO-SELF-PROBES.lastz



phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
    --fasta Actinopteri-v1-master-probe-list.fasta \
    --lastz Actinopteri-v2-master-probe-list-TO-SELF-PROBES.lastz \
    --probe-prefix=uce-


In silico tests 

phyluce_probe_run_multiple_lastzs_sqlite \
    --db ScoJapo-200test.sqlite \
    --output Actinopteri-genome-lastz \
    --probefile ../probe-design200/Actinopteri-v1-master-probe-list-DUPE-SCREENED.fasta \
    --scaffoldlist AciRuth BorAnta CheRost EleAnta LampInc MugCeph PhyEque SipTubi ThaSept AntMacu BosSine CheUndu GymJava MegCypr NeoGrae PleLeop SerAure SynVerr XipGlad ScoJapo \
    --genome-base-path ../genomes \
    --identity 50 \
    --cores 32

phyluce_probe_slice_sequence_from_genomes \
    --conf ~/bait_project/probe-design/Actinopteri-genome.conf \
    --lastz Actinopteri-genome-lastz\
    --output Actinopteri-genome-fasta \
    --flank 500 \
    --name-pattern "Actinopteri-v1-master-probe-list-DUPE-SCREENED.fasta_v_{}.lastz.clean"

# match contigs to baits
phyluce_assembly_match_contigs_to_probes \
    --contigs Actinopteri-genome-fasta \
    --probes ~/bait_project/probe-design200/Actinopteri-v1-master-probe-list-DUPE-SCREENED.fasta \
    --output in-silico-genomes-lastz \
    --min-coverage 67 \
    --log-path log

# get match counts and extract fasta informaiton 

phyluce_assembly_get_match_counts \
    --locus-db in-silico-genomes-lastz/probe.matches.sqlite \
    --taxon-list-config ~/bait_project/probe-design-test/taxon-set.conf \
    --taxon-group 'all' \
    --output taxon-sets/all/all-taxa-incomplete.conf \
    --log-path log \
    --incomplete-matrix

#extract fasta info for each  locus into a monolithic FASTA file

phyluce_assembly_get_fastas_from_match_counts \
    --contigs ~/bait_project/probe-design200-test/Actinopteri-genome-fasta \
    --locus-db ../../in-silico-genomes-lastz/probe.matches.sqlite \
    --match-count-output all-taxa-incomplete.conf \
    --output all-taxa-incomplete.fasta \
    --incomplete-matrix all-taxa-incomplete.incomplete \
    --log-path log

#explode monolithic fasta file 

phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas-uce \
    --by-taxon

# get summary stats on the FASTAS
for i in exploded-fastas-uce/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done


#concatenate the files
cat exploded-fastas-uce/* > insilico-incomplete.fasta

#explode the fastas 

phyluce_assembly_explode_get_fastas_file \
    --input insilico-incomplete.fasta \
    --output exploded-uce-loci

#make mafft directory 
mkdir mafft_alns log

# remove single seq from exploded_uce_loci 
grep -c ">" *.fasta | awk -F ':' '$2==1{print $1}' | xargs -I {} mv {} ../exploded_uce_loci_single_seq


for i in *.unaligned.fasta; do
mafft --auto --adjustdirection --thread 20 ${i} > ../mafft_alns/${i%.*}.aligned.fasta;
done

#  rename the alignments before cleaning and subsettings - remove "unaligned" from the file name 
## CONVERT FROM FASTA TO NEXUS

cd mafft_alns
for file in *unaligned.aligned.fasta ; do mv $file ${file//unaligned.aligned.fasta/fasta} ; done


### GET INFORATIVE SITES 
# merge exon + uce and do combo
phyluce_align_get_informative_sites \
                --alignments total \
                --cores 20 \
                --input-format nexus \
                --log-path log \
                --output total-edge-50-informative-sites 


phyluce_align_get_smilogram_from_alignments --alignments  --output-file  --output-database data --log log --output-missing data
