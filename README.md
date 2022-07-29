# Transponson Element
My own transponson element analysis pipeline
## transposon element identification
Using [EDTA](https://github.com/oushujun/EDTA) identification
~~~shell
singularity exec -B /mnt:/mnt /home/assembly/tools/EDTA.sif EDTA.pl \
--genome/path/to/genome.fasta \
--cds /path/to/genome.cds.fa \
--sensitive 1 --anno 1
~~~
Get the [RepeatMasker](https://www.repeatmasker.org/) output from EDTA out put,the path is `*.mod.EDTA.final/*.mod.EDTA.raw.fa.preFriJul152239272022.RMoutput/*.mod.EDTA.raw.fa.out`
~~~shell
awk '!/\*/' RM.out > RM_nostar.out
~~~
using [RM_TRIPS](https://github.com/clbutler/RM_TRIPS) clean RM results
and you will get `RM_tidy.out` 

## TE classification (optional)
If you have lots of unknown TEs, I recommend use [DeepTE](https://github.com/LiLabAtVT/DeepTE) to reannotate your TE sequences(`genome.fasta.mod.EDTA.TElib.fa` in EDTA results)
~~~shell
DeepTE.py -d working_dir -o output_dir -i TElib.fa -sp F -m_dir /path/to/fungi/model
~~~
if your organism is plant, I also recommend usng [TEsorter](https://github.com/zhangrengang/TEsorter)

## Genomic distribution
Calculate LTR genomic distribution. I split the genomic region into four parts:intergenetic, upstream of downstream 2kb of gene, upstream of downstream 2~5kb of gene, inside gene.
~~~shell
python ./TE/Distribution.py -i RM_tidy.out -o RM_tidy_distribution.out -g /path/to/genomic.gff
~~~
*ps:* Make ensure the gff file is correct  

## Insert time
Only intact LTR could calculate insert time. Using EDTA results.
~~~shell
python Insert_time.py -i /path/to/EDTA/genome.fasta.mod.EDTA.raw/genome.fasta.mod.LTR.intact.gff3 -o insert_time.lst
~~~

## Intact, truncated and solo LTR
Please refer to reference 1 if you don't know what is intact, truncated and solo LTR and why calculate their ration. You need download sequence from [GyDB](https://gydb.org/index.php?title=Main_Page). The url for the download data is https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip.   
**Requierment:** `BLAST` and `silix`
~~~shell
unzip GyDB_collection.zip
cat GyDB_collection/consensus/GAG_* > gag.fa
python solo_detect.py -g genome.fasta -i ../genome.fasta.mod.EDTA.raw/genome.fasta.mod.LTR.intact.fa --gag gag.fa -t 10 -o out.tsv
~~~
*ps:* EDTA will release a version that can calculate solo LTR directly soon.

## Reference
1. Lyu, H., He, Z., Wu, C.-I. and Shi, S. (2018), Convergent adaptive evolution in marginal environments: unloading transposable elements as a common strategy among mangrove genomes. New Phytol, 217: 428-438. https://doi.org/10.1111/nph.14784
2. Miele, V., Penel, S. & Duret, L. Ultra-fast sequence clustering from similarity networks with SiLiX. BMC Bioinformatics 12, 116 (2011). https://doi.org/10.1186/1471-2105-12-116
3. Ou S., Su W., Liao Y., Chougule K., Agda J. R. A., Hellinga A. J., Lugo C. S. B., Elliott T. A., Ware D., Peterson T., Jiang N.✉, Hirsch C. N.✉ and Hufford M. B.✉ (2019). Benchmarking Transposable Element Annotation Methods for Creation of a Streamlined, Comprehensive Pipeline. Genome Biol. 20(1): 275.
4. Zhang RG, Li GL, Wang XL et. al. TEsorter: an accurate and fast method to classify LTR retrotransposons in plant genomes. Horticulture Research, 2022, 9: uhac017 https://doi.org/10.1093/hr/uhac017
5. Llorens C, Futami R, Covelli L et. al. The Gypsy Database (GyDB) of mobile genetic elements: release 2.0. Nucleic Acids Research, 2011, 39: 70–74 https://doi.org/10.1093/nar/gkq1061
6. Bell, E. A., Butler, C. L., Oliveira, C., Marburger, S., Yant, L., & Taylor, M. I. (2022). Transposable element annotation in non-model species: The benefits of species-specific repeat libraries using semi-automated EDTA and DeepTE de novo pipelines. Molecular Ecology Resources, 22, 823– 833. https://doi.org/10.1111/1755-0998.13489