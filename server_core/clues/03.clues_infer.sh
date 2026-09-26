#/bin/sh
source /storage/public/home/2021050457/software/01.conda_env/py39_tskit_bft/bin/activate
#source /storage/public/home/2021050457/conda_initialize.sh

chrom=${1}
pos=${2}

python /storage/public/home/2020110005/liuanguo/biosoftware/clues_bft/inference.py \
        --times clues_infer/${chrom}_${pos}/${chrom}_${pos} \
        --coal /storage/public/home/2021050457/01.Cattle/01.ARG_beef/04.556-samples_paint_Ne/11.reinfer_Holstein/01.reinfer_holstein/Chr${chrom}_Holstein_reinfer.coal \
        --out clues_infer/${chrom}_${pos}/${chrom}_${pos} \
	--burnin 100 \
	--timeBins timeBins.txt
