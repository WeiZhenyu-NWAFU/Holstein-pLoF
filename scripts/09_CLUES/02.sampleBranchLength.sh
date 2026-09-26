#/bin/sh

chrom=${1}
pos=${2}

/storage/public/home/2020110005/liuanguo/biosoftware/relate_v1.2.1/scripts/SampleBranchLengths/SampleBranchLengths.sh \
        -i /storage/public/home/2021050457/01.Cattle/01.ARG_beef/04.556-samples_paint_Ne/11.reinfer_Holstein/01.reinfer_holstein/Chr${chrom}_Holstein_reinfer \
        -o clues_infer/${chrom}_${pos}/${chrom}_${pos} \
        -m 1.26e-8 \
        --coal /storage/public/home/2021050457/01.Cattle/01.ARG_beef/04.556-samples_paint_Ne/11.reinfer_Holstein/01.reinfer_holstein/Chr${chrom}_Holstein_reinfer.coal \
        --format b \
        --first_bp ${pos} \
        --last_bp ${pos} \
        --seed 1 \
        --num_samples 200
