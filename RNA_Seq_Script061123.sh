#!/bin/bash

# Tmpfs


# User input
while true;
do 

read -p "enter input parameters file path: " -i "/data/" -e par

if [[ -f "$par" ]]; then
	echo "Located input parameters file: $par"
	break
fi
echo Could not find file
done

# Get parameters
readarray -t smp < $par 
# Get the input directory
inDr=${smp[1]}
echo "Input directory: $inDr"
# Get output directory
outDr=${smp[3]}
echo "Output directory: $outDr"
# BAM out 
bamOut=${smp[5]}
echo "Output BAM file: $bamOut"
# Adaptor file
adaptFn=${smp[7]}
echo "Adaptor file name: $adaptFn"
# Sample extension
smpExt=${smp[9]}
smpExtNoGz=${smpExt/.gz/}
echo smpExtNoGz = $smpExtNoGz
echo "First Sample path: $inDr${smp[11]}1$smpExt"

# Check that first sample exists

if [[ -f "$inDr${smp[11]}1$smpExt" ]]; then
	echo "Located First sample"
	else
	echo "Could not find First Sample, Exiting"
exit 1

fi


# Should probably have something that checks for directories in the file name and creates them in the results

# Make a reports directory
mkdir -p $outDr/reports
mkdir -p $outDr/reports/raw

# Loop through samples and apply steps
for ((i=11;i<${#smp[@]};++i)); do

unset dr
dr=($(echo ${smp[i]} | tr "/" "\n")) # need to split directory for fastqc

if [[ ${#dr[@]} > 1 ]]; then

 mkdir -p $outDr/${dr[0]}
	
 	exp=${dr[1]}
	tDr=${dr[0]}
	else
	exp=${dr[0]}
	tDr=
fi 
echo $exp
pth=$inDr
R1=${smp[i]}1$smpExt
R2=${smp[i]}2$smpExt
echo R1=$R1
echo R2=$R2

# Fastqc raw data
docker run -v $inDr/:/data/ --rm -v $outDr/reports/raw/:/rep/ biocontainers/fastqc:v0.11.5_cv3 fastqc --outdir=/rep/ /data/$R1
docker run -v $inDr/:/data/ --rm -v $outDr/reports/raw/:/rep/ biocontainers/fastqc:v0.11.5_cv3 fastqc --outdir=/rep/ /data/$R2


# scythe and sickle

if [[ -z "$adaptFn" ]]; then

echo "No adpator file go to sickle trimming"


pth=$inDr

else 

echo -e "Found adaptor file performing scythe then sickle trimming \n\n"

# Scythe on first read pair (remove adpator contamination)
docker run --rm -v /$inDr/:/data/ \
-v /$outDr/:/out/ \
sickle_scythe_dj scythe \
-a /data/$adaptFn \
-o /out/${R1/$smpExt/_Scythe$smpExtNoGz} \
-m /out/${exp}matches.txt \
/data/$R1 > $outDr/reports/${exp}1_scythe_out.txt

# Scythe on second read pair
docker run --rm -v /$inDr/:/data/ \
-v /$outDr/:/out/ \
sickle_scythe_dj scythe \
-a /data/$adaptFn \
-o /out/${R2/$smpExt/_Scythe$smpExtNoGz} \
-m /out/${exp}matches.txt \
/data/$R2 > $outDr/reports/${exp}2_scythe_out.txt

pigz -p12 $outDr${R1/$smpExt/_Scythe$smpExtNoGz}
pigz -p12 $outDr${R2/$smpExt/_Scythe$smpExtNoGz}

R1=${R1/$smpExt/_Scythe$smpExt}
R2=${R2/$smpExt/_Scythe$smpExt}

pth=$outDr

echo -e "\n R1 after adaptor=$R1 \n"

fi

echo -e "\n Run sickle \n"

## Run sickle (remove low quality bases)

docker run --user $(id -u):1001  --rm -it \
-v $pth:/data/ \
-v $outDr:/out/ \
sickle_scythe_dj sickle pe \
-f /data/$R1 \
-r /data/$R2 \
-o /out/${R1/$smpExt/_Sickle$smpExtNoGz} \
-p /out/${R2/$smpExt/_Sickle$smpExtNoGz} \
-s /out/${R1/$smpExt/_SickleSfile$smpExtNoGz} \
-t sanger \
-n > $outDr/reports/${exp}_sickle_out.txt

RS=${R1/$smpExt/_SickleSfile$smpExtNoGz}
R1=${R1/$smpExt/_Sickle$smpExtNoGz}
R2=${R2/$smpExt/_Sickle$smpExtNoGz}

## Remove the lines from the sickle out whcih mu
sed '1,6d' $outDr/reports/${exp}_sickle_out.txt > $outDr/reports/${exp}_sickle_out.log


#echo $RS
#&&
#sickle
## Run fastqc
docker run -v $outDr/:/data/ --rm -v $outDr/reports/:/rep/ biocontainers/fastqc:v0.11.5_cv3 fastqc --outdir=/rep/ /data/$R1
docker run -v $outDr/:/data/ --rm -v $outDr/reports/:/rep/ biocontainers/fastqc:v0.11.5_cv3 fastqc --outdir=/rep/ /data/$R2


pth=$outDr
# zip output files
pigz -p12 $pth$R1
pigz -p12 $pth$R2
pigz -p12 $pth$RS
R1=${R1}.gz
R2=${R2}.gz

echo -e "\n Run STAR \n"

# STAR

if [[ $bamOut = "yes" ]]; then 
	wrtBam="Full"
else
	wrtBam="None"
fi

echo "Starting STAR \n"
echo "Reads 1: $R1 Reads 2 : $R2 SAMmode : $wrtBam"

docker run --rm --user $(id -u):1001 -it  \
-v /userhome/data2_Backup/Jason/Genome_Directory_STAR/:/alignfiles/ \
-v $pth:/inputfiles/ stevetsa/star STARlong \
--runThreadN 21 \
--runMode alignReads \
--outFilterMultimapNmax 1 \
--outFilterMatchNmin 35 \
--quantMode GeneCounts \
--twopassMode Basic \
--readFilesIn /inputfiles/$R1 \
/inputfiles/$R2 \
--readFilesCommand zcat \
--outFileNamePrefix /inputfiles/$tDr/$exp \
--sjdbGTFfile /alignfiles/gencode.v30.chr_patch_hapl_scaff.annotation.gff3 \
--genomeDir /alignfiles/ \
--outSAMmode $wrtBam \
--outSAMtype BAM SortedByCoordinate


done

conda activate multiqc
multiqc $outDr $outDr/reports -o $outDr/reports
multiqc $outDr/reports/raw -o $outDr/reports/raw

