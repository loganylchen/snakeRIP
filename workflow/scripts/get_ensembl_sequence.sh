#!/bin/bash
set -x
set -e

exec > "${snakemake_log[0]}" 2>&1

species="${snakemake_params[species]}"
release="${snakemake_params[release]}"
build="${snakemake_params[build]}"
datatype="${snakemake_params[datatype]}"
threads="${snakemake[threads]}"
output_file="${snakemake_output[0]}"



branch=""
if [ "$release" -ge 81 ] && [ "$build" = "GRCh37" ]; then
    branch="grch37/"
elif [ -n "${snakemake_params[branch]}" ]; then
    branch="${snakemake_params[branch]}/"
fi


if [ "$release" -gt 75 ]; then
    spec="$build"
else
    spec="${build}.${release}"
fi


if [[ "$output_file" == *.gz ]]; then
    decompress="cat "
else
    decompress="gzip -dc"
fi


url="https://ftp.ensembl.org/pub"
url_prefix="${url}/${branch}release-${release}/fasta/${species}/${datatype}/${species^}.${spec}"


tmpdir=$(mktemp -d -p ./)
trap 'rm -rf "$tmpdir"' EXIT



declare -a suffixes
if [ "$datatype" = "dna" ]; then
    suffixes=("dna.primary_assembly.fa.gz" "dna.toplevel.fa.gz")
elif [ "$datatype" = "cdna" ]; then
    suffixes=("cdna.all.fa.gz")
elif [ "$datatype" = "cds" ]; then
    suffixes=("cds.all.fa.gz")
elif [ "$datatype" = "ncrna" ]; then
    suffixes=("ncrna.fa.gz")
elif [ "$datatype" = "pep" ]; then
    suffixes=("pep.all.fa.gz")
else
    echo "ERROR:  dna、cdna、cds、ncrna or pep " 
    exit 1
fi

success=false
for suffix in "${suffixes[@]}"; do
    url_https="${url_prefix}.${suffix}"
    url_ftp="${url_https/https:\/\//ftp:\/\/}"
    
    # 检查文件是否存在
    if curl --location --head "$url_https" 2>/dev/null | grep -q 'Content-Length'; then
        lftp -c "pget -n ${threads} ${url_https} -o ${tmpdir}/${suffix}";
        $decompress ${tmpdir}/${suffix} >> "$output_file"
        success=true

    elif curl --location --head "$url_ftp" 2>/dev/null | grep -q 'Content-Length'; then
        lftp -c "pget -n ${threads} ${url_ftp} -o ${tmpdir}/${suffix} ";
        $decompress ${tmpdir}/${suffix} >> "$output_file"
        success=true
    fi
done

if ! $success; then
    echo "ERROR"
    exit 1
fi