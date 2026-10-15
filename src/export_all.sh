mkdir -p exported

for f in *.qza; do
    base=$(basename "$f" .qza)
    qiime tools export \
        --input-path "$f" \
        --output-path "exported/$base"
done
