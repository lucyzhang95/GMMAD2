# dump.sh: download all URLs passed as arguments
set -euo pipefail

urls=(
  "http://gepa.org.cn/static2/file/micro_metabolic.csv"
  "http://gepa.org.cn/static2/file/disease_species.csv"
  "http://gepa.org.cn/static2/file/meta_gene_net.csv"
)

for url in "${urls[@]}"; do
  echo "Downloading $url …"
  wget --continue --timestamping --progress=bar:force "$url"

  echo
done

echo "All GMMAD2 downloads are complete."

# Usage:
# chmod +x dump.sh
# ./dump.sh