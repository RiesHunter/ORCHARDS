##################################################
## Concatenate all HA .fa files
for file in *.fa
do
  # Remove the .fa extension from the filename to get the base name
  f=${file%%.fa}

  ## Extract HA sequences and modify names
  # Extract lines containing ">HA" and save to a temporary file
  sed -n '/>*HA/p' ./${f}.fa > name
  
  # Replace "HA" with "HA|<filename>" in the temporary name file
  sed "s/HA/HA|${f}/g" name > name.temp
  
  # Append modified names to the concatenated file
  cat name.temp >> segmented_compiled-HA
  
  ## Extract sequences
  # Extract sequences between the header ">HA" and the next header
  sed -n '/>*HA/,/>/p' ./${f}.fa > seq
  
  # Remove the last line (which is likely a header) from the seq file
  sed '$ d' seq > seq.temp
  
  # Remove the first line (the header) from seq.temp to get only the sequence
  sed '1d' seq.temp > seq
  
  # Append the extracted sequence to the concatenated file
  cat seq >> segmented_compiled-HA
  
  # Show progress in the terminal
  echo -ne "${f}\r"
  echo ""
done

# Rename the concatenated file
mv segmented_compiled-HA segmented_compiled-HA.fasta

# Clean up temporary files
rm seq
rm name
rm name.temp
rm seq.temp

### Remove everything before HA in the concatenated file
# Replace specific headers and identifiers in the final file
# H3
sed -i.bak 's|>A/Singapore/INFIMH-16-0019/2016|>|g' segmented_compiled-HA.fasta
sed -i.bak 's|EPI_ISL_296168||g' segmented_compiled-HA.fasta
sed -i.bak 's|>A/Hong_Kong/4801/2014|>|g' segmented_compiled-HA.fasta
sed -i.bak 's|EPI_ISL_233740||g' segmented_compiled-HA.fasta
# H1
sed -i.bak 's|>A/Michigan/45/2015|>|g' segmented_compiled-HA.fasta
sed -i.bak 's|EPI_ISL_336680||g' segmented_compiled-HA.fasta
# Remove "HA" and other identifiers
sed -i.bak 's|HA||g' segmented_compiled-HA.fasta
sed -i.bak 's|_consensus||g' segmented_compiled-HA.fasta
sed -i.bak 's/|||//g' segmented_compiled-HA.fasta

# Remove the backup file created by sed
rm segmented_compiled-HA.fasta.bak

### Align concatenated HA .fa file
# Run the alignment using ClustalW
clustalw \
  -ALIGN \
  -INFILE=./segmented_compiled-HA.fasta \
  -OUTFILE=./clustalw_HA.fasta \
  -OUTPUT=FASTA

#### RAxML-NG for Maximum Likelihood estimation
# Create .phy file to check for MSA
raxml-ng --check --msa ./clustalw_*.fasta --model GTR+G
# Parse the alignment and format if necessary
raxml-ng --parse --msa ./clustalw_*.fasta --model GTR+G
# Infer the maximum likelihood tree
raxml-ng --msa ./clustalw_*.fasta --model GTR+G --prefix clustalw --threads 8 --seed 920
##Generate bootstrap replicates
#raxml-ng --bootstrap --msa ./clustalw_*.fasta --model GTR+G --prefix clustalw --seed 920 --bs-trees 2000 --threads 8
##Perform tree search using the bootstrap replicates
#raxml-ng --bsconverge --bs-trees clustalw.raxml.bootstraps --prefix clustalw --seed 920 --threads 8
##Combine the bootstrap trees with the maximum likelihood tree
#raxml-ng --msa ./clustalw_*.fasta --model GTR --prefix clustalw --threads 8 --seed 920 --tree clustalw.raxml.bestTree --bs-trees clustalw.raxml.bootstraps --support


#### Clean up
# Get the second-to-last directory part (e.g., "18-19_H1N1")
second_last_dir=$(basename "$(dirname "$PWD")")

# Define the target directory (using quotes to handle spaces)
target_dir="/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs/dfs/"

# Ensure the target directory exists
if [ ! -d "$target_dir" ]; then
  echo "Target directory does not exist. Creating directory: $target_dir"
  mkdir -p "$target_dir"
fi

# Construct the new filename with the second-to-last directory part
new_file_name="${second_last_dir}-clustalw.raxml.bestTree"

# Copy and rename the file in one step
cp "clustalw.raxml.bestTree" "${target_dir}${new_file_name}"

# Confirm the file has been copied
if [ $? -eq 0 ]; then
  echo "File successfully copied and renamed to ${new_file_name} in ${target_dir}"
else
  echo "Failed to copy the file."
fi

# Construct the new filename for clustalw_HA.fasta with the second-to-last directory part
new_ha_file_name="${second_last_dir}-clustalw_HA.fasta"

# Copy and rename the clustalw_HA.fasta file in one step
cp "clustalw_HA.fasta" "${target_dir}${new_ha_file_name}"

# Confirm the file has been copied
if [ $? -eq 0 ]; then
  echo "File successfully copied and renamed to ${new_ha_file_name} in ${target_dir}"
else
  echo "Failed to copy the file."
fi

# Define the ML directory (create it if it doesn't exist)
ml_dir="ML"
if [ ! -d "$ml_dir" ]; then
  echo "Creating directory: $ml_dir"
  mkdir "$ml_dir"
fi

# List of files to move
files_to_move=(
  "clustalw.raxml.bestModel"
  "clustalw.raxml.bestTree"
  "clustalw.raxml.bestTreeCollapsed"
  "clustalw.raxml.log"
  "clustalw.raxml.mlTrees"
  "clustalw.raxml.rba"
  "clustalw.raxml.reduced.phy"
  "clustalw.raxml.startTree"
  "clustalw_HA.fasta"
  "clustalw_HA.fasta.raxml.log"
  "clustalw_HA.fasta.raxml.rba"
  "clustalw_HA.fasta.raxml.reduced.phy"
  "segmented_compiled-HA.dnd"
  "segmented_compiled-HA.fasta"
)

# Move each file to the ML directory if it exists
for file in "${files_to_move[@]}"; do
  if [ -f "$file" ]; then
    echo "Moving $file to $ml_dir/"
    mv "$file" "$ml_dir/"
  else
    echo "$file does not exist in the current directory."
  fi
done

# Confirm all files have been moved
echo "All files moved to $ml_dir."
