#===================================================================
# Commands to process script output files for QIIME
#
# Files needed:
# - hitchip_data.txt
# - qiime_sample_mapping_file.txt
# - HITChip16S_tax_assignments.txt
#
# Output directory including HITChip data in QIIME format:
# - outputdir
#
#==================================================================

# Assign taxonomy to HITChip FASTA sequences. Needs to be done just once.
# ~/data/HITChip/QIIME/taxonomicAssignment.sh
# assign_taxonomy.py -i ~/data/HITChip/QIIME/HITChip16S.fasta -m rdp –o hitchip_taxonomic_assignments

# Convert HITChip species matrix into QIIME abundance matrix
~/bin/R-2.15.2/bin/R CMD BATCH HITChip2QIIME.R

# Check mapping file correctness
check_id_map.py -m qiime_sample_mapping_file.txt -o mapping_output

# Creat BIOM abundance table
trflp_file_to_otu_table.py -i hitchip_data.txt -o abundance_otu.table.biom

# Link observations and metadata
add_metadata.py -i abundance_otu.table.biom -o abundance_with_taxa.biom --sample_mapping_fp qiime_sample_mapping_file.txt --observation_mapping_fp hitchip_taxonomic_assignments/HITChip16S_tax_assignments.txt --observation_header OTUID,taxonomy,confidence --sc_separated taxonomy

# Get the different summary levels from QIIME for HITChip
rm -rf outputdir
summarize_taxa_through_plots.py -m qiime_sample_mapping_file.txt -i abundance_with_taxa.biom -o outputdir

