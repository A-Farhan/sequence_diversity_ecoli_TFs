# given a project accession, download metadata, filter and transfer to server for read files downloading

# command line argument
prj=$1

echo "creating directory $prj, if it doesn't exist already"
if [[ ! -d $prj ]]; then
    mkdir $prj
fi

# paths
fl="$prj/${prj}.txt"

echo "downloading metadata table for the project $prj"
wget "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$prj&result=read_run&download=txt&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,library_layout,library_strategy,library_source,library_selection,base_count,first_public,last_updated,fastq_ftp,sample_alias" -O $fl

echo "filtering data based on relevance and cut-off"
python filter_data.py $fl 50

echo "transferring metadata to the storage server to download read sequences"
scp -r $prj <storage_server_address>
