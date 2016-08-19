update variant set started=0, started_time=NULL, finished=0, finished_time=NULL, comments=NULL, username=NULL, n_available_samples=0, readviz_bam_paths=NULL;
# drop table sample;

#  reset all sample records 
# update sample set started=0, started_time=NULL, finished=0, finished_time=NULL, original_bam_path=NULL, original_gvcf_path=NULL, output_bam_path=NULL, hc_succeeded=0, hc_error_text=NULL, hc_error_code=NULL, hc_command_line=NULL, sample_id='', comments=NULL;

