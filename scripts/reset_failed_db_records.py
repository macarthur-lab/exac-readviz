"""
This script updates database records  in the interval, variant, sample tables to
recover from various kinds of transient errors fixable errors
such as cluster filesystem problems, some previously-generated files being deleted, etc.
"""
import os
import sys
from mysql.connector import MySQLConnection
from utils.constants import DB_HOST, DB_PORT, DB_USER

INTERVALS_TABLE = 'python2_pipeline_i1000000'

# initialize flags that control which sections are actually executed
reset_samples_with_transient_error = 0
reset_variants_with_transient_errors = 0
reset_variants_with_fewer_than_expected_available_samples = 0
reset_variants_with_original_bams_marked_missing_due_to_transient_error = 0
reset_variants_with_bams_in_db_but_not_on_disk = 0
reset_intervals_that_had_error_code = 0
reset_variants_that_contain_unfinished_samples = 0
reset_intervals_that_contain_unfinished_variants = 0
reset_intervals_that_contain_unfinished_samples = 0
reset_unfinished_intervals_to_clear_job_id = 0
reset_unfinished_samples_in_finished_chroms = 0
run_stat_queries = 0
reset_unfinished_samples_in_finished_chroms = 0
set_intervals_where_all_contained_variants_have_finished = 0
reset_unfinished_intervals_in_important_genes = 0
reset_samples_with_original_bams_marked_missing_due_to_transient_error = 0

# set flags to execute particular sections of code
#reset_variants_with_transient_errors = 1
#reset_variants_with_fewer_than_expected_available_samples = 1
#reset_variants_with_original_bams_marked_missing_due_to_transient_error = 1
#reset_variants_with_bams_in_db_but_not_on_disk = 1
#reset_variants_that_contain_unfinished_samples = 1
#reset_intervals_that_contain_unfinished_variants = 1
#reset_intervals_that_had_error_code = 1

reset_samples_with_transient_error = 1
reset_samples_with_original_bams_marked_missing_due_to_transient_error = 1
reset_unfinished_samples_in_finished_chroms = 1
reset_intervals_that_contain_unfinished_samples = 1

#reset_unfinished_intervals_to_clear_job_id = 1
#run_stat_queries = 1


print("connecting to db")
conn = MySQLConnection(user=DB_USER, host=DB_HOST, port=DB_PORT, database='exac_readviz')
c = conn.cursor(buffered=True)

def print_query(q):
    print(q)


def run_query(q):
    print('---------')
    print(q)
    #if not q.lower().strip().startswith("select"):
    #    return

    c.execute(q)
    if not q.lower().strip().startswith("select"):
        conn.commit()
        print("%d rows updated" % c.rowcount)
    else:
        print("%d rows returned" % c.rowcount)
    return c

"""
+----------------------------+--------------+------+-----+---------+----------------+
| Field                      | Type         | Null | Key | Default | Extra          |
+----------------------------+--------------+------+-----+---------+----------------+
| id                         | int(11)      | NO   | PRI | NULL    | auto_increment |
| chrom                      | varchar(5)   | NO   | MUL | NULL    |                |
| pos                        | int(11)      | NO   |     | NULL    |                |
| ref                        | varchar(500) | NO   |     | NULL    |                |
| alt                        | varchar(500) | NO   |     | NULL    |                |
| het_or_hom                 | varchar(3)   | NO   |     | NULL    |                |
| finished                   | tinyint(1)   | NO   |     | NULL    |                |
| sample_id                  | varchar(50)  | NO   |     | NULL    |                |
| sample_i                   | int(11)      | YES  |     | NULL    |                |
| original_bam_path          | longtext     | YES  |     | NULL    |                |
| original_gvcf_path         | longtext     | YES  |     | NULL    |                |
| output_bam_path            | longtext     | YES  |     | NULL    |                |
| is_missing_original_gvcf   | tinyint(1)   | NO   |     | NULL    |                |
| calling_interval_start     | int(11)      | YES  |     | NULL    |                |
| calling_interval_end       | int(11)      | YES  |     | NULL    |                |
| hc_succeeded               | tinyint(1)   | NO   |     | NULL    |                |
| hc_error_code              | int(11)      | YES  |     | NULL    |                |
| hc_error_text              | longtext     | YES  |     | NULL    |                |
| hc_n_artificial_haplotypes | int(11)      | YES  |     | NULL    |                |
| hc_started_time            | datetime     | YES  |     | NULL    |                |
| hc_finished_time           | datetime     | YES  |     | NULL    |                |
| hc_command_line            | longtext     | YES  |     | NULL    |                |
+----------------------------+--------------+------+-----+---------+----------------+
"""


if reset_variants_with_original_bams_marked_missing_due_to_transient_error:
    print("=== reset_variants_with_original_bams_marked_missing_due_to_transient_error ===")
    print("step 1: Find bams that caused a file-doesn't-exist error, but actually do exist")
    c.execute("select distinct original_bam_path from sample where hc_error_code=1000")
    all_original_bam_paths = list(c.fetchall())
    print("Found %d distinct bams that caused 1000 error" % len(all_original_bam_paths))
    found_bam_paths = tuple([p[0] for p in all_original_bam_paths if os.path.isfile(p[0])])
    print("Of these, %d actually exist on disk. Reset records with missing-bam errors to finished=0 for bams in this list" % len(found_bam_paths))
    if found_bam_paths:
        pass
        #run_query(("update variant as v join sample as s on "
        #      "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom_or_hemi=s.het_or_hom_or_hemi "
        #      "set v.finished=0, v.comments=NULL, n_available_samples=NULL, n_expected_samples=NULL, readviz_bam_paths=NULL "
        #      "where v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples and "
        #       "s.hc_error_code IN (1000, 1010) and s.original_bam_path IN %s") % str(found_bam_paths).replace(',)', ')'))

        #run_query(("update sample as s join variant as v on "
        #       "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom_or_hemi=s.het_or_hom_or_hemi "
        #       "set s.finished=0, s.comments=NULL, hc_succeeded=0, hc_error_code=NULL, hc_error_text=NULL, sample_i=NULL, original_bam_path=NULL, original_gvcf_path=NULL, output_bam_path=NULL, hc_command_line=NULL "
        #       "where v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples and "
        #       "s.hc_error_code IN (1000, 1010) and s.original_bam_path IN %s") % str(found_bam_paths).replace(',)', ')'))

if reset_samples_with_original_bams_marked_missing_due_to_transient_error:
    print("=== reset_samples_with_original_bams_marked_missing_due_to_transient_error ===")
    print("step 1: Find bams that caused a file-doesn't-exist error, but actually do exist")
    c.execute("select distinct original_bam_path from sample where hc_error_code=1000")
    all_original_bam_paths = list(c.fetchall())
    print("Found %d distinct bams that caused 1000 error" % len(all_original_bam_paths))
    found_bam_paths = tuple([p[0] for p in all_original_bam_paths if os.path.isfile(p[0])])
    print("Of these, %d actually exist on disk. Reset records with missing-bam errors to finished=0 for bams in this list" % len(found_bam_paths))
    if found_bam_paths:
        run_query("update sample set started=0, started_time=NULL, finished=0, finished_time=NULL, hc_succeeded=0, hc_error_text=NULL, hc_error_code=NULL, comments=NULL "
                  "where hc_error_code=1000 and original_bam_path IN %s" % str(found_bam_paths).replace(',)', ')').replace("u'", "'"))

if reset_variants_with_transient_errors:
    print("=== reset_variants_with_transient_errors ===")
    print("For *samples* with transient errors, reset them to finished=0")
    run_query(("update sample as s join variant as v on "
               "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom_or_hemi=s.het_or_hom_or_hemi "
               "set s.comments=NULL, s.finished=0, hc_succeeded=0, hc_error_code=NULL, hc_error_text=NULL, sample_i=NULL, original_bam_path=NULL, original_gvcf_path=NULL, output_bam_path=NULL, hc_command_line=NULL "
               "where (s.hc_error_code IN (2001, 2011, 2009, 2019, 2021, 4000) or (s.hc_error_code is NULL and s.hc_succeeded=0)) "    # 3001,
               "and v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples"))

    print("For *variants* with transient errors, reset them to finished=0")
    run_query(("update variant as v join sample as s on "
               "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom_or_hemi=s.het_or_hom_or_hemi "
               "set v.comments=NULL, v.finished=0, n_available_samples=NULL, n_expected_samples=NULL, readviz_bam_paths=NULL  "
               "where (s.hc_error_code IN (2001, 2011, 2009, 2019, 2021, 4000) or (s.hc_error_code is NULL and s.hc_succeeded=0)) "    # 3001,
               "and v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples"))

ALL_CHROMS = list(map(str, [1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9, 'X','Y']))
FINISHED_CHROMS = list(map(str, [12,14,15,7,"X"])) #21,22,3,4,5,6,7,8,9, 'X','Y']))
FINISHED_CHROMS_STRING = str(tuple(map(str, FINISHED_CHROMS))).replace(",)", ")") #


# reset unfinished intervals not in top genes and not in genomic regions
if reset_unfinished_intervals_in_important_genes:
    print("=== reset_unfinished_intervals_in_important_genes ===")
    print("Reset intervals in important genes")
    sql_is_overlapping = []
    with open(os.path.join(os.path.abspath(os.path.dirname(__file__)),  "data/top_intervals.txt")) as f:
        for line in f:
            i_name, i_chrom, i_start, i_end = line.strip("\n").split("\t")
            sql_is_overlapping.append( "(chrom='%(i_chrom)s' and start_pos <= %(i_end)s and end_pos>=%(i_start)s )" % locals() )
            if len(sql_is_overlapping) > 300:
                break

    run_query(("select * from " + INTERVALS_TABLE +
               "where finished=0 and ( %s )") % " or ".join(sql_is_overlapping))

    # enable intervals that are overlapping the intervals of interest
    print("Reset intervals overlapping intervals of interest")
    run_query(("update " + INTERVALS_TABLE + " set job_id=NULL, comments=NULL, task_id=NULL, unique_id=NULL "
               "where (job_id is NULL or job_id = -1 ) and finished=0 and ( %s )") % " or ".join(sql_is_overlapping))



if run_stat_queries or set_intervals_where_all_contained_variants_have_finished or reset_intervals_that_contain_unfinished_variants or reset_intervals_that_contain_unfinished_samples:
    """
    +---------------------+--------------+------+-----+---------+----------------+
    | Field               | Type         | Null | Key | Default | Extra          |
    +---------------------+--------------+------+-----+---------+----------------+
    | id                  | int(11)      | NO   | PRI | NULL    | auto_increment |
    | chrom               | varchar(5)   | NO   | MUL | NULL    |                |
    | pos                 | int(11)      | NO   |     | NULL    |                |
    | ref                 | varchar(500) | NO   |     | NULL    |                |
    | alt                 | varchar(500) | NO   |     | NULL    |                |
    | het_or_hom          | varchar(3)   | NO   |     | NULL    |                |
    | finished            | tinyint(1)   | NO   |     | NULL    |                |
    | n_expected_samples  | int(11)      | YES  |     | NULL    |                |
    | n_available_samples | int(11)      | YES  |     | NULL    |                |
    | readviz_bam_paths   | longtext     | YES  |     | NULL    |                |
    +---------------------+--------------+------+-----+---------+----------------+
    """

    # python3_4_generate_HC_bams_py_i200
    # get all intervals
    from collections import defaultdict
    all_intervals = defaultdict(set)
    c = run_query("select chrom, start_pos, end_pos, started, finished, error_code from %(INTERVALS_TABLE)s" % locals())
    for chrom, start_pos, end_pos, started, finished, error_code in c:
        all_intervals[chrom].add((chrom, int(start_pos), int(end_pos), int(started), int(finished), int(error_code)))

    all_finished_intervals = defaultdict(set)
    c = run_query("select chrom, start_pos, end_pos, started, finished, error_code from %(INTERVALS_TABLE)s where started=1 and finished=1" % locals())
    for chrom, start_pos, end_pos, started, finished, error_code in c:
        all_finished_intervals[chrom].add((chrom, int(start_pos), int(end_pos), int(started), int(finished), int(error_code)))

    total = finished_count = 0
    for chrom in all_intervals:
        finished_count += len(all_finished_intervals[chrom])
        total += len(all_intervals[chrom])
        print("%s: %s finished out of %s intervals" % (chrom, len(all_finished_intervals[chrom]), len(all_intervals[chrom])))
    print("total: %s finished out of %s total intervals" % (finished_count, total))



# set intervals as finished if all variants in them are finished
if set_intervals_where_all_contained_variants_have_finished:
    # unlike the other recovery code
    # this is useful for rebuilding the intervals table from scratch by setting
    # each interval record to finished if all variants it contains are finished
    print("=== set_intervals_where_all_contained_variants_have_finished ===")
    for current_chrom in ALL_CHROMS:
        #all_intervals[chrom].add((chrom, int(start_pos), int(end_pos), int(started), int(finished), int(error_code)))
        for _, start_pos, end_pos, _, _, _ in all_intervals[current_chrom]:
            c = run_query("select chrom, pos from variant as v where chrom='%(current_chrom)s' and pos >= %(start_pos)s and pos <= %(end_pos)s and finished=0" % locals())
            if c.rowcount > 0:
                print("Found %s unfinished variants in %s:%s-%s. Skipping" % (c.rowcount, current_chrom, start_pos, end_pos))
            else:
                run_query(("update " + INTERVALS_TABLE + " set job_id=1, started=1, finished=1 "
                           "where chrom='%(current_chrom)s' and start_pos=%(start_pos)s and end_pos=%(end_pos)s") % locals())

if reset_variants_with_fewer_than_expected_available_samples:
    print("=== reset_variants_with_fewer_than_expected_available_samples ===")
    # Reset samples to finished = 0 where hc_u
    run_query(("update sample as s join variant as v on "
               "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom_or_hemi=s.het_or_hom_or_hemi "
               "set s.finished=0, s.comments=NULL, hc_succeeded=0, hc_error_code=NULL, hc_error_text=NULL, sample_i=NULL, original_bam_path=NULL, original_gvcf_path=NULL, output_bam_path=NULL, hc_command_line=NULL "
               "where v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples and "
               "v.finished=0 and s.finished=1 and s.hc_succeeded=0"))

    # Reset variants to finished = 0 where number of records in the sample table that either succeeded or had an error is < n_expected_samples

    run_query(("update variant as v "
              "set v.finished=0, v.comments=NULL, n_available_samples=NULL, n_expected_samples=NULL, readviz_bam_paths=NULL  "
              "where n_available_samples<n_expected_samples and "
              "n_expected_samples > ("
              " select count(*) from sample as s where chrom=v.chrom and pos=v.pos and alt=v.alt and het_or_hom_or_hemi=v.het_or_hom_or_hemi and (hc_succeeded=1 or hc_error_code>0)"
              ")"))

    # Reset variants to finished = 0 where the variant.n_available_samples < records in the sample table that have hc_succeeded=1
    run_query(("update variant as v "
              "set v.finished=0, v.comments=NULL, n_available_samples=NULL, n_expected_samples=NULL, readviz_bam_paths=NULL "
              "where n_available_samples<n_expected_samples and n_available_samples < ("
              "  select count(*) from sample as s where chrom=v.chrom and pos=v.pos and ref=v.ref and alt=v.alt and het_or_hom_or_hemi=v.het_or_hom_or_hemi and hc_succeeded=1"
              ")"))

if reset_variants_with_bams_in_db_but_not_on_disk:
    print("=== reset_variants_with_bams_in_db_but_not_on_disk ===")
    import glob
    os.chdir("/broad/hptmp/exac_readviz_backend/")
    for current_chrom in ALL_CHROMS:
        print("globbing for all bam files in chr%s" % current_chrom)
        actual_files_on_disk = set(glob.glob(current_chrom +"/*/chr*.bam"))
        for t in run_query("select readviz_bam_paths, chrom, pos, ref, alt, het_or_hom_or_hemi from variant "
                           "where chrom='%(current_chrom)s' and finished=1 and n_available_samples>0" % locals()).fetchall():
            cached_filenames_list = t[0].split('|')
            cached_filenames_set = set(cached_filenames_list)
            duplicates_found = len(cached_filenames_set) < len(cached_filenames_list)
            some_cached_files_not_found_on_disk = len(cached_filenames_set - actual_files_on_disk) > 0
            if duplicates_found or some_cached_files_not_found_on_disk:
                if some_cached_files_not_found_on_disk:
                    print('readviz_bam_paths that are no longer found on disk: %s ' % str(cached_filenames_set - actual_files_on_disk))
                else:
                    print('duplicates_found: %s' % str(cached_filenames_list))
                run_query("update variant as v "
                          "set v.finished=0, v.comments=NULL, n_available_samples=NULL, n_expected_samples=NULL, readviz_bam_paths=NULL "
                          "where chrom='%s' and pos=%s and ref='%s' and alt='%s' and het_or_hom_or_hemi='%s' " % t[1:])
                run_query("update sample set started=0, started_time=NULL, finished=0, finished_time=NULL, hc_succeeded=0, hc_error_text=NULL, hc_error_code=NULL, comments=NULL "
                          "where chrom='%s' and pos=%s and ref='%s' and alt='%s' and het_or_hom_or_hemi='%s' " % t[1:])

if reset_variants_that_contain_unfinished_samples:
    print("=== reset_variants_that_contain_unfinished_samples ===")
    run_query(("update sample as s join variant as v on "
               "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom_or_hemi=s.het_or_hom_or_hemi "
               "set s.finished=0, s.comments=NULL, hc_succeeded=0, hc_error_code=NULL, hc_error_text=NULL, sample_i=NULL, original_bam_path=NULL, original_gvcf_path=NULL, output_bam_path=NULL, hc_command_line=NULL "
               "where v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples and "
               "s.finished=0"))

    run_query(("update variant as v join sample as s on "
               "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom_or_hemi=s.het_or_hom_or_hemi "
               "set v.finished=0, v.comments=NULL, n_available_samples=NULL, n_expected_samples=NULL, readviz_bam_paths=NULL "
               "where s.finished=0"))
               
if reset_intervals_that_contain_unfinished_variants:
    print("=== reset_intervals_that_contain_unfinished_variants ===")
    for current_chrom in FINISHED_CHROMS:
        c = run_query("select chrom, pos from variant as v where chrom='%(current_chrom)s' and v.finished=0 order by pos asc" % locals()) 
        all_unfinished_variants = c.fetchall()

        unfinished_intervals = set()
        current_interval = None
        print("Checking for unfinished_intervals using %s unfinished variants in chr%s" % (len(all_unfinished_variants), current_chrom))

        for chrom, pos in all_unfinished_variants:
            pos = int(pos)
            if current_interval is None or current_interval[0] != chrom or pos < current_interval[1] or pos > current_interval[2]:
                for i in all_intervals[chrom]:
                    if i[1] <= pos and pos <= i[2]:
                        current_interval = i
                        sys.stdout.write("Found matching interval %s for variant: %s" % (i, "%s:%s" % (chrom, pos)))
                        if i[3] > 0 or i[4] > 0:
                            unfinished_intervals.add(i)
                            sys.stdout.write(". Will reset it..\n")
                        else:
                            sys.stdout.write(". It's already marked as not started.\n")
                            #print("%s: %s" % (len(unfinished_intervals), i))
                        #print("%(chrom)s %(pos)s is in interval %(i)s" % locals())
                        break
                else:
                    raise ValueError("%(chrom)s-%(pos)s is not in any intervals" % locals())

            #else:
                #print("%(chrom)s %(pos)s is in same interval %(current_interval)s" % locals())

        print("Found %s unfinished intervals in %s" % (len(unfinished_intervals), current_chrom))

        print("Updating intervals")
        for i in unfinished_intervals:
            chrom = i[0]
            start_pos = i[1]
            end_pos = i[2]
            run_query(("update " + INTERVALS_TABLE +
                       "set job_id=NULL, comments=NULL, unique_id=NULL, task_id=NULL, started=0, started_date=NULL, error_message=NULL, error_code=0, finished=0, error_message=NULL "
                       "where chrom='%(chrom)s' and start_pos=%(start_pos)s and end_pos=%(end_pos)s") % locals())

        #print_query("update python3_4_generate_HC_bams_py_i200 set "
        #            "job_id=NULL, task_id=NULL, unique_id=NULL, started=0, "
        #            "started_date=NULL, finished=0, finished_date=NULL, "
        #            "error_code=500, error_message=NULL where finished=0 and started_date <")

if reset_samples_with_transient_error:
        run_query(("update sample set started=0, started_time=NULL, finished=0, finished_time=NULL, hc_succeeded=0, hc_error_text=NULL, hc_error_code=NULL, comments=NULL "
                   "where hc_error_code >= 2000 and hc_error_code < 3000 and chrom in %(FINISHED_CHROMS_STRING)s") % locals())

if reset_unfinished_samples_in_finished_chroms:
    run_query("update sample set started=0, started_time=NULL, finished=0, finished_time=NULL, hc_succeeded=0, hc_error_text=NULL, hc_error_code=NULL, comments=NULL "
              "where chrom in %(FINISHED_CHROMS_STRING)s and started in (0, 1) and finished=0" % locals())


if reset_intervals_that_contain_unfinished_samples:
    print("=== reset_intervals_that_contain_unfinished_samples ===")
    for current_chrom in FINISHED_CHROMS:
        c = run_query("select chrom, pos from sample as s where chrom='%(current_chrom)s' and s.started in (0, 1) and s.finished=0 order by pos asc" % locals())
        all_unfinished_samples = c.fetchall()

        unfinished_intervals = set()
        current_interval = None
        print("Checking for unfinished_intervals using %s unfinished samples in chr%s" % (len(all_unfinished_samples), current_chrom))

        for chrom, pos in all_unfinished_samples:
            pos = int(pos)
            if current_interval is None or current_interval[0] != chrom or pos < current_interval[1] or pos > current_interval[2]:
                for i in all_intervals[chrom]:
                    if i[1] <= pos and pos <= i[2]:
                        current_interval = i
                        sys.stdout.write("Found matching interval %s for variant: %s" % (i, "%s:%s" % (chrom, pos)))
                        if not (i[3] == 0 and i[4] == 0):  # reset intervals that are not started or not finished
                            unfinished_intervals.add(i)
                            sys.stdout.write(". Will reset it..\n")
                        else:
                            sys.stdout.write(". It's already marked as not started.\n")

                        break
                else:
                    raise ValueError("%(chrom)s-%(pos)s is not in any intervals" % locals())

            #else:
                #print("%(chrom)s %(pos)s is in same interval %(current_interval)s" % locals())

        print("Found %s unfinished intervals in %s" % (len(unfinished_intervals), current_chrom))

        print("Updating intervals")
        for i in unfinished_intervals:
            chrom = i[0]
            start_pos = i[1]
            end_pos = i[2]
            run_query(("update %(INTERVALS_TABLE)s "
                       "set job_id=null, task_id=null, unique_id=null, started=0, started_date=null, finished=0, finished_date=null, "
                       "error_code=0, error_message=null, priority=null, username=null, machine_hostname=null, machine_average_load=null, comments=null " 
                       "where chrom='%(chrom)s' and start_pos=%(start_pos)s and end_pos=%(end_pos)s") % locals())

        #print_query("update python3_4_generate_HC_bams_py_i200 set "
        #            "job_id=NULL, task_id=NULL, unique_id=NULL, started=0, "
        #            "started_date=NULL, finished=0, finished_date=NULL, "
        #            "error_code=500, error_message=NULL where finished=0 and started_date <")


if reset_intervals_that_had_error_code:
    print("=== reset_intervals_that_had_error_code ===")
    run_query("update " + INTERVALS_TABLE +
              "set job_id=null, task_id=null, unique_id=null, started=0, started_date=null, finished=0, finished_date=null, "
              "error_code=0, error_message=null, priority=null, username=null, machine_hostname=null, machine_average_load=null, comments=null " 
              "where error_code > 0")

if reset_unfinished_intervals_to_clear_job_id:
    print("=== reset_unfinished_intervals_to_clear_job_id ===")
    run_query("update " + INTERVALS_TABLE +
              "set job_id=null, task_id=null, unique_id=null, started=0, started_date=null, finished=0, finished_date=null, "
              "error_code=0, error_message=null, priority=null, username=null, machine_hostname=null, machine_average_load=null, comments=null " 
              "where finished=0")


print("Done")


"""
+----------------------+--------------+------+-----+---------+----------------+
| Field                | Type         | Null | Key | Default | Extra          |
+----------------------+--------------+------+-----+---------+----------------+
| id                   | int(11)      | NO   | PRI | NULL    | auto_increment |
| chrom                | varchar(5)   | NO   | MUL | NULL    |                |
| start_pos            | int(11)      | NO   |     | NULL    |                |
| end_pos              | int(11)      | NO   |     | NULL    |                |
| job_id               | int(11)      | YES  | MUL | NULL    |                |
| task_id              | int(11)      | YES  |     | NULL    |                |
| unique_id            | int(11)      | YES  |     | NULL    |                |
| started              | tinyint(1)   | NO   | MUL | NULL    |                |
| started_date         | datetime     | YES  | MUL | NULL    |                |
| finished             | tinyint(1)   | NO   | MUL | NULL    |                |
| finished_date        | datetime     | YES  |     | NULL    |                |
| error_code           | int(11)      | NO   | MUL | NULL    |                |
| error_message        | longtext     | YES  |     | NULL    |                |
| machine_hostname     | varchar(100) | YES  |     | NULL    |                |
| machine_average_load | float        | YES  |     | NULL    |                |
+----------------------+--------------+------+-----+---------+----------------+
"""

c.close()
conn.close()

