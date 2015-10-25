import gzip
import os
import sys
from mysql.connector import MySQLConnection

print("Connecting to db")
conn = MySQLConnection(user='root', host='exac-dev', port=3307, database='exac_readviz')
c = conn.cursor(buffered=True)


def print_query(q):
    print(q)


def run_query(q):
    print('---------')
    print(q)
    c.execute(q)
    if q.lower().strip().startswith("update"):
        conn.commit()
    print("%d rows" % c.rowcount)
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

reset_missing_bams_that_actually_exist = False
if reset_missing_bams_that_actually_exist:
    print("step 1: Find bams that caused a file-doesn't-exist error, but actually do exist")
    c.execute("select distinct original_bam_path from sample where hc_error_code=1000")
    all_original_bam_paths = list(c.fetchall())
    print("Found %d distinct bams that caused 1000 error" % len(all_original_bam_paths))
    found_bam_paths = tuple([p[0] for p in all_original_bam_paths if os.path.isfile(p[0])])


    print("Found %d such bams. Reset records with missing-bam errors to finished=0 for bams in this list" % len(found_bam_paths))
    run_query(("update variant as v join sample as s on "
              "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom=s.het_or_hom "
              "set v.finished=0, s.hc_error_code=s.hc_error_code+10 "
              "where v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples and s.hc_error_code=1000 and s.original_bam_path IN %s") % str(found_bam_paths))

    run_query(("update sample as s join variant as v on "
               "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom=s.het_or_hom "
               "set s.finished=0 where v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples and s.hc_error_code=1000 and s.original_bam_path IN %s") % str(found_bam_paths))


reset_samples_with_transient_errors = True
if reset_samples_with_transient_errors:

    print("For *samples* with transient errors, reset them to finished=0")
    run_query(("update sample as s join variant as v on "
               "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom=s.het_or_hom "
               "set s.finished=0, s.hc_error_code=NULL, s.hc_error_text=NULL "
               "where (s.hc_error_code IN (2001, 2011, 2009, 2019, 2021, 4000) or s.hc_error_text like '%ermiss%') "    # 3001,
               "and v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples"))

    print("For *variants* with transient errors, reset them to finished=0")
    run_query(("update variant as v join sample as s on "
               "v.chrom=s.chrom and v.pos=s.pos and v.ref=s.ref and v.alt=s.alt and v.het_or_hom=s.het_or_hom "
               "set v.finished=0 "
               "where (s.hc_error_code IN (2001, 2011, 2009, 2019, 2021, 4000) or s.hc_error_text like '%ermiss%') "    # 3001,
               "and v.n_available_samples>=0 and v.n_available_samples<v.n_expected_samples"))

ALL_CHROMS = list(map(str, [10,11,12,13,14,15,16,17,18,19,20,21,22, 1,2,3,4,5,6,7,8,9, 'X','Y']))
#FINISHED_CHROMOSOMES = str(tuple(map(str, ALL_CHROMS))).replace(",)", ")") #


# reset unfinished intervals not in top genes and not in genomic regions
reset_unfinished_intervals_in_important_genes = False
if reset_unfinished_intervals_in_important_genes:

    print("Reset intervals in important genes")
    sql_is_overlapping = []
    with open(os.path.join(os.path.abspath(os.path.dirname(__file__)),  "data/top_intervals.txt")) as f:
        for line in f:
            i_name, i_chrom, i_start, i_end = line.strip("\n").split("\t")
            sql_is_overlapping.append( "(chrom='%(i_chrom)s' and start_pos <= %(i_end)s and end_pos>=%(i_start)s )" % locals() )
            if len(sql_is_overlapping) > 300:
                break



    run_query(("select * from parallelize.python3_4_generate_HC_bams_py_i200 "
               "where finished=0 and ( %s )") % " or ".join(sql_is_overlapping))


    # enable intervals that are overlapping the intervals of interest
    print("Reset intervals overlapping intervals of interest")
    run_query(("update parallelize.python3_4_generate_HC_bams_py_i200 set job_id=NULL, task_id=NULL, unique_id=NULL "
               "where (job_id is NULL or job_id = -1 ) and finished=0 and ( %s )") % " or ".join(sql_is_overlapping))



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
# parallelize.python3_4_generate_HC_bams_py_i200
# get all intervals
from collections import defaultdict
all_intervals = defaultdict(set)
c = run_query("select chrom, start_pos, end_pos, started, finished, error_code from parallelize.python3_4_generate_HC_bams_py_i200") #" where chrom in %(FINISHED_CHROMOSOMES)s" % locals())
for chrom, start_pos, end_pos, started, finished, error_code in c:
    all_intervals[chrom].add((chrom, int(start_pos), int(end_pos), int(started), int(finished), int(error_code)))


total = 0
for chrom in all_intervals:
    total += len(all_intervals[chrom])
    print("%s: %s intervals" % (chrom, len(all_intervals[chrom])))
print("total: %s intervals" % total)




# reset intervals marked as finished
reset_intervals_that_contain_unfinished_variants = True
if reset_intervals_that_contain_unfinished_variants:
    for current_chrom in ALL_CHROMS:
        c = run_query("select chrom, pos from variant as v where chrom='%(current_chrom)s' and finished=0" % locals()) #" and chrom in %(FINISHED_CHROMOSOMES)s" % locals())
        all_uninished_variants = c.fetchall()

        unfinished_intervals_marked_as_finished = set()
        current_interval = None
        print("Checking for unfinished_intervals_marked_as_finished using %s unfinished variants in chr%s" % (len(all_uninished_variants), current_chrom))

        for chrom, pos in all_uninished_variants:
            pos = int(pos)
            if current_interval is None or current_interval[0] != chrom or pos < current_interval[1] or pos > current_interval[2]:
                for i in all_intervals[chrom]:
                    if i[1] <= pos and pos <= i[2]:
                        current_interval = i
                        if i[4] > 0:
                            unfinished_intervals_marked_as_finished.add(i)
                            #print("%s: %s" % (len(unfinished_intervals_marked_as_finished), i))
                        #print("%(chrom)s %(pos)s is in interval %(i)s" % locals())
                        break
                else:
                    raise ValueError("%(chrom)s-%(pos)s is not in any intervals" % locals())
            #else:
                #print("%(chrom)s %(pos)s is in same interval %(current_interval)s" % locals())

        print("Found %s unfinished intervals marked as finished=1 in %s" % (len(unfinished_intervals_marked_as_finished), current_chrom))

        print("Updating intervals")
        for i in unfinished_intervals_marked_as_finished:
            chrom = i[0]
            start_pos = i[1]
            end_pos = i[2]
            run_query(("update parallelize.python3_4_generate_HC_bams_py_i200 set finished=0 "
                       "where chrom='%(chrom)s' and start_pos=%(start_pos)s and end_pos=%(end_pos)s") % locals())

        #print_query("update parallelize.python3_4_generate_HC_bams_py_i200 set "
        #            "job_id=NULL, task_id=NULL, unique_id=NULL, started=0, "
        #            "started_date=NULL, finished=0, finished_date=NULL, "
        #            "error_code=500, error_message=NULL where finished=0 and started_date <")

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

