"""
This script takes the name of another python script (let's say some_script.py)
and launches N parallel instances of some_script.py on the short queue as an
array job. Each instance of some_script.py will get passed
 --chrom, --start, and --end args which will define the genomic
region it should operate on. This region should be small enough for
some_script.py to finish in < 4 hours (the short queue's runtime limit) in the
worst case. It's up to some_script.py to avoid redoing the same work if
it is run multiple times on the same genomic interval (eg. if the 1st run fails).

some_script.py should return 0 if it completed successfully, or not 0 otherwise.
"""

import argparse
import datetime
import os
import peewee
import random


# table for keeping track of which intervals some_script.py has already processed
class ParallelIntervals(peewee.Model):
    chrom = peewee.CharField(max_length=5, null=False)
    start_pos = peewee.IntegerField(null=False)
    end_pos = peewee.IntegerField(null=False)
    #current_pos = peewee.IntegerField(default=0)  # not currently used because some_script.py has no way of feeding back this info

    job_id = peewee.IntegerField(default=0)  # job that is processing this interval
    finished_successfully = peewee.BooleanField(default=False)
    finished_with_error = peewee.IntegerField(default=0)

    error_message = peewee.TextField()

    started_date = peewee.DateTimeField(default=datetime.datetime.now)
    finished_date = peewee.DateTimeField(default=None)

    # execution environment stats
    machine_hostname = peewee.CharField(max_length=30)
    machine_average_load = peewee.Float(default=0)

    class Meta:
        database = db
        db_table = "parallelintervals_%s" % ( os.path.basename(script_name), )


indexes = (
    (('chrom', 'start_pos', 'end_pos'), True), # True means unique index
    (('started', 'finished_successfully', 'finished_with_error'), False),
)




# if this instance of parallelize.py is already running as one of many array
# job tasks, keep launching some_script.py on the next unprocessed interval
is_array_job_task = os.getenv('SGE_TASK_ID', False)  # the actual task id number doesn't matter
if is_array_job_task:
    unique_job_id = random.random()  # don't use actual job id to avoid collisions in case parallel.py was restarted
    hostname = os.getenv('HOSTNAME')
    while True:
        stared_date = datetime.datetime.now()

        # get next interval
        with db.atomic() as txn:
            # claim an interval
            interval_found = ParallelIntervals.update(
                job_id=unique_job_id,
                started_date=started_date)
            .where(
                job_id==0)
            .limit(1)

            if not interval_found:
                print("Finished all intervals. Exiting..")
                break

            interval_to_work_on_next = ParallelIntervals.get(
                job_id == unique_job_id,
                started_date == started_date)

        return_code = run("python " + script_name + i.chrom + i.start_pos + i.end_pos)
        if return_code == 0:
            interval_to_work_on_next.finished_successfully = 1
        else:
            interval_to_work_on_next.finished_with_error = return_code


else:
    # if this is being run for the first time, create all intervals
    db_utils.create_table(db, ParallelIntervals, indexes, safe=True)

    for interval in intervals:
        row, created = ParallelIntervals.get_or_create(
            chrom=chrom, start_pos=start_pos, end_pos=end_pos)

        if created:
            # only set current_pos if the table didn't exist previously
            #row.current_pos = start_pos
            row.save()
    # start array job with N jobs, running self




p = argparse.ArugmentParser()
p.add_argument("-L", "--interval-list", help="An interval file")
p.add_argument("-n", "--num-jobs", help="Number of array job tasks to launch")
p.add_argument("command", help="The command to parallelize. The command must work with --chrom, --start, --end appended")

# generate database of intervals


#p.add_argument("--chrom", help="Genomic region - chromosome")
#p.add_argument("--start", "--start", help="Genomic region - start (1-based)")
#p.add_argument("--end", "--end", help="Genomic region - end (1-based)")




script_name = "some_script.py"
interval_size_per_job =2*1000*1000 # bp

def run(cmd):
    print(cmd)
    os.system(cmd)

chrom_sizes = {
"1":249250621,
"2":243199373,
"3":198022430,
"4":191154276,
"5":180915260,
"6":171115067,
"7":159138663,
"8":146364022,
"9":141213431,
"10":135534747,
"11":135006516,
"12":133851895,
"13":115169878,
"14":107349540,
"15":102531392,
"16":90354753,
"17":81195210,
"18":78077248,
"19":59128983,
"20":63025520,
"21":48129895,
"22":51304566,
"X":155270560,
#"Y":59373566,
#"MT":16569,
}

total_bp = sum(chrom_sizes.values())
total_jobs = total_bp/interval_size_per_job
jobs_per_chrom = { chrom : int((total_jobs*chrom_size)/float(total_bp))+1 for chrom, chrom_size in chrom_sizes.items() }

# run parse_calling_intervals - no need to parallelize
run("qsub -q short -cwd -o /broad/hptmp/exac_readviz_backend/step1_step1_parse_calling_intervals.log -j y -V ./run_python.sh -u step1_parse_calling_intervals.py" % locals())


# run choose_samples in parallel by chromosome
counter = 0
for chrom in chrom_sizes:
    N = jobs_per_chrom[chrom]
    run("qsub -q short -cwd -t 1-%(N)s -j y -V ./run_python.sh -u step1_choose_samples.py -c %(chrom)s -n %(N)s" % locals())
    counter += N

print("%s jobs total" % counter)


            runlog.started = datetime.datetime.now()
            runlog.save()

            start_pos = runlog.current_pos   # pick up where the previous job left off

            runlog.current_pos = minrep_pos
            runlog.heartbeat = datetime.datetime.now()
            runlog.save()


    if args.chrom:
        runlog.succeeded = 1
        runlog.save()

    db.close()


if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--info-table",        help="Path of ExAC info table",
        default='/humgen/atgu1/fs03/lek/resources/ExAC/ExAC.r0.3_meta_Final.tsv')
    p.add_argument("--full-vcf",
        help="Path of the ExAC full vcf with genotypes",
        default='/humgen/atgu1/fs03/konradk/exac/gqt/exac_all.vcf.gz')
    p.add_argument("-c", "--chrom",
        help="If specified, only this chromosome will be processed. Can start with 'chr' or not.")
    p.add_argument("-i", "--job-i", type=int,
        help="If specified, the chromomse will be divided into n equally-sized " \
                   "intervals and only variants in interval i will be processed. " \
                   "Ignored unless --chrom is also specified.", default=int(os.getenv("SGE_TASK_ID", "1").replace("undefined", "1")))
    p.add_argument("-n", "--total-jobs", type=int,
        help="If specified, the chromomse will be divided into n equally-sized " \
                   "intervals and only variants in interval i will be processed. " \
                   "Ignored unless --chrom is also specified.")
    args = p.parse_args()

    # print out settings
    print("Running with settings: ")
    for argname in filter(lambda n: not n.startswith("_"), dir(args)):
        print("   %s = %s" % (argname, getattr(args, argname)))
    print("\n")

    if args.chrom:
        chrom_sizes = {
            "1":249250621,
            "2":243199373,
            "3":198022430,
            "4":191154276,
            "5":180915260,
            "6":171115067,
            "7":159138663,
            "8":146364022,
            "9":141213431,
            "10":135534747,
            "11":135006516,
            "12":133851895,
            "13":115169878,
            "14":107349540,
            "15":102531392,
            "16":90354753,
            "17":81195210,
            "18":78077248,
            "19":59128983,
            "20":63025520,
            "21":48129895,
            "22":51304566,
            "X":155270560,
            "Y":59373566,
            "MT":16569,
        }


        args.chrom = args.chrom.replace("chr", "").upper()
        if args.chrom not in chrom_sizes:
            p.error("Invalid chromosome name: " + args.chrom)

        if args.job_i is not None and args.total_jobs is not None:
            interval_size = int(chrom_sizes[args.chrom] / args.total_jobs) # round up
            args.start_pos = (args.job_i - 1)*interval_size
            args.end_pos = min(1 + args.job_i*interval_size, chrom_sizes[args.chrom] + 1)
        else:
            args.start_pos = 0
            args.end_pos = chrom_sizes[args.chrom] + 1

    main(args)
