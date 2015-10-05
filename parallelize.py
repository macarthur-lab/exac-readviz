#!/usr/bin/env python3.4

"""
This script takes the name of another python script (let's say some_script.py)
and launches N parallel instances of some_script.py on the short queue as an
array job. Each instance of some_script.py will get passed
 --chrom, --start, and --end args which will define the genomic
region it should operate on. This region should be small enough for
some_script.py to finish in < 4 hours (the short queue's runtime limit) in the
worst case. It's up to some_script.py to avoid redoing the same work if
it is run multiple times on the same genomic interval (eg. if the 1st run fails).

some_script.py should return 0 if it completed successfully, or something
other than 0 if it fails.

Example:

python parallelize.py -L /seq/references/Homo_sapiens_assembly19/v1/variant_calling/exome_calling_regions.v1.interval_list -n 1000 python3.4 generate_HC_bams.py
"""


import argparse
import datetime
import os
import peewee
import playhouse.pool
import random
import slugify
import subprocess

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
#logging.getLogger('peewee').setLevel(logging.DEBUG)


p = argparse.ArgumentParser()
p.add_argument("-L", "--interval-list", help="An interval file")
p.add_argument("-n", "--num-jobs", help="Number of array job tasks to launch")
p.add_argument("-isize", "--interval-size", help="Max interval size", default=200)
p.add_argument("--log-dir", help="Logging directory", default="/broad/hptmp/exac_readviz_backend2/logs")
p.add_argument("-bsub", "--run-on-LSF", help="Submit to LSF", action="store_true")
p.add_argument("-local", "--run-local", help="Run locally instead of submitting array jobs", action="store_true")
p.add_argument("-f", "--regenerate-intervals-table", help="Regenerate intervals table from scratch", action="store_true")
p.add_argument("command", nargs="+", help="The command to parallelize. The command must work with --chrom, --start-pos, --end-pos")

args = p.parse_args()

db_table_name = "%s_i%d" % ("_".join(args.command[0:2]), args.interval_size)
db_table_name = slugify.slugify(db_table_name).replace("-", "_")  # remove special chars

args.command = " ".join(args.command)
#logging.info("args: command: " + args.command)
logging.info("db_table_name: " + db_table_name)

db = playhouse.pool.PooledMySQLDatabase('parallelize', user='root', host='dmz-exac-dev.broadinstitute.org', port=3307)

print("Starting..")

# table for keeping track of which intervals some_script.py has already processed
class ParallelIntervals(peewee.Model):
    chrom = peewee.CharField(max_length=5, null=False)
    start_pos = peewee.IntegerField(null=False)
    end_pos = peewee.IntegerField(null=False)
    #current_pos = peewee.IntegerField(default=0)  # not currently used because some_script.py has no way of feeding back this info

    job_id = peewee.IntegerField(null=True)  # cluster array job id
    task_id = peewee.IntegerField(null=True)  # cluster array job task id
    unique_id = peewee.IntegerField(null=True)  # unique id for this job task just in case the cluster re-uses the job_id and task_id

    started = peewee.BooleanField(default=0, index=True)
    started_date = peewee.DateTimeField(null=True, index=True)

    finished = peewee.BooleanField(default=0, index=True)
    finished_date = peewee.DateTimeField(null=True)

    error_code = peewee.IntegerField(default=0, index=True)
    error_message = peewee.TextField(null=True)

    # execution environment stats
    machine_hostname = peewee.CharField(null=True, max_length=100)
    machine_average_load = peewee.FloatField(null=True)

    class Meta:
        db_table = db_table_name
        database = db
        indexes=(
            (('chrom', 'start_pos', 'end_pos'), True), # True means unique index
            (('job_id', 'task_id', 'unique_id'), False),  # not unique because a given task can process multiple intervals
        )

array_job_task_id = os.getenv('SGE_TASK_ID', -1)

is_startup = array_job_task_id == -1 or array_job_task_id == "undefined"
if is_startup:
    # this instance of parallelize.py is being run for the first time
    if not args.interval_list:
        p.error("-L arg required")

    if not args.num_jobs and not args.run_local:
        p.error("--num-jobs arg required")

    # create intervals
    intervals = []
    with open(args.interval_list) as interval_list_file:
        for line in interval_list_file:
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            intervals.append({"chrom": chrom, "start_pos": start, "end_pos": end})

    logging.info("Parsed %s intervals from %s" % (len(intervals), args.interval_list))

    # split intervals so they are no bigger than args.interval_size
    final_intervals = []
    for interval in intervals:
        chrom, start, end = interval["chrom"], interval["start_pos"], interval["end_pos"]
        if end - start > args.interval_size:
            #print("Breaking up large interval %s:%s-%s" % (chrom, start, end))
            while end - start > args.interval_size:
                #print("   %s:%s-%s" % (chrom, start, start + args.interval_size-1))
                final_intervals.append({"chrom": chrom,
                                        "start_pos": start,
                                        "end_pos": start + args.interval_size-1})
                start = start + args.interval_size
            #print("   %s:%s-%s" % (chrom, start, end))

        final_intervals.append({"chrom": chrom, "start_pos": start, "end_pos": end})

    logging.info("Broke the %s intervals into %s intervals of size %s or less" % (len(intervals), len(final_intervals), args.interval_size))

    # populate database table if it doesn't exist already
    if ParallelIntervals.table_exists() and ParallelIntervals.select().count() == len(final_intervals) and not args.regenerate_intervals_table:
        logging.info("%s: %s intervals loaded previously" % (ParallelIntervals._meta.db_table, len(final_intervals)))
    else:
        if ParallelIntervals.table_exists():
            logging.info("Dropping existing table: %s" % ParallelIntervals._meta.db_table)
            ParallelIntervals.drop_table()
        logging.info("Creating table: " + str(ParallelIntervals._meta.db_table))
        ParallelIntervals.create_table()

        insert_batch_size = 25000  # a batch size that's too large will cause query to fail
        for batch_start in range(0, len(final_intervals), insert_batch_size):
            logging.info("Inserting records %s through %s" % (batch_start, batch_start + insert_batch_size))
            batch = final_intervals[batch_start: batch_start + insert_batch_size]
            ParallelIntervals.insert_many(batch).execute()


        assert ParallelIntervals.select().count() == len(final_intervals), \
            "Number of intervals in database (%s) != expected number (%s)" % (
                ParallelIntervals.select().count(), len(final_intervals))

    # start array job with N jobs, running self
    if not args.run_local:
        if not os.path.isdir(args.log_dir):
            os.system("mkdir -m 777 -p %s" % args.log_dir)

        if args.run_on_LSF:
            launch_array_job_cmd = (
                "bsub -N -J prog[1-%(num_jobs)s] -o %(log_dir)s -q hour "
                    "python3.4 parallelize.py %(command)s"
            )
        else:
            launch_array_job_cmd = ("qsub -q short "
                "-t 1-%(num_jobs)s "
                "-cwd "
                "-o %(log_dir)s "
                "-e %(log_dir)s "
                "-j y -V "
                "./run_python.sh python3.4 parallelize.py %(command)s")

        launch_array_job_cmd = launch_array_job_cmd  % {
            "num_jobs": args.num_jobs,
            "log_dir" : args.log_dir,
            "command" : args.command}
        logging.info("Running: %s" % launch_array_job_cmd)
        subprocess.check_call(launch_array_job_cmd, shell=True)

        # run several times
        subprocess.check_call(launch_array_job_cmd, shell=True)
        subprocess.check_call(launch_array_job_cmd, shell=True)
        subprocess.check_call(launch_array_job_cmd, shell=True)

        # TODO run loop that restarts array jobs

if not is_startup or args.run_local:
    # this instance of parallelize.py is running as one of many array job tasks.
    # run a loop that continually launches some_script.py on the next unprocessed
    # interval until times runs out for this task

    logging.info("USER: %s" % os.getenv('USER', ''))

    job_id = os.getenv('JOB_ID', os.getenv('LSB_JOBID', -1))
    if job_id != -1 and job_id != "undefined" and not args.run_local:
        logging.info("parallelize.py - job id: %s" % job_id)
    else:
        job_id = os.getpid()
        array_job_task_id = 0
        logging.info("parellelize.py - running as local process - id: %s" % job_id)

    #if not args.run_local:
    #    time.sleep(random.randint(1, 30)) # sleep between 0 and 60 seconds to avoid all tasks trying to aquire intervals at the same time

    unique_8_digit_id = random.randint(10**8, 10**9 - 1)  # don't use actual job id to avoid collisions in case this script has been restarted and the same job id is reused.
    while True:
        started_date = datetime.datetime.now()

        # get next interval
        current_interval = None
        with db.atomic() as txn:
            # claim an interval
            unprocessed_intervals = list(ParallelIntervals.select()\
                .where((
                    ParallelIntervals.job_id >> None) & (
                    ParallelIntervals.task_id >> None) & (
                    ParallelIntervals.unique_id >> None))\
                .limit(1) #.order_by(peewee.fn.Rand())\
            )

            if len(unprocessed_intervals) == 0:
                logging.info("Finished all intervals. Exiting..")
                break

            current_interval = unprocessed_intervals[0]
            current_interval.job_id = job_id
            current_interval.task_id = array_job_task_id
            current_interval.unique_id = unique_8_digit_id
            current_interval.started = 1
            current_interval.started_date = started_date
            current_interval.save()

        cmd = "%s --chrom %s --start-pos %s --end-pos %s" % (args.command,
            current_interval.chrom, current_interval.start_pos, current_interval.end_pos)
        logging.info("interval: %s:%s-%s - launching %s" % (
            current_interval.chrom, current_interval.start_pos, current_interval.end_pos, cmd))
        try:
            cmd_output = subprocess.check_output(cmd.split(" "), stderr=subprocess.STDOUT).decode()
            for line in cmd_output.split("\n"):
                logging.info("      %s" % line.strip())
            if "generate_HC_bams finished" not in cmd_output:
                raise subprocess.CalledProcessError(100, cmd, cmd_output)
        except subprocess.CalledProcessError as e:
            error_message = ("%s\n"
                             "return code: %s\n"
                             "output: %s") % (e.cmd, e.returncode, e.output.strip())
            current_interval.error_code = e.returncode
            current_interval.error_message = error_message
            current_interval.save()
            logging.info("interval: %s:%s-%s - failed: %s" % (current_interval.chrom, current_interval.start_pos, current_interval.end_pos, error_message))
        else:
            # execution environment stats
            current_interval.machine_hostname = os.getenv('HOSTNAME', '')[0:100]
            current_interval.machine_average_load = os.getloadavg()[-1]

            # finished
            current_interval.finished = 1
            current_interval.finished_date = datetime.datetime.now()
            current_interval.save()

            logging.info("interval: %s:%s-%s - succeeded!" % (current_interval.chrom, current_interval.start_pos, current_interval.end_pos))


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

#args.chrom = args.chrom.replace("chr", "").upper()
#if args.chrom not in chrom_sizes:
#    p.error("Invalid chromosome name: " + args.chrom)


# run parse_calling_intervals - no need to parallelize
#run("qsub -q short -cwd -o /broad/hptmp/exac_readviz_backend/step1_step1_parse_calling_intervals.log -j y -V ./run_python.sh -u step1_parse_calling_intervals.py" % locals())

