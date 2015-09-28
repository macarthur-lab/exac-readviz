import logging
import peewee
import playhouse.pool
from utils.constants import MAX_ALLELE_SIZE, MAX_VCF_SAMPLE_ID_SIZE

# disable peewee warning messages
logging.getLogger('peewee').setLevel(logging.ERROR)

# define database
_readviz_db = playhouse.pool.PooledMySQLDatabase(
    'exac_readviz',
    user='root',
    host='dmz-exac-dev.broadinstitute.org',
    port=3307)
#_readviz_db = peewee.MySQLDatabase('exac_readviz2', user='root', host='dmz-exac-dev.broadinstitute.org', port=3307, autocommit=False)
#_readviz_db = peewee.SqliteDatabase('exac_readviz.db', autocommit=False)
#_readviz_db = peewee.SqliteDatabase(':memory:')

logging.info("db: %s(%s), autocommit=%s" % (
    type(_readviz_db).__name__,
    _readviz_db.database,
    _readviz_db.get_autocommit()))


# common Meta parameters for all models
class _SharedMeta(peewee.Model):
    class Meta:
        database = _readviz_db


# create table for exac calling intervals
class ExacCallingInterval(_SharedMeta):
    chrom = peewee.CharField(max_length=5, null=False, index=True)  # without 'chr'
    start = peewee.IntegerField()
    end = peewee.IntegerField()
    strand = peewee.CharField(max_length=1)
    name = peewee.CharField(max_length=500)

    def __str__(self):
        return "%s:%s-%s" % (self.chrom, self.start, self.end)

    class Meta:
        indexes=(
            (('chrom', 'start', 'end', 'strand'), True),  # True means unique index
        )


# fields shared by Variant and Sample tables
class _SharedVariantFields(_SharedMeta):
    chrom = peewee.CharField(max_length=5, null=False, index=True)
    pos = peewee.IntegerField()
    ref= peewee.CharField(max_length=MAX_ALLELE_SIZE)
    alt = peewee.CharField(max_length=MAX_ALLELE_SIZE)
    het_or_hom = peewee.CharField(max_length=3)
    # finished is 1 if all processing steps for this variant have finished (either successfully or with non-transient errors)
    finished = peewee.BooleanField(default=0)


# create table for non-sensitive variant-level info for all ExAC variants -
# it will later be exported to a web-accessible sqlite db that can be queried by
# the readviz backend to generate the readviz frontend for each variant page
class Variant(_SharedVariantFields):
    # number of samples the user expects (this will be 5 or fewer)
    n_expected_samples = peewee.IntegerField(default=-1, index=True)
    # number of samples available (this may be < n_expected_samples if bams are missing or GVCF calls couldn't be reproduced for some samples)
    n_available_samples = peewee.IntegerField(default=-1, index=True)
    # list of bam paths to show (separated by '|') of size = n_available_samples
    readviz_bam_paths = peewee.TextField(default=None, null=True)


    class Meta:
        indexes = (
            (('chrom', 'pos', 'ref', 'alt', 'het_or_hom'), True), # True means unique index
        )


# create table for per-variant-sample info
# WARNING: this table contains sensitive info (eg. sample ids) and should not
# be publicly available
class Sample(_SharedVariantFields):
    sample_id = peewee.CharField(index=True, max_length=MAX_VCF_SAMPLE_ID_SIZE)
    sample_i = peewee.IntegerField(default=-1, null=True)

    original_bam_path = peewee.TextField(null=True)
    original_gvcf_path = peewee.TextField(null=True)
    output_bam_path = peewee.TextField(null=True)

    is_missing_original_gvcf = peewee.BooleanField(default=0, index=True)

    calling_interval_start = peewee.IntegerField(default=None, null=True)
    calling_interval_end = peewee.IntegerField(default=None, null=True)

    hc_succeeded = peewee.BooleanField(default=0, index=True)
    hc_error_code = peewee.IntegerField(default=None, index=True, null=True)
    hc_error_text = peewee.TextField(default=None, null=True)
    hc_n_artificial_haplotypes = peewee.IntegerField(default=None, index=True, null=True)

    hc_started_time = peewee.DateTimeField(default=None, null=True)
    hc_finished_time = peewee.DateTimeField(default=None, null=True)
    hc_command_line = peewee.TextField(default=None, null=True)

    #screenshot_started = peewee.BooleanField(default=0)
    #screenshot_finished = peewee.BooleanField(default=0)

    class Meta:
        indexes = (
            (('chrom', 'pos', 'ref', 'alt', 'het_or_hom', 'sample_id'), True), # True means unique index
        )


def _create_table(model, fail_silently=True):
    """Utility method for creating a database table and indexes that is a
    work around for unexpected behavior by the peewee ORM module. Specifically,
    peewee create_table doesn't create compound indexes as expected.

    Args:
      db: an active peewee database connection
      model: subclass of peewee.Model
      indexes: a tuple of indexes that would normally be specified in
         the peewee.Mode's class Meta. Example:

         indexes = (
            (('chrom', 'start', 'end'), True),  # True means unique index
          )
      fail_silently: if True, no error will be raised if the table already exists
    """
    # create table as a compressed TokuDB table
    db = model._meta.database
    indexes = model._meta.indexes

    raw_query = db.compiler().create_table(model, safe=fail_silently)
    raw_query = list(raw_query)

    raw_query[0] = raw_query[0] + " engine=TokuDB, compression='tokudb_zlib', charset=latin1"
    db.execute_sql(*raw_query)
    logging.debug(raw_query[0])

    # create indexes
    safe_str = "IF NOT EXISTS" if fail_silently else ""
    model_name = model.__name__.lower()
    for i, (columns, unique) in enumerate(indexes):
        columns_str = ",".join(map(lambda c: "`"+c+"`", columns))
        unique_str = "UNIQUE" if unique else ""
        q = "ALTER TABLE `%(model_name)s` ADD %(unique_str)s KEY %(safe_str)s `index%(i)s`(%(columns_str)s);" % locals()
        logging.debug(q)
        db.execute_sql(q)


def init_db():
    _create_table(ExacCallingInterval, fail_silently=True)
    _create_table(Variant, fail_silently=True)
    _create_table(Sample, fail_silently=True)

    #_readviz_db.connect()

    # print info about created tables
    #for table in (ExacCallingInterval._meta.db_table, Variant._meta.db_table, Sample._meta.db_table):
    #    logging.info("%s indexes: %s" % (table, _readviz_db.get_indexes(table),))

    #_readviz_db.close()

    return _readviz_db