import peewee
import playhouse.pool

readviz_db = playhouse.pool.PooledMySQLDatabase('exac_readviz', user='root', host='dmz-exac-dev.broadinstitute.org', port=3307)
#_calling_intervals_db = peewee.MySQLDatabase('exac_readviz2', user='root', host='dmz-exac-dev.broadinstitute.org', port=3307, autocommit=False)
#_calling_intervals_db = peewee.SqliteDatabase('exac_readviz.db', autocommit=False)
#_calling_intervals_db = peewee.SqliteDatabase(':memory:')

# define database model for storing the calling intervals
class ExacCallingInterval(peewee.Model):
    chrom = peewee.CharField(max_length=5, null=False, index=True)  # without 'chr'
    start = peewee.IntegerField()
    end = peewee.IntegerField()
    strand = peewee.CharField(max_length=1)
    name = peewee.CharField(max_length=500)

    def __str__(self):
        return "%s:%s-%s" % (self.chrom, self.start, self.end)

    class Meta:
        database = readviz_db
        #primary_key = peewee.CompositeKey('chrom', 'start', 'end', 'strand', 'name')

indexes = (
    (('chrom', 'start', 'end', 'strand'), True),  # True means unique index
)

