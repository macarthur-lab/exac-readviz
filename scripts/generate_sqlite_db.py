"""
Copies data from the MySQL database variant table into a SQLite db that can be
used by the web backend to locate the bams to display
"""

import peewee

from utils.database import Variant

sqlite_db = peewee.SqliteDatabase('exac_readviz.db', autocommit=False)
class T(Variant):
    class Meta:
        database = sqlite_db
        indexes = (
            (('chrom', 'pos', 'ref', 'alt', 'het_or_hom'), True), # True means unique index
        )

#if not T.table_exists():
#    T.delete().execute()

T.create_table(fail_silently=True)

# copy the table
sqlite_db.connect()
with sqlite_db.atomic():
    for v in Variant.select().dicts():
        T.insert(**v).execute()
sqlite_db.close()

#CREATE TABLE t(
# chrom text,
# minrep_pos integer,
# minrep_ref text,
# minrep_alt text,
# n_het integer,
# n_hom_alt integer,
# reassembled_bams_het text,
# reassembled_bams_hom text,
# finished bool);

# CREATE UNIQUE INDEX variant_idx ON t(
#   chrom,
#   minrep_pos,
#   minrep_ref,
#   minrep_alt);
