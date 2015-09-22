# work around for unexpected / incorrect behavior by peewee ORM module
# peewee create_table doesn't seem to create compound indexes as expected, 
# and also doesn't allow engine to be set when creating a table (eg. engine=TokuDB)

import logging
import peewee


def create_table(db, model, indexes=tuple(), safe=True):
    """Utility method for creating a database table and indexes.

    Args:
      db: an active peewee database connection
      model: subclass of peewee.Model
      indexes: a tuple of indexes that would normally be specified in 
         the peewee.Mode's class Meta. Example:

         indexes = (
            (('chrom', 'start', 'end'), True),  # True means unique index
          )
      safe: if True, no error will be raised if the table already exists
    """
    # create table as a compressed TokuDB table
    raw_query = db.compiler().create_table(model, safe=safe)
    raw_query = list(raw_query)
    #raw_query[0] = raw_query[0] + " engine=TokuDB, compression='tokudb_zlib', charset=latin1"
    db.execute_sql(*raw_query)
    logging.debug(raw_query[0])

    # create any indexes
    safe_str = "IF NOT EXISTS" if safe else ""
    model_name = model.__name__.lower()
    for i, (columns, unique) in enumerate(indexes):
        columns_str = ",".join(map(lambda c: "`"+c+"`", columns))
        unique_str = "UNIQUE" if unique else ""
        q = "ALTER TABLE `%(model_name)s` ADD %(unique_str)s KEY %(safe_str)s `index%(i)s`(%(columns_str)s);" % locals()
        logging.debug(q)
        db.execute_sql(q)
