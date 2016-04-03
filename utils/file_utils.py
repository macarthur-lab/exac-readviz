import os
import time


_missing_files_cache = set()  # used to speed up the does_file_exist(..) method


def retry_if_IOError(f, *args, **kwargs):
    num_retries = 6
    for retry_counter in range(num_retries):
        try:
            return f(*args, **kwargs)
        except OSError as e:
            err = e
            time.sleep(1)
        except Exception as e:
            raise e  # some other error occurred, so just raise it
    raise err  # all retries failed


def does_file_exist(file_path, num_retries=3, use_cache=True):
    """Returns true if the given file exists, using num_retries to allow for
    possible transient NFS issues"""
    if use_cache and file_path in _missing_files_cache:
        return False

    for retry_counter in range(num_retries):
        if os.access(file_path, os.R_OK):
            return True
        time.sleep(1)
    else:
        if use_cache:
            _missing_files_cache.add(file_path)
        return False
