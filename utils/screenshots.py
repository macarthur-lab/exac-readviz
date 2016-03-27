
import igv_api
import logging
import os
from utils.constants import EXAC_CALLING_INTERVALS_BED_PATH, GENCODE_BED_PATH, \
    SELF_CHAIN_BED_PATH, IGV_JAR_PATH, IGV_SCREEN_WIDTH, IGV_TRACK_HEIGHT


assert os.path.isfile(IGV_JAR_PATH), "Couldn't find IGV jar: %s" % IGV_JAR_PATH


def take_screenshots(chrom, pos, file_paths, file_output_dir=None):
    """
    Launches IGV and takes screenshots.

    Args:
      chrom: chromosome
      pos: pos
      file_paths: List of file paths to load.
    """

    logging.info("%s:%s - taking igv screenshots of %s files: %s" % (
        chrom, pos, len(file_paths), ", ".join(file_paths)))

    r = igv_api.IGVCommandLineRobot(verbose=False,
                                    igv_window_width=1200,
                                    igv_window_height=1600,
                                    igv_jar_path=IGV_JAR_PATH)

    r.new_session()
    r.max_panel_height(IGV_TRACK_HEIGHT*len(file_paths))
    r.load(file_paths)
    for path in [EXAC_CALLING_INTERVALS_BED_PATH, SELF_CHAIN_BED_PATH, GENCODE_BED_PATH]:
        if path: r.load([path])

    r.goto("%s:%s-%s" % (chrom, pos - IGV_SCREEN_WIDTH, pos + IGV_SCREEN_WIDTH))

    png_filename = ".".join(os.path.basename(file_paths[0]).split(".")[0:-1])+".png"
    r.screenshot(os.path.join(file_output_dir, png_filename))
    r.exit_igv()

    r.execute()

