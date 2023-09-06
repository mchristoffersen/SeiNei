import seinei.ui
import seinei.fetch
import seinei.tilt
import seinei.output

__all__ = ["ui", "fetch", "tilt", "output"]

import logging

log = logging.getLogger("seinei")
log.setLevel("INFO")
log.propagate = 0
ch = logging.StreamHandler()
FORMAT = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
formatter = logging.Formatter(FORMAT, datefmt="%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
log.addHandler(ch)
