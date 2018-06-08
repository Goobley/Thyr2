import multiprocessing
import sys

try:
    with open('CoreCount', 'w') as f:
        f.write('%d\n' % multiprocessing.cpu_count())
        f.flush()
except IOError:
    sys.exit(1)

