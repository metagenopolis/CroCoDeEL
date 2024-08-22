import os
# Disable implicit parallelism in numpy
os.environ["OMP_NUM_THREADS"] = "1"
