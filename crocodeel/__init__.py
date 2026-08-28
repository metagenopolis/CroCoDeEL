"""Configure the CroCoDeEL runtime environment."""

import os

# Prevent oversubscription when CroCoDeEL uses multiprocessing.
# NumPy/OpenMP is restricted to one thread per process.
os.environ.setdefault("OMP_NUM_THREADS", "1")
