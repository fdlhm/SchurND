# matlab_rss.py
# If you use it, please cite the corresponding paper:
# Carlota M. Cuesta, Francisco de la Hoz,
# A non-recursive Schur-Decomposition Algorithm for $N$-Dimensional Matrix Equations,
# arXiv:2412.15840, (2026).
import psutil, sys, time
MATLAB_process_ID = sys.argv[1]
MATLAB_process = psutil.Process(int(MATLAB_process_ID))
try:
    while MATLAB_process.is_running():
        print(f"{time.time():.2f} {MATLAB_process.memory_info().rss}", flush = True)
#        time.sleep(0.001)
except psutil.NoSuchProcess:
    pass
