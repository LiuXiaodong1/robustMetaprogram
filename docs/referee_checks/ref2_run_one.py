import sys, importlib.util, numpy as np
a = sys.argv[1:]
label, d, ke, seed, R = a[0], float(a[1]), float(a[2]), int(a[3]), int(a[4])
sys.argv = ['x', str(seed)]
spec = importlib.util.spec_from_file_location('m', 'ref2_invarmp_evalue.py'); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
rows = [m.one_rep(d, kappa_err=ke) for _ in range(R)]
m.summarise(label, rows)
