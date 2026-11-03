
import sys, json, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import activation_validation as av
from multiprocessing import Pool
import numpy as np
av.SCHEMES = ["current", "FN", "MBN"]
if __name__ == "__main__":
    cs = av.cases()
    with Pool(int(sys.argv[2])) as p:
        res = p.map(av.case_eval, cs, chunksize=4)
    json.dump(res, open(sys.argv[1], "w"))
    print(len(res))
