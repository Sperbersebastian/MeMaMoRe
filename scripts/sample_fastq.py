#!/usr/bin/env python3
import argparse, gzip, random, sys, os
def open_auto(p, m): 
    # Use larger buffer size for better I/O performance
    bufsize = 1024 * 1024  # 1MB buffer
    if p.endswith(".gz"):
        return gzip.open(p, m, compresslevel=6)
    return open(p, m, buffering=bufsize)

def it(f):
    while True:
        h=f.readline()
        if not h: return
        s=f.readline(); p=f.readline(); q=f.readline()
        if not q: return
        yield h,s,p,q

def wr(o, r): 
    # Write all four lines at once to reduce I/O calls
    o.write(r[0])
    o.write(r[1])
    o.write(r[2])
    o.write(r[3])

def single(i1,o1,prop,seed):
    random.seed(seed); t=k=0
    with open_auto(i1,"rb") as f1, open_auto(o1,"wb") as o:
        for r in it(f1):
            t+=1
            if random.random()<prop: k+=1; wr(o,r)
    print(f"[single] kept {k}/{t} ({(k/max(1,t)):.2%})", file=sys.stderr)

def paired(i1,i2,o1,o2,prop,seed):
    random.seed(seed); t=k=0
    with open_auto(i1,"rb") as f1, open_auto(i2,"rb") as f2, open_auto(o1,"wb") as o1h, open_auto(o2,"wb") as o2h:
        for r1,r2 in zip(it(f1), it(f2)):
            t+=1
            if random.random()<prop: k+=1; wr(o1h,r1); wr(o2h,r2)
    print(f"[paired] kept {k}/{t} ({(k/max(1,t)):.2%})", file=sys.stderr)
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--in1", required=True); ap.add_argument("--in2")
    ap.add_argument("--out1", required=True); ap.add_argument("--out2")
    ap.add_argument("-p","--proportion", type=float, default=0.10)
    ap.add_argument("-s","--seed", type=int, default=42)
    a=ap.parse_args()
    out1_dir = os.path.dirname(a.out1)
    if out1_dir:
        os.makedirs(out1_dir, exist_ok=True)
    if a.in2:
        if not a.out2: ap.error("When --in2 is provided, --out2 is required.")
        paired(a.in1,a.in2,a.out1,a.out2,a.proportion,a.seed)
    else:
        single(a.in1,a.out1,a.proportion,a.seed)
