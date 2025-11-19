#!/usr/bin/env python3
import argparse, sys, json, os
try:
    import yaml
except ImportError:
    sys.stderr.write("Install PyYAML in your tools env (e.g., env_sra_tools).\n"); sys.exit(2)

def deep_update(a,b):
    for k,v in (b or {}).items():
        if isinstance(v,dict) and isinstance(a.get(k),dict): deep_update(a[k],v)
        else: a[k]=v
    return a

def parse_set(items):
    out={}
    for kv in items or []:
        k,sep,v=kv.partition("=")
        if not sep: raise SystemExit(f"--set expects key=value, got: {kv}")
        try: v=json.loads(v)
        except: pass
        cur=out
        parts=k.split(".")
        for p in parts[:-1]: cur=cur.setdefault(p,{})
        cur[parts[-1]]=v
    return out

def load_yaml(p): 
    if not p or not os.path.isfile(p): return {}
    with open(p, 'r') as f: 
        return yaml.safe_load(f) or {}

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--module-defaults")
    ap.add_argument("--defaults")
    ap.add_argument("--local")
    ap.add_argument("--set", action="append")
    ap.add_argument("--get")
    ap.add_argument("--out")
    a=ap.parse_args()

    merged={}
    for src in (a.module_defaults,a.defaults,a.local):
        deep_update(merged, load_yaml(src))
    deep_update(merged, parse_set(a.set))

    if a.get:
        cur=merged
        for p in a.get.split("."): cur=cur[p]
        print(json.dumps(cur) if isinstance(cur,(dict,list)) else cur); sys.exit(0)

    s=yaml.safe_dump(merged, sort_keys=False)
    if a.out:
        with open(a.out, "w") as f:
            f.write(s)
    else:
        print(s)
