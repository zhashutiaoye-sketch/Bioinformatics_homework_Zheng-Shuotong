import json, urllib.request, urllib.error, os, time, base64

API = "http://127.0.0.1:8010"
OUT = "C:/emp_hw6/out"
os.makedirs(OUT, exist_ok=True)

def post(path, body=None, timeout=1800):
    data = json.dumps(body or {}).encode()
    req = urllib.request.Request(API + path, data=data,
                                 headers={"Content-Type": "application/json"}, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            txt = r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return {"http_error": e.code, "body": e.read().decode("utf-8", "replace")[:2000]}
    try:
        return json.loads(txt)
    except Exception:
        return {"raw": txt[:2000]}

def get(path, timeout=1800):
    try:
        with urllib.request.urlopen(API + path, timeout=timeout) as r:
            txt = r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return {"http_error": e.code, "body": e.read().decode("utf-8", "replace")[:2000]}
    try:
        return json.loads(txt)
    except Exception:
        return {"raw": txt[:2000]}

def dump(name, obj, sub=""):
    d = os.path.join(OUT, sub) if sub else OUT
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, name)
    with open(p, "w", encoding="utf-8") as f:
        if isinstance(obj, str):
            f.write(obj)
        else:
            json.dump(obj, f, ensure_ascii=False, indent=1)
    return p

def save_png(b64, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(base64.b64decode(b64))
    return path

def wait_job(job_id, label="job", max_wait=3600, verbose=True):
    t0 = time.time()
    last = None
    while time.time() - t0 < max_wait:
        r = get("/api/jobs/%s" % job_id)
        j = r.get("job", r)
        st = j.get("status")
        if st != last:
            if verbose:
                print("   %s: %s %s%% %s" % (label, st, j.get("progress"), j.get("message", "")), flush=True)
            last = st
        if st in ("done", "failed", "error", "cancelled"):
            return j
        time.sleep(4)
    return {"status": "timeout", "job_id": job_id}

def run_job(path, body, label="job"):
    r = post(path, body)
    if "job_id" not in r:
        return r, None
    j = wait_job(r["job_id"], label)
    return r, j
