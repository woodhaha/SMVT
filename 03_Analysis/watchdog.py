#!/usr/bin/env python3
"""Self-healing watchdog for the SMVT docking pipeline.

Stages (each resumable via checkpoint):
  1. main run:      dock.py --input combined_candidates.csv      (target 10569)
  2. anticancer:    dock.py --input anticancer_new_input.csv      (target = rows-1)
  3. merge:         merge_anticancer.py                          (output marker file)

Loops every 5 min; restarts any dead-but-not-done stage; exits when merge done.
"""
import subprocess, json, os, time

BASE = "/root/docking"
MAIN_CK = f"{BASE}/results.csv.checkpoint.json"
AC_CK = f"{BASE}/anticancer_results.csv.checkpoint.json"
AC_IN = f"{BASE}/anticancer_new_input.csv"
MERGE_OUT = f"{BASE}/anticancer_combined_results.csv"
WD_LOG = f"{BASE}/watchdog.log"
MAIN_TARGET = 10569


def log(msg):
    line = f"{time.strftime('%F %T')} {msg}"
    with open(WD_LOG, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def ck_done(path):
    try:
        c = json.load(open(path))
        return len(c.get("completed", [])) + len(c.get("failed", []))
    except Exception:
        return 0


def alive(pattern):
    r = subprocess.run(["pgrep", "-f", pattern], capture_output=True)
    return r.returncode == 0


def launch(cmd, logfile):
    with open(logfile, "a") as out:
        subprocess.Popen(cmd, cwd=BASE, stdout=out, stderr=out,
                         start_new_session=True, close_fds=True)


def ac_target():
    try:
        return sum(1 for _ in open(AC_IN)) - 1
    except Exception:
        return 2478


def main():
    log(f"watchdog started (main target {MAIN_TARGET})")
    AC_TARGET = ac_target()
    while True:
        main = ck_done(MAIN_CK)

        if main < MAIN_TARGET:
            # stage 1: main run
            if not alive("combined_candidates"):
                log(f"[recover] main run dead at {main}/{MAIN_TARGET} -> restarting")
                launch(["python3", "dock.py", "--input", "combined_candidates.csv",
                        "--output", "results.csv"], f"{BASE}/dock.log")
        else:
            ac = ck_done(AC_CK)
            if ac < AC_TARGET:
                # stage 2: anticancer run
                if not alive("anticancer_new_input"):
                    log(f"[recover] anticancer run dead at {ac}/{AC_TARGET} -> restarting")
                    launch(["python3", "dock.py", "--input", "anticancer_new_input.csv",
                            "--output", "anticancer_results.csv"], f"{BASE}/anticancer.log")
            else:
                # stage 3: merge
                if not os.path.exists(MERGE_OUT):
                    log("[recover] anticancer done, running merge")
                    with open(f"{BASE}/anticancer.log", "a") as out:
                        subprocess.run(["python3", "merge_anticancer.py"],
                                       cwd=BASE, stdout=out, stderr=out)
                else:
                    log("ALL DONE (merge output exists). watchdog exiting.")
                    break
        time.sleep(300)


if __name__ == "__main__":
    main()
