import subprocess, sys, os, glob

# Find gaff-2.11.xml
paths = [
    "/usr/local/lib/python3.12/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml",
    "/usr/local/lib/python3.12/dist-packages/openmmforcefields/",
]

print("=== openmmforcefields location ===")
for p in glob.glob("/usr/local/lib/python*/dist-packages/openmmforcefields/ffxml/amber/gaff/ffxml/gaff-2.11.xml"):
    print(f"FOUND: {p}")

# Also check if openmmforcefields is installed
try:
    import openmmforcefields
    print(f"openmmforcefields: {openmmforcefields.__file__}")
except ImportError:
    print("openmmforcefields NOT installed")
    # Try to install it
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "openmmforcefields"])
    import openmmforcefields
    print(f"Installed openmmforcefields: {openmmforcefields.__file__}")

# Find GAFF XML
import openmmforcefields
base = os.path.dirname(openmmforcefields.__file__)
for root, dirs, files in os.walk(base):
    for f in files:
        if "gaff" in f.lower() and f.endswith(".xml"):
            print(f"GAFF XML: {os.path.join(root, f)}")
