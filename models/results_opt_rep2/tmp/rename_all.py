import os

prefix = "04_"

for f in os.listdir():
    if f.endswith(".py"):
        continue
    if (not f.startswith("01")) and (not f.startswith("02")) and (not f.startswith("03")):
        os.rename(f, prefix + f)
