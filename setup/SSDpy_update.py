import shutil
import os

## this code copies the SSDpy package from the local folder to the venv folder so that it can be used in other projects

# define the source and destination directories
src_dir = "/Users/vincentdenoel/git/SSDpy"
dst_dir = "/opt/anaconda3/envs/SSDpy/lib/python3.9/site-packages/SSDpy"

print("\nMoving the SSD package files from\n"+src_dir+"\nto\n"+dst_dir+"\n")

# run SSDpy
with open(os.path.join(src_dir, "setup", "SSDpy.py")) as f:
    exec(f.read())

# remove contents of existing directory
if os.path.exists(dst_dir):
    shutil.rmtree(dst_dir)

# copy the entire directory tree from source to destination
shutil.copytree(src_dir, dst_dir)

print('Successfully done.')