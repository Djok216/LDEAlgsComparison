import subprocess
import sys


res = subprocess.Popen("python3 stressTest.py %s" % ' '.join(sys.argv[1:]), shell=True)
if res.wait() != 0:
  print("Error")