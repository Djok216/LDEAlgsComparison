from os import system
import sys

system("gcc -static slopesV7i.c -std=c11 -O3 -o slopesV7i")
system("g++ -static -lm -s -x c++ -std=c++17 -O3 -o ldegraphmain ldegraphmain.cpp ../src/ldegraphalg.cpp ../src/ldealg.cpp")
system("python3 generate.py %s" % ' '.join(sys.argv[1:]))