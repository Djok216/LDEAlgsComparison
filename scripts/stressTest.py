import time
import os
import subprocess
import sys
from random import randint, seed


exes = ["./ldegraphmain", "./slopesV7i"]


def checkSkip(n, m, mxVal):
  if len(sys.argv) > 1:
    mx1, mx2, n1, m1, n2, m2 = map(int, sys.argv[1:])
    if mxVal < mx1 or mxVal > mx2:
      return True
    if n < n1 or n > n2:
      return True
    if m < m1 or m > m2:
      return True
  if mxVal == 1021:
    if n == 1 and m > 5: return True
    if n == 2 and m > 4: return True
    if n > 2: return True
  if mxVal == 503:
    if n == 1 and m > 7: return True
    if n == 2 and m > 5: return True
    if n == 3 and m > 3: return True
    if n == 4: return True
  if mxVal == 107:
    if n == 1 and m == 9: return True
    if n == 2 and m > 6: return True
    if n == 3 and m > 5: return True
    if n == 4 and m == 5: return True
  return False


def genTest(n, m, mxVal):
  return ' '.join(map(str, sorted(randint(1, mxVal) for _ in range(n))[::-1])) + ' 0 ' + ' '.join(map(str, sorted([mxVal] + [randint(1, mxVal) for _ in range(m - 1)]))) + ' 0'


nms = [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9)]
nms += [(2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8)]
nms += [(3, 3), (3, 4), (3, 5), (3, 6), (4, 4), (4, 5)]
mxVals = [2, 3, 5, 13, 29, 39, 107, 503, 1021]

seed(216)

header = """
\\begin{figure}
\\begin{center}
\\begin{tabular}{c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c@{\hspace{0.5cm}}c}
\hline
& A & 2 & 3 & 5 & 13 & 29 & 39 & 107 & 503 & 1021 \\\\
N & M \\\\
\hline\n"""

finNoEps = open("tableNoEps.tex", "w")
finWithEps = open("tableWithEps.tex", "w")
finTimeouts = open("tableTimeouts.tex", "w")
finSumTime = open("tableSumTime.tex", "w")

for fin in [finNoEps, finWithEps, finTimeouts, finSumTime]:
  fin.write(header)

totalClasses, totalTime = 0, 0
graphWinNoEps, slopesWinNoEps = 0, 0
graphWinWithEps, slopesWinWithEps = 0, 0
graphTotalTimeouts, slopesTotalTimeouts = 0, 0
graphTotalSumTime, slopesTotalSumTime = 0, 0
for n, m in nms:
  if n == m:
    for fin in [finNoEps, finWithEps, finTimeouts, finSumTime]:
      fin.write("\\hline\n")
  for fin in [finNoEps, finWithEps, finTimeouts, finSumTime]:
    fin.write("%d & %d" % (n, m))
  for mxVal in mxVals:
    if checkSkip(n, m, mxVal):
      for fin in [finNoEps, finWithEps, finTimeouts, finSumTime]:
        fin.write(" & ")
      continue
    graphNoEps, slopesNoEps = 0, 0
    graphWithEps, slopesWithEps = 0, 0
    graphTimeouts, slopesTimeouts = 0, 0
    graphSumTime, slopesSumTime = 0, 0
    for __ in range(10):
      test = genTest(n, m, mxVal)
      rt, sols = [], []
      for exe in exes:
          times = []
          for __ in range(5):
              cmd = "echo %s | %s 1>temp.txt" % (test, exe)
              start = time.perf_counter()
              try:
                output = subprocess.check_output(cmd, timeout=600, shell=True)
              except:
                open("temp.txt", "w").write("-1")
              end = time.perf_counter()
              times.append(end - start)
              sols.append(int(open('temp.txt', 'r').readline()))
              sys.stderr.write("%s %s %d %.4lf\n" % (test, exe, sols[-1], times[-1]))
              sys.stderr.flush()
              if times[-1] > 15 or sols[-1] == -1:
                if 'graph' in exe and sols[-1] == -1:
                  graphTimeouts += 1
                if 'slopes' in exe and sols[-1] == -1:
                  slopesTimeouts += 1
                break
          if times[-1] > 15:
            rt.append(sum(times) / len(times))
          else:
            times.sort()
            rt.append(sum(times[1:-1]) / 3)
      if rt[0] < rt[1]:
        graphNoEps += 1
      elif rt[1] < rt[0]:
        slopesNoEps += 1
      if abs(rt[0] - rt[1]) < 1e-2:
        graphWithEps += 0.5
        slopesWithEps += 0.5
      elif rt[0] < rt[1]:
        graphWithEps += 1
      else:
        slopesWithEps += 1
      graphSumTime += rt[0]
      slopesSumTime += rt[1]
    totalClasses += 1
    totalTime += graphSumTime + slopesSumTime
    if graphNoEps > 7:
      graphWinNoEps += 1
    elif slopesNoEps > 7:
      slopesWinNoEps += 1
    if int(graphWithEps) > 7:
      graphWinWithEps += 1
    elif int(slopesWithEps) > 7:
      slopesWinWithEps += 1
    graphTotalTimeouts += graphTimeouts
    slopesTotalTimeouts += slopesTimeouts
    graphTotalSumTime += graphSumTime
    slopesTotalSumTime += slopesSumTime
    finNoEps.write(" & %d:%d " % (graphNoEps, slopesNoEps))
    finWithEps.write(" & %d:%d " % (int(graphWithEps), int(slopesWithEps)))
    finTimeouts.write(" & %d:%d " % (graphTimeouts, slopesTimeouts))
    finSumTime.write(" & %.1lf:%.1lf " % (graphSumTime, slopesSumTime))
  for fin in [finNoEps, finWithEps, finTimeouts, finSumTime]:
    fin.write('\\\\\n')

header = """
\hline
\end{tabular}
\caption{\label{%s}%s}
\end{center}
\end{figure}\n"""
finNoEps.write(header % ("fig1", "Comparison of Graph and Slopes algorithms. %d %d %d" % (graphWinNoEps, slopesWinNoEps, totalClasses)))
finWithEps.write(header % ("fig2", "Comparison of Graph and Slopes algorithms using an Epsilong of 0.01. %d %d %d" % (graphWinWithEps, slopesWinWithEps, totalClasses)))
finTimeouts.write(header % ("fig3", "Number of timeouts (10 minutes) for Graph and Slopes algorithms. %d %d" % (graphTotalTimeouts, slopesTotalTimeouts)))
finSumTime.write(header % ("fig4", "The total time spent for each algorithm (Graph and Slopes) on each test class. %.02lf %.02lf" % (graphTotalSumTime, slopesTotalSumTime)))
