import numpy as np
import matplotlib.pyplot as plt
from process_dump_file import Processor

TOTAL = 7500000
INTERVAL = 5000
NO_MOLECULES = 500
IGNORE = 3800000
NO_REPEATS = 5

def process(linking, seeds, energy):
  p = Processor()
  output = {}
  for lk in linking:
    rs = []
    ts = []
    ws = []
    timesteps = []
    for seed in seeds:
      string = str(lk) + "_" + str(seed) + "_" + str(energy)
      p.set("dumps/dump_" + string + ".DNA", NO_MOLECULES, "outputs/output_" + string + ".dat")
      timestep, r, t, w = p.process()
      rs.append(r)
      ts.append(t)
      ws.append(w)
      timesteps = timestep
    output[lk] = [timesteps, rs, ts, ws]
  np.save("dict_" + str(energy) + ".npy", output)

def get_average(array):
  return np.mean(array, axis = 0)

def flat(arrays):
  a_filter = []
  for a in arrays:
    a_filter.append(a[int(IGNORE/INTERVAL):])
  return np.array(a_filter).flatten()

def t_cor(repeats):
  cs = []
  for r in repeats:
    cs.append(correlation(TOTAL, r[int(IGNORE/INTERVAL):]))
  cs = get_average(cs)
  cs = cs[:4]
  slope, _ = np.polyfit(np.arange(0, INTERVAL * len(cs), INTERVAL), np.log(cs), 1)
  return -1 / slope

def calculate_averages(energy):
  data = np.load("dict_" + str(energy) + ".npy", allow_pickle=True).item()

  out = open("averages_" + str(energy) + ".dat", 'w')
  out.write("# linking number, radius of gyration, gyration error, twist, twist error, writhe error,\n")
  for lk in data:
    out.write(str(lk))
    for repeats in data[lk][1:]:
      out.write(", " + str(np.mean(flat(repeats))) + ", " + str(error(repeats)))
    out.write("\n")
  out.close()

def error(repeats):
  stdDev = np.mean(np.square(flat(repeats))) - np.mean(flat(repeats)) ** 2
  N = NO_REPEATS * TOTAL / t_cor(repeats)
  return np.sqrt(stdDev) / np.sqrt(N)

def correlation(T, x):
  cs = []
  sq_mean = np.mean(x) ** 2
  mean_sq = np.mean(np.square(x))
  for t in range(0, T - IGNORE + INTERVAL, INTERVAL):
    c = 0
    counter = 0
    for t0 in range(0, T - t - IGNORE + INTERVAL, INTERVAL):
      c += x[int(t0 / INTERVAL)] * x[int((t + t0) / INTERVAL)]
      counter += 1
    c /= counter
    c -= sq_mean
    c /= (mean_sq - sq_mean)
    cs.append(c)
  return cs


def gyration_fluctuation(energy):
  data = np.load("dict_" + str(energy) + ".npy", allow_pickle=True).item()
  rGs = get_average(data[2][1])[int(IGNORE/INTERVAL):]
  out = open('gyration.dat', 'w')
  for rG in rGs:
    out.write(str(rG) + "\n")
  out.close()
  plt.hist(rGs, bins = 50)
  plt.show()


process([2,3,4,5,6], [12345, 23456, 34567, 45678, 56789], 50)
process([2,3,4,5,6], [12345, 23456, 34567, 45678, 56789], 30)
process([2,3,4,5,6], [12345, 23456, 34567, 45678, 56789], 40)
calculate_averages(30)
calculate_averages(40)
calculate_averages(50)
#gyration_fluctuation(50)