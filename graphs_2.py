import numpy as np
import matplotlib.pyplot as plt
from process_dump_file import Processor

TOTAL = 8000000
INTERVAL = 5000
NO_MOLECULES = 500
IGNORE = 7000000
NO_REPEATS = 20

def process(linking, seeds, energy):
  p = Processor()
  output = {}
  for lk in linking:
    rs = []
    ts = []
    ws = []
    timesteps = []
    for seed in seeds:
      print(seed)
      string = str(lk) + "_" + str(seed) + "_" + str(energy)
      p.set("dumps_2/dump_" + string + ".DNA", NO_MOLECULES, "outputs_2/output_" + string + ".dat")
      timestep, r, t, w = p.process()
      rs.append(r)
      ts.append(t)
      ws.append(w)
      timesteps = timestep
    output[lk] = [timesteps, rs, ts, ws]
  np.save("dict_" + str(energy) + ".npy", output)

def average_dynamics(energy):
  output = np.load("dict_" + str(energy) + ".npy", allow_pickle=True).item()
  for lk in output:
    out = open('average_dynamics_' + str(lk) + '_' + str(energy), 'w')
    times = output[lk][0]
    rs = np.mean(output[lk][1], axis = 0)
    ts = np.mean(output[lk][2], axis = 0)
    ws = np.mean(output[lk][3], axis = 0)
    for i in range(len(times)):
      out.write(str(times[i]) + " " + str(rs[i]) + " " + str(ts[i]) + " " + str(ws[i]) + "\n")
    out.close()
    if (lk == 8):
      print(np.var(ts[350:]))
      
  


#process([8, 11], ['01234', 12345, 23456, 34567, 45678, 56789, 67890, 78901, 89012, 90123, 101234, 112345, 123456, 134567, 145678, 156789, 167890, 178901, 189012, 190123], 50)
#process([11], ['01234', 12345, 23456, 123456], 30)
average_dynamics(30)
average_dynamics(40)
average_dynamics(50)