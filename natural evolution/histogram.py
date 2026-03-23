import matplotlib.pyplot as plt
import numpy as np

dt = 0.1
dx = 0.5
nx = 128

length_threshold = 2*dx*nx

def name(seed, field_type):
    return f'Data\\valsPerLoop_2d_{field_type}_nx_0_128_nt4800_seed{seed}.txt'

    
def new_func(seed, field_type):
    # Read data, skip the header line
    data = np.genfromtxt(name(seed, field_type), skip_header=1)
    
    # The last column contains the values of interest
    values = data[:, -1]
    
    # We need the time (j) corresponding to each row
    # j = i * dt, starting from 0
    # Use enumerate to get both index and value
    times = [i * dt for i in range(len(values))]
    
    # Find the first time when the value is less than 128
    # Use next() with a generator; if not found, default to None
    initial_time = next((t for t, v in zip(times, values) if v < length_threshold), None)
    
    # Find the first time when the value is equal to 0
    final_time = next((t for t, v in zip(times, values) if v == 0), None)
    
    # If both events were found, return the difference (initial_time - final_time)
    # Otherwise return -1 (as in the original logic when one list was empty)
    if initial_time is not None and final_time is not None:
        lif = final_time - initial_time
    else:
        lif = -1
    
    # plt.plot(data[125:,-1],label=str(seed), linewidth=0.1)
    
    return lif

# plt.show()



# Collect results for seeds 1 to 15
hist_2HDM = []
hist_4fields = []
hist_2fields = []
for i in range(1, 301):
    a = new_func(i,"2hdm")
    b = new_func(i,"4fields")
    c = new_func(i,"2fields")
    if a>0:
        hist_2HDM.append(a)
    if b>0:
        hist_4fields.append(b)
    if c>0:
        hist_2fields.append(c)

# Plot histogram of the computed differences
plt.hist(hist_2HDM, log=False, range=[0, 400], density=True,  bins=50, color='skyblue', alpha=0.8,label="2HDM")
plt.hist(hist_4fields, log=False, range=[0, 400], density=True,  bins=50, color='orange', alpha=0.8,label="4fields")
plt.hist(hist_2fields, log=False, range=[0, 400], density=True,  bins=50, color='green', alpha=0.8,label="2fields")
plt.ylabel('counts')
plt.xlabel('ball-period')
plt.legend()
plt.show()
