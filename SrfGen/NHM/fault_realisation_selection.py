import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

model_Mw = [5.5, 6.0, 8.0, 8.5]
model_NumSim = [10, 10, 50, 50]


input_file_name = "Fault_LeonardMag_Length_Width_Rake.txt"

cs_list_file_name = "CS_list_v18p6.txt"

output_file_name = "fault_selection_c18p6.txt"

output_image_file_name = "fig_NumSim_vs_Mw.png"


assert len(model_Mw) == len(model_NumSim), "model_MW and model_NumSim must be the same length"
assert len(model_Mw) > 1, "model_MW and model_NumSim must have length greater than 1"

# Name, magnitude
fault = []

for line in open(input_file_name):
    parts = line.split()

    fault.append((parts[0], float(parts[1])))


# Name
srf = []

for line in open(cs_list_file_name):
    srf.append(line.split()[0])

# Name, Magnitude, Num_sims
data = []

for event in fault:
    if event[0] in srf:
        data.append(list(event[0:2]))

for event in data:
    fault_Mw = event[1]
    num_sim = 0
    for i in range(0, len(model_Mw)-1):
        if (fault_Mw >= model_Mw[i]) and (fault_Mw < model_Mw[i+1]):
            num_sim = ((model_NumSim[i+1] - model_NumSim[i]) / (model_Mw[i+1] - model_Mw[i])) * \
                      (fault_Mw - model_Mw[i]) + model_NumSim[i]

    event.append(math.ceil(num_sim))

print("Total number of realisations =" + str(sum([event[2] for event in data])))


with open(output_file_name, 'w') as output_file:
    for event in data:
        output_file.write("{0} {1}r \n".format(event[0], event[2]))

data.sort(key=lambda x: (x[2], x[1]))

plt.figure()
plt.plot([x[1] for x in data], [y[2] for y in data], linewidth=1.5, markersize=8)
plt.xlabel("Rupture magnitude, Mw")
plt.ylabel("Number of realisations")
plt.axis([model_Mw[0], model_Mw[-2], 0,1.2*max([event[2] for event in data])])
plt.grid(True, "both")
#plt.axes().yaxis.set_major_locator(MultipleLocator(5))
plt.savefig(output_image_file_name, format="png")