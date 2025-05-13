import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


df = pd.read_csv('/Users/xiangzhewei/Desktop/ASPHY/Course_Computational_Astrophysics/Target/gamer/bin/new_problem_HSES/Record__Conservation', comment='#', delim_whitespace=True, header=None)
# index 0 => col 1
# index 1 => col 2
# index 2 => col 3, and so on ...
time = df[0]
step = df[1]
mass_gas = df[2]
new_df = pd.DataFrame({
    'Time': time,
    'Step': step,
    'Mass_Gas': mass_gas
})
# output file in csv format, I don't want index
new_df.to_csv('mass_conservation.csv', index=False)
print(new_df.head())

print(mass_gas)

# --- draw figure ---
plt.figure(figsize=(8, 5))
plt.plot(time, mass_gas, marker='o', linestyle='-')
plt.xlabel('Time')
plt.ylabel('Mass_Gas')
plt.ylim(0.9297, 0.9299)
plt.title('Mass_Gas vs Time unstable')
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.6f'))
plt.grid(True)
plt.tight_layout()
plt.savefig('mass_gas_vs_time_unstable.png') 
plt.show()  

