import subprocess

scripts = ["momentum_plot.py", "ke_plot.py"] 

for script in scripts:
    for i in range(1, 13):
        arg = f"../data_run_{i}"
        subprocess.run(["python3", script, arg])
