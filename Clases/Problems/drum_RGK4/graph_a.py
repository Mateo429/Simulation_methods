import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_table("data.txt", header = None)
plt.plot(data[0], data[1], label = "R(r)")
plt.title("Bessel function with $\lambda = 1$")
plt.xlabel("r")
plt.ylabel("R(r)")
plt.xlim(0, 10)
plt.legend()
plt.grid()
plt.savefig("Bessel_lambda_1.pdf")
