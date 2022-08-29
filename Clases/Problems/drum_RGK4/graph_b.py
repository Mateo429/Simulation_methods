import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


data=pd.read_table("data_b.txt", header=None)
plt.plot(data[0], data[1], label = "$f(\lambda)$")
plt.title("$R(\lambda, r=1)=f(\lambda)$")
plt.xlabel("$\lambda$")
plt.ylabel("$f(\lambda)$")
plt.xlim(0, 15)
plt.legend()
plt.grid()
plt.savefig("Bessel_lambda_b.pdf")
