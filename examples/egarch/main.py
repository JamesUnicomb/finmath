import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from egarch import loss, egarch_eval
from ewma import ewma_eval


df = pd.read_csv("data.csv", parse_dates=["date"], dayfirst=True).sort_values("date")
dates = df["date"].values
spx = df["spx"].values

deltas = np.diff(spx)

x0 = [0.0, 0.0, 0.0, 0.0, 0.0]
print(loss(0.1*np.random.randn(5), deltas))

res = minimize(loss, 0.1*np.random.randn(5), method='nelder-mead', args=(deltas),
               options={'disp': True})

# plt.plot(dates, spx)
# plt.show()

x = res.x

plt.plot([np.sqrt(v) for v in egarch_eval(x, deltas)])
plt.plot([1.0/np.sqrt(2.0)*np.sqrt(v) for v in ewma_eval(0.95, deltas)])
plt.show()