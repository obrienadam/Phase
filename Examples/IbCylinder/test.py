import numpy as np
from mat import A

print(np.sum(A[9,:]))
print(A.shape)

with open('coeffs.dat', 'w') as f:
  print(np.array(A[13,:]).flatten().shape)

  for i, coeff in enumerate(np.array(A[13,:]).flatten()):
    f.write('{}, {}\n'.format(i, coeff))
