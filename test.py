import modern_robotics as mr
import numpy as np 

w1 = np.array([2, -6, 4])
R1 = mr.MatrixExp3(mr.VecToso3(w1))
p1 = np.array([3, 5, 7])
T1 = mr.RpToTrans(R1, p1)
print(T1)

w2 = np.array([-9, 7, 5])
R2 = mr.MatrixExp3(mr.VecToso3(w2))
p2 = np.array([-3, 6, -9])
T2 = mr.RpToTrans(R2, p2)
print(T2)

T = np.matmul(T1, T2)
print(T)


