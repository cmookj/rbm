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

R3 = mr.MatrixExp3(mr.VecToso3(np.array([1, 3, 5])))
print('R3')
print(R3)
R34 = mr.MatrixExp3(mr.VecToso3(np.array([-1, -2, -3])))
print('R34')
print(R34)
R4 = np.matmul(R3, R34)
print(R4)

R34_calc = np.matmul(mr.RotInv(R3), R4)
print(R34)
print(R34_calc)

