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

print('Matrix Projection to SO3')
R5 = mr.MatrixExp3(mr.VecToso3(np.array([3, 5, 7])))
print('R5: {}', np.linalg.det(R5))
print(R5)
M6 = R5 + np.array([[0.0001, -0.0001, -0.0001],
                    [-0.0001, 0.0001, 0.0001],
                    [0.0001, 0.0001, -0.0001]])
print('M6: {}', np.linalg.det(M6))
print(M6)
R6 = mr.ProjectToSO3(M6)
print('R6: {}', np.linalg.det(R6))
print(R6)

