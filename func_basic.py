import numpy as np

class JacobiSolver:

    def __init__(self, A, b, x0, iterMax=100, errMax=1E-8):
        
        self.A = A
        self.b = b
        self.x = x0
        self.iterMax = iterMax
        self.errMax = errMax
    
    def solve(self):
        
        self.D = np.diagonal(self.A, offset=0)
        self.LU = self.A.copy()
        for i in range(len(self.D)):
            self.LU[i][i] -=  self.D[i]
        
        while 1:
            temp = self.b - np.matmul(self.LU, self.x).astype(np.float32)
            x = self.x.copy().astype(np.float32)
            for i in range(len(self.D)):
                print(x)
                x[i] = 1 / self.D[i] * temp[i]
                print(1 / self.D[i] * temp[i])
                print(x)
            # x = np.array([1 / self.D[i] * temp[i] for i in range(len(self.D))])
            print(x)
            if np.max(np.abs(x-self.x)) < self.errMax:
                self.x = x
                break
            else:
                self.x = x

if __name__ == '__main__':
    
    A = np.array([[2, 1, 1], 
                  [0, 2, 1], 
                  [0, 1, 1]])
    b = np.array([2, 1, 1])
    x0 = np.array([0, 1, 0])
    
    J = JacobiSolver(A, b, x0)
    J.solve()
    print(J.x)
    