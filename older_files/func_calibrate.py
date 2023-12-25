import numpy as np

def mse_fit(x, y, dx=None, dy=None, order=1):

    if dx is None:
        dx = np.zeros(x.shape)
    if dy is None:
        dy = np.zeros(y.shape)
    dX = np.diag( dx ) 
        
    x_x, y_order = np.meshgrid(x, np.arange(order+1))
    x_x, y_order = x_x.T, y_order.T
    X = x_x ** y_order
    JX = x_x ** np.maximum(y_order-1, 0) * y_order
    
    coefs = np.linalg.inv(X.T@X) @ X.T @ y
    
    db1 = np.linalg.inv(X.T@X) @ X.T @ dy
    db2 = np.linalg.inv(X.T@X) @ (dX@JX).T @ y 
    db3 = np.linalg.inv(X.T@X) @ ( (dX@JX).T@X + X.T@(dX@JX) ) @ np.linalg.inv(X.T@X) @ X.T @ y
    db = db1 + db2 + db3
    
    var1 = (np.linalg.inv(X.T@X) @ X.T)**2 @ dy
    var2 = (np.linalg.inv(X.T@X) @ JX.T)**2 @ dX.T @ y**2
    var3 =  (np.linalg.inv(X.T@X) @ (JX.T+X.T))**2 @ dX.T @ ((X+JX) @ np.linalg.inv(X.T@X) @ X.T @ y)**2
    var = (var1 + var2 + var3)**0.5
    
    return coefs, db, var

if __name__ == '__main__':
    
    np.random.seed(100)
    x = np.random.normal(1, 0.02, (20,)) * np.linspace(0, 10, 20)
    y = np.random.normal(1, 0.002, (20,)) * (np.linspace(0, 10, 20)*5 + 4)
    dx = abs(np.random.normal(0, 0.05, (20,)))
    dy = abs(np.random.normal(0, 0.05, (20,)))
    coefs, db, var = mse_fit(x, y, dx, dy)
    print(coefs, db, var)

