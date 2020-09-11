import numpy as np


def div(y1, y2, y):
    def wdiv(o1, o2, t):
        res = {0: 2, 1: 1, 2: 0}
        a = res[np.abs(o1-t)+np.abs(o2-t)]
        print((o1, o2, t, a))
        return a
    a = np.array([y1, y2, y]).transpose()
    return 1/(2*a.shape[0]) * sum([wdiv(*a[i, :]) for i in range(a.shape[0])])


a = np.array([[1, 1, 1], [0, 1, 1], [1, 0, 0], [0, 0, 0], [1, 0, 1]])
# a = np.array([[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [1, 1, 0]])
print(a)

print(div(a[:, 0], a[:, 1], a[:, 2]))
