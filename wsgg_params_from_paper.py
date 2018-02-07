import numpy as np

#----------------------------------------------------------

d = np.array([
[0.0340429,  0.0652305, -0.0463685,  0.0138684, -0.0014450],
[0.3509457,  0.7465138, -0.5293090,  0.1594423, -0.0166326],
[4.5707400,  2.1680670, -1.4989010,  0.4917165, -0.0542999],
[109.81690, -50.923590,  23.432360, -5.1638920,  0.4393889]
])

c = np.empty((4,5,5))

k = 0
c[:,:,k] = [
[0.7412956, -0.9412652, 0.8531866, -0.3342806, 0.0431436],
[0.1552073, 0.6755648, -1.1253940, 0.6040543, -0.1105453],
[0.2550242, -0.6065428, 0.8123855, -0.4532290, 0.0869309],
[-0.0345199, 0.4112046, -0.5055995, 0.2317509, -0.0375491] 
]

k = 1
c[:,:,k] = [
[-0.5244441, 0.2799577, 0.0823075, 0.1474987, -0.0688622],
[-0.4862117, 1.4092710, -0.5913199, -0.0553385, 0.0464663],
[0.3805403, 0.3494024, -1.1020090, 0.6784475, -0.1306996],
[0.2656726, -0.5728350, 0.4579559, -0.1656759, 0.0229520]
]

k = 2
c[:,:,k] = [
[0.5822860, -0.7672319, 0.5289430, -0.4160689, 0.1109773],
[0.3668088, -1.3834490, 0.9085441, -0.1733014, -0.0016129],
[-0.4249709, 0.1853509, 0.4046178, -0.3432603, 0.0741446],
[-0.1225365, 0.2924490, -0.2616436, 0.1052608, -0.0160047]
]

k = 3
c[:,:,k] = [
[-0.2096994, 0.3204027, -0.2468463, 0.1697627, -0.0420861],
[-0.1055508, 0.4575210, -0.3334201, 0.0791608, -0.0035398],
[0.1429446, -0.1013694, -0.0811822, 0.0883088, -0.0202929],
[0.0300151, -0.0798076, 0.0764841, -0.0321935, 0.0050463]
]

k = 4
c[:,:,k] = [
[0.0242031, -0.0391017, 0.0310940, -0.0204066, 0.0049188],
[0.0105857, -0.0501976, 0.0384236, -0.0098934, 0.0006121],
[-0.0157408, 0.0130244, 0.0062981, -0.0084152, 0.0020110],
[-0.0028205, 0.0079966, -0.0079084, 0.0033870, -0.0005364]
]

#----------------------------------------------------------

def get_k_a(T, P, XCO2, XH2O):

    nGG = 5            # 4 gray gases and 1 clear gas

    Mr = XH2O/XCO2
    Tr = T/1200

    K = np.zeros(nGG)
    a = np.zeros(nGG)

    b = np.zeros((4,5))
    for i in range(0,4):
        for j in range(0,5):
            b[i,j] = np.sum(c[i,j,:]*Mr**np.arange(5))

    for i in range(1, nGG):
        a[i] = np.sum(b[i-1,:]*Tr**np.arange(5))
    a[0] = 1-np.sum(a[1:])

    for i in range(1,nGG):
        K[i] = np.sum(d[i-1,:]*Mr**np.arange(5))

    return K*(P/101325)*(XH2O+XCO2), a


#----------------------------------------------------------

K, a = get_k_a(1000,101325,0.1,0.1)

print("K =", K)
print("a =", a)


