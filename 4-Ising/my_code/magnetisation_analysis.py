import numpy as np
import matplotlib.pyplot as plt

# power = np.linspace(0, 1, 21)
# for x in power :
#     plt.clf()
#     for i in range(10, 61, 10) :
#         data = np.loadtxt("data_L_{}.txt".format(i))
#         plt.plot(data[:,0], (i**x)*abs(data[:,2]), label=r"{} $\times$ {}".format(i, i))
#     plt.legend()
#     plt.grid(True)
#     plt.title(r"M $\times$ $L^{% .2f}$ vs $\beta$J" % x + "\n")
#     plt.xlabel(r"$\beta$ J")
#     plt.ylabel(r"M $\times$ L$^{% .2f}$" % x)
#     plt.savefig("M_times_L/magL_500_x_{p:.2f}.png".format(p=x))

power = np.linspace(0, 1, 21)
# for x in power :
#     plt.clf()
#     for i in range(20, 61, 10) :
#         data = np.loadtxt("data_L_modified_{}.txt".format(i))
#         plt.plot(data[:,0], (i**x)*abs(data[:,2]), label=r"{} $\times$ {}".format(i, i))
#     plt.legend()
#     plt.grid(True)
#     plt.title(r"M $\times$ $L^{% .2f}$ vs $\beta$J" % x + "\n")
#     plt.xlabel(r"$\beta$ J")
#     plt.ylabel(r"M $\times$ L$^{% .2f}$" % x)
#     plt.savefig("M_times_L/magL_2000_x_{p:.2f}.png".format(p=x))

for x in power :
    plt.clf()
    for i in range(20, 61, 10) :
        data = np.loadtxt("my_code/data_L_{}.txt".format(i))
        plt.plot(i*(np.reciprocal(data[:,0]) - 1/0.46), (i**x)*abs(data[:,2]), label=r"{} $\times$ {}".format(i, i))
    plt.legend()
    plt.grid(True)
    plt.title(r"M $\times$ $L^{% .2f}$ vs $L^{1/1}(T - T'_c)$" % x + "\n")
    plt.xlabel(r"$L^{1/1}(T - T'_c)$")
    plt.ylabel(r"M $\times$ L$^{% .2f}$" % x)
    plt.savefig("M_times_L/vmagL_500_x_{p:.2f}.png".format(p=x))

for x in power :
    plt.clf()
    for i in range(20, 61, 10) :
        data = np.loadtxt("my_code/data_L_modified_{}.txt".format(i))
        plt.plot(i*(np.reciprocal(data[:,0]) - 1/0.46), (i**x)*abs(data[:,2]), label=r"{} $\times$ {}".format(i, i))
    plt.legend()
    plt.grid(True)
    plt.title(r"M $\times$ $L^{% .2f}$ vs $L^{1/1}(T - T'_c)$" % x + "\n")
    plt.xlabel(r"$L^{1/1}(T - T'_c)$")
    plt.ylabel(r"M $\times$ L$^{% .2f}$" % x)
    plt.savefig("M_times_L/vmagL_2000_x_{p:.2f}.png".format(p=x))

# for v in range(1, 12) :
#     plt.clf()
#     for i in range(30, 61, 10) :
#         data = np.loadtxt("data_L_modified_{}.txt".format(i))
#         plt.plot((i**(1/v))*(data[:,0] - 0.46), abs(data[:,2]), label=r"{} $\times$ {}".format(i, i))
#     plt.legend()
#     plt.grid(True)
#     plt.title(r"M $\times$ $L^0$ vs $L^{1/\, %d}\beta$J" % v + "\n")
#     plt.xlabel(r"$L^{1/%d}\beta$J" % v)
#     plt.ylabel(r"M $\times$ $L^0$")
#     plt.savefig("critical_exponent/critical_exponent_2000_{p:.2f}.png".format(p=v))
