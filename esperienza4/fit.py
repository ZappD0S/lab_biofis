import numpy as np
import matplotlib.pyplot as plt
from scipy import fft, optimize, integrate, interpolate, stats, signal
from collections import defaultdict
from itertools import pairwise
from iminuit import Minuit
from labellines import labelLines


plt.rcParams["text.usetex"] = True


def compute_shifted_diff(y, b, offset):
    if offset > 0:
        y2 = y[offset:] - b[:-offset]
    elif offset < 0:
        y2 = y[:offset] - b[-offset:]
    else:
        y2 = y - b

    return y2


t0s = [-1e-7, -1e-6, -1e-4]
data = defaultdict(dict)


for pres in ["1", "01"]:
    base_path = f"1e01_20C/{pres}atmCO/"
    # base_path = f"2022dic19_Mb/{pres}atmCO/T20/"
    for i, t_scale in enumerate(["1us", "100us", "10ms"]):
        t_b, b = np.loadtxt(base_path + f"b{t_scale}.dat", skiprows=5, delimiter=",").T
        t_y, y = np.loadtxt(base_path + f"i{t_scale}.dat", skiprows=5, delimiter=",").T

        print(f"{t_scale}: {np.diff(t_y).mean()}")
        print(f"{t_scale}: {t_y.size}")

        if t_scale == "10ms":
            t_b2, b2 = np.loadtxt(
                base_path + f"b2_{t_scale}.dat", skiprows=5, delimiter=","
            ).T

            assert np.all(t_y == t_b2)

            # b2 /= b2[0]
            # y /= b2

            res = stats.linregress(t_b2, b2)
            b2_fit = res.intercept + res.slope * t_b2
            b2_fit /= b2_fit[0]

            y /= b2_fit

            # b2 -= b2.mean()
            # y -= b2

        loss = []
        offsets = np.arange(-30, 15)

        for offset in offsets:
            y2 = compute_shifted_diff(y, b, offset)
            y_filt = signal.savgol_filter(y2, 1501, 2)

            loss.append(np.std(y2 - y_filt))

        offset = offsets[np.argmin(loss)]
        y2 = compute_shifted_diff(y, b, offset)
        # y_filt = signal.savgol_filter(y2, 1501, 2)

        print(offset)

        # plt.plot(offsets, np.array(loss))
        # plt.show()

        if offset > 0:
            t_y2 = t_y[offset:]
        elif offset < 0:
            t_y2 = t_y[:offset]
        else:
            t_y2 = t_y

        t0 = t0s[i]
        ind = np.argmin(np.abs(t_y2 - t0))

        print("y2 mean:", y2[:ind].mean())
        print("y2 std:", y2[:ind].std())

        y2 = np.log(y2[:ind].mean() / y2)
        y2 /= y2.max()

        data[pres][t_scale] = [t_y2, y2]


eps = np.finfo(float).eps


def rescale(t_short, y_short, t_long, y_long):
    mask = (t_short.min() <= t_long) & (t_long <= t_short.max())
    # mask = (t_long.min() <= t_short) & (t_short <= t_long.max())

    t0 = -min(t_short[0], t_long[0]) + eps

    # logt = np.log(t_short + t0)

    f = interpolate.interp1d(np.log(t_short + t0), y_short)
    y_short_resamp = f(np.log(t_long[mask] + t0))

    # f = interpolate.interp1d(t_long, y_long)
    # y_long_resamp = f(t_short[mask])

    # y_short_resamp = np.interp(t_long[mask], t_short, y_short)

    res = optimize.least_squares(lambda k: y_short_resamp - k * y_long[mask], x0=(1,))
    # res = optimize.least_squares(lambda k: k * y_long_resamp - y_short[mask], x0=(1,))

    if not res.success:
        raise Exception()

    k = res.x.item()

    # plt.plot(np.log(t_long[mask] + t0), y_short_resamp)
    # plt.plot(np.log(t_long[mask] + t0), y_long[mask])
    # plt.plot(np.log(t_long[mask] + t0), k * y_long[mask])
    # plt.show()

    return k


for pres, pres_data in data.items():
    for t_scale_short, t_scale_long in pairwise(pres_data):
        [t_short, y_short] = pres_data[t_scale_short]
        [t_long, y_long] = pres_data[t_scale_long]

        k = rescale(t_short, y_short, t_long, y_long)
        print(k)

        mask = t_long > t_short.max()
        pres_data[t_scale_long] = [t_long[mask], k * y_long[mask]]


final_data = {}
# t_min = np.inf
# t_max = -np.inf

for pres, pres_data in data.items():
    ts, ys = [], []
    for t, y in pres_data.values():
        ts.append(t)
        ys.append(y)

    t, y = [np.concatenate(xs) for xs in (ts, ys)]

    # t_min = min(t[0], t_min)
    # t_max = max(t[-1], t_max)
    final_data[pres] = (t, y)


for pres, (t, y) in final_data.items():
    # t += -t_min + eps

    # logt = np.log(t)

    # y_offset = -y.min() + eps
    # logy = np.log(y + y_offset)

    idx = y.argmax()
    # y /= y.max()
    y /= 0.985

    t = t[idx:]
    y = y[idx:]

    final_data[pres] = (t, y)


def funMb1(t, y, *args):
    k_1, kout, kin, kc, k_c, MbCO, Mb, tra1, CO, MbCO2, Mb2, tra2, CO2 = args

    return [
        k_c * (tra1 + y[2])
        + (-k_1 - kc - kout) * (MbCO + y[0])
        + kin * (Mb + y[1]) * (CO + y[1]),
        kout * (MbCO + y[0]) - kin * (CO + y[1]) * (Mb + y[1]),
        kc * (MbCO + y[0]) - k_c * (tra1 + y[2]),
    ]


def funMb2(t, y, *args):
    k_1, kout, kin, kc, k_c, MbCO, Mb, tra1, CO, MbCO2, Mb2, tra2, CO2 = args

    return [
        k_c * (tra2 + y[2])
        + (-k_1 - kc - kout) * (MbCO2 + y[0])
        + kin * (Mb2 + y[1]) * (CO2 + y[1]),
        kout * (MbCO2 + y[0]) - kin * (CO2 + y[1]) * (Mb2 + y[1]),
        kc * (MbCO2 + y[0]) - k_c * (tra2 + y[2]),
    ]


args0 = [0.3e6, 0.65e7, 1.5e7, 1e7, 1e7, 17e-6, 0, 0, 9.7342e-4, 17e-6, 0, 0, 9.7342e-5]
# args0 = [299999.68750000006,
#  6499978.303914831,
#  14999971.389751548,
#  10000000.0,
#  10000004.768372059,
#  1.710507949615988e-05,
#  0.0,
#  0.0,
#  0.0009145161049738802,
#  1.7167775970350916e-05,
#  0.0,
#  0.0,
#  0.00017992550546792787]

args0 = [300003.0,
 6499863.624485266,
 14999999.99999659,
 10000004.768372059,
 10000000.0,
 1.710576771055941e-05,
 0.0,
 0.0,
 0.0009124997051676579,
 1.7167924891173963e-05,
 0.0,
 0.0,
 0.00017959201069608487]

def fit_fun(args):
    res1 = integrate.solve_ivp(
        funMb1,
        # t_span=[t1[0], t1[-1]],
        t_span=[t_min, t_max],
        y0=[0, 0, 0],
        method=method,
        t_eval=t1,
        vectorized=True,
        args=args,
        first_step=1e-10,
        max_step=1e-3,
    )

    if not res1.success:
        raise Exception()

    # y_fit1 = res1.y.sum(axis=0) + args[5]
    y_fit1 = res1.y.sum(axis=0) + sum(args[5:8])

    res2 = integrate.solve_ivp(
        funMb2,
        # t_span=[t2[0], t2[-1]],
        t_span=[t_min, t_max],
        y0=[0, 0, 0],
        method=method,
        t_eval=t2,
        vectorized=True,
        args=args,
        first_step=1e-10,
        max_step=1e-3,
    )

    if not res2.success:
        raise Exception()

    # y_fit2 = res2.y.sum(axis=0) + args[9]
    y_fit2 = res2.y.sum(axis=0) + sum(args[9:12])

    return np.sum((y1 - y_fit1) ** 2) + np.sum((y2 - y_fit2) ** 2)


(t1, y1), (t2, y2) = final_data.values()

t_min, t_max = min(t1[0], t2[0]), max(t1[-1], t2[-1])
y1 = y1 * 17e-6
y2 = y2 * 17e-6

method = "BDF"

# res = optimize.minimize(
#     fit_fun,
#     x0=args0,
#     bounds=[(0, None) for _ in range(len(args0))],
#     tol=1e-6,
#     # args=(t1, y1, t2, y2),
# )
# args_best = res.x

m = Minuit(fit_fun, args0)
m.limits = (0, None)
# m.tol = 1e-6
m.migrad()
args_best = list(m.values)


res1 = integrate.solve_ivp(
    funMb1,
    # t_span=[t1[0], t1[-1]],
    t_span=[t_min, t_max],
    y0=[0, 0, 0],
    method=method,
    t_eval=t1,
    vectorized=True,
    args=args_best,
    first_step=1e-10,
    max_step=1e-3,
)

# y_fit1 = res1.y.sum(axis=0) + args_best[5]
y_fit1 = res1.y.sum(axis=0) + sum(args_best[5:8])


res2 = integrate.solve_ivp(
    funMb2,
    # t_span=[t2[0], t2[-1]],
    t_span=[t_min, t_max],
    y0=[0, 0, 0],
    method=method,
    t_eval=t2,
    vectorized=True,
    args=args_best,
    first_step=1e-10,
    max_step=1e-3,
)

# y_fit2 = res2.y.sum(axis=0) + args_best[9]
y_fit2 = res2.y.sum(axis=0) + sum(args_best[9:12])

plt.plot(
    t1,
    y1,
    ".",
    label=r"dati sperimentali $p_{\mathrm{CO}} = 0.1\,\mathrm{atm}$",
    zorder=0,
)
plt.plot(t1, y_fit1, label=r"fit $p_{\mathrm{CO}} = 0.1\,\mathrm{atm}$", zorder=1)

plt.plot(
    t2,
    y2,
    ".",
    label=r"dati sperimentali $p_{\mathrm{CO}} = 1\,\mathrm{atm}$",
    zorder=0,
)
plt.plot(t2, y_fit2, label=r"fit $p_{\mathrm{CO}} = 1\,\mathrm{atm}$", zorder=1)


plt.xlabel("tempo (s)")
plt.ylabel("concentrazione [Mb:CO] + [Mb] + [Trap] (M)")
plt.xscale("log")

plt.legend()
plt.tight_layout()
plt.savefig("fit2.png", dpi=300, bbox_inches="tight")
plt.show()

plt.plot(t1, y_fit1, label="[Mb:CO] + [Mb] + [Trap]")
plt.plot(t1, res1.y[1] + args_best[6], label="[Mb]")
plt.plot(t1, res1.y[2] + args_best[7], label="[Trap]")


plt.xlabel("tempo (s)")
plt.ylabel("concentrazione (M)")
plt.xscale("log")

plt.legend()
plt.tight_layout()
plt.savefig("tot_mb_trap.png", dpi=300, bbox_inches="tight")
plt.show()
