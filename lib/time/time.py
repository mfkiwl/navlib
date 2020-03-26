# Time from reference epoch
def dt(t, t0):

    t = t - t0

    # Correction for beginning or end of week crossovers
    if t > 302400:
        t = t - 604800
    elif t < -302400:
        t = t + 604800

    return t
