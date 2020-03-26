# Time from reference epoch
def dt(t, tref):

    dt = t-tref

    # Correction for beginning or end of week crossovers
    if dt > 302400:
        dt = dt - 604800
    elif dt < -302400:
        dt = dt + 604800

    return dt
