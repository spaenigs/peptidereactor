def is_struc_based(e):
    if "asa" in e:
        return True
    elif "delaunay" in e:
        return True
    elif "delaun" in e:
        return True
    elif "disorder" in e:
        return True
    elif "disord" in e:
        return True
    elif "elect" in e:
        return True
    elif "hull" in e:
        return True
    elif "qsar" in e:
        return True
    elif "sse" in e:
        return True
    elif "ta" in e:
        return True
    else:
        return False


def path(paths, prop1, prop2): \
    return [p for p in paths if prop1 in p and prop2 in p][0]