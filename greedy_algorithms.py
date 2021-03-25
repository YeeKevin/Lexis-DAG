smalls = [1]
smallf = [0, 4]

start = [0, 1, 3, 0, 5, 3, 5, 6, 8, 8, 2, 12]
finish = [0, 4, 5, 6, 7, 9, 9, 10, 11, 12, 14, 16]

def recursive_activity_selector(s, f, k, n):
    """
        Args:
            s: a list of start times
            f: a list of finish times
            k: current position in
            n: total possible activities
        Returns:
            A maximal set of activities that can be scheduled.
            (We use a list to hold the set.)
    """
    m = k + 1
    while m < n and s[m] < f[k]:  # find an activity starting after our last
                                   # finish
        m = m + 1
    if m < n:
        print("Adding activity " + str(m) + " that finishes at "
              + str(f[m]))
        return [m] + recursive_activity_selector(s, f, m, n)
    else:
        return []


def greedy_activity_selector(s, f):
    """
        Args:
            s: a list of start times
            f: a list of finish times
        Returns:
            A maximal set of activities that can be scheduled.
            (We use a list to hold the set.)
    """
    assert(len(s) == len(f))  # each start time must match a finish time
    n = len(s)  # could be len f as well!
    a = []
    k = 0
    for m in range(1, n):
        if s[m] >= f[k]:
            a.append(m)
            k = m
    return a