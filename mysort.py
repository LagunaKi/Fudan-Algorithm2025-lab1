def my_sort(lst, key=lambda x: x):
    _quicksort(lst, 0, len(lst)-1, key)

def _quicksort(lst, low, high, key):
    if low < high:
        pivot_idx = _partition(lst, low, high, key)
        _quicksort(lst, low, pivot_idx-1, key)
        _quicksort(lst, pivot_idx+1, high, key)
def _partition(lst, low, high, key):
    pivot = lst[high]
    i = low - 1
    for j in range(low, high):
        if key(lst[j]) <= key(pivot):
            i += 1
            lst[i], lst[j] = lst[j], lst[i]
    lst[i+1], lst[high] = lst[high], lst[i+1]
    return i + 1
