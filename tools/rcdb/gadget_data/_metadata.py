"""Account for retained Python catalog containers without encoding their values."""
import sys
from collections import OrderedDict


_ENTRY_BYTES = sys.getsizeof(OrderedDict(((None, None),))) - sys.getsizeof(OrderedDict())


def retained_metadata_bytes(key, value, maximum):
    """Count owned builtin containers and stop once retention exceeds admission.

    Shared objects count once within an entry. Unknown objects are not cached.
    Allocator overhead, transient decoding and process RSS remain separate.
    """
    from itertools import chain
    pending, seen = [iter((key, value))], set()
    total = _ENTRY_BYTES + sys.getsizeof((None, 0))
    while pending and total <= maximum:
        try:
            item = next(pending[-1])
        except StopIteration:
            pending.pop()
            continue
        identity = id(item)
        if identity in seen:
            continue
        seen.add(identity)
        if type(item) not in (dict, tuple, list, str, bytes, int, bool, float, type(None)):
            return maximum + 1
        total += sys.getsizeof(item)
        if type(item) is dict:
            pending.append(iter(chain(item.keys(), item.values())))
        elif type(item) in (tuple, list):
            pending.append(iter(item))
    return total
